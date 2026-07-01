#!/usr/bin/env python3
"""
prepare_references.py
=====================
Prepare EBV reference resources for the typing pipeline:
  1. Convert NCBI GFF3 gene annotations to GTF (for featureCounts).
  2. Align EBV-1 vs EBV-2 reference genomes with minimap2 (via mappy),
     extract all SNP positions from the pairwise alignment.
  3. Annotate type-specific SNPs by gene region (EBNA2, EBNA3A/3B/3C,
     BZLF1, BZLF2, BLLF1, etc.).
  4. Compile EBV gene categories table (latent vs lytic), with
     normalised gene names so the dictionary keys match the GFF3 names.

Usage:
  python prepare_references.py \
      --ebv1-fa resources/ebv1_reference.fa \
      --ebv2-fa resources/ebv2_reference.fa \
      --ebv1-gff3 resources/ebv1_reference.gff3 \
      --ebv2-gff3 resources/ebv2_reference.gff3 \
      --outdir resources/
"""

import argparse
import os
import re
import sys
from collections import defaultdict

import mappy  # minimap2 Python bindings


# ---------------------------------------------------------------------------
# Gene-name normalisation
# ---------------------------------------------------------------------------
# NCBI GFF3 uses hyphenated names (EBNA-2, EBNA-3A, LMP-2A, EBER-1).
# We normalise to the canonical literature names (EBNA2, EBNA3A, LMP2A, EBER1)
# so the gene-category dictionary and downstream code can use a single
# naming convention.

GENE_NAME_NORMALISATION = {
    'EBNA-1': 'EBNA1',
    'EBNA-2': 'EBNA2',
    'EBNA-3A': 'EBNA3A',
    'EBNA-3B': 'EBNA3B',
    'EBNA-3C': 'EBNA3C',
    'EBNA-3B/EBNA-3C': 'EBNA3BC',  # merged annotation in NC_007605
    'EBNA-LP': 'EBNALP',
    'LMP-1': 'LMP1',
    'LMP-2A': 'LMP2A',
    'LMP-2B': 'LMP2B',
    'EBER-1': 'EBER1',
    'EBER-2': 'EBER2',
    'BART': 'BART',
    'RPMS1': 'RPMS1',
}


def normalise_gene_name(name):
    return GENE_NAME_NORMALISATION.get(name, name)


# ---------------------------------------------------------------------------
# 1. GFF3 -> GTF conversion (gene-level, for featureCounts)
# ---------------------------------------------------------------------------

def gff3_to_gtf(gff3_path, gtf_path, seqid):
    """Convert NCBI GFF3 to a minimal GTF with 'gene' features.

    Gene names are normalised to canonical literature names.
    For EBV-2 (AG876), the GFF3 gene features only have locus_tag names
    (e.g. HHV4tp2_gp06), so we scan CDS/transcript features for
    'product=' attributes (e.g. product=EBNA-2) and map them back to
    the parent gene via the ID= / Parent= linkage.
    """
    # First pass: build gene_id → locus_tag and CDS product → parent gene
    gene_id_to_locus = {}   # gene-ID → locus_tag
    parent_to_product = {}  # gene-ID → product name (from CDS/transcript)

    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            ftype = fields[2]
            attrs = fields[8]

            if ftype == 'gene':
                # Extract ID and locus_tag
                gid = None
                locus = None
                for kv in attrs.split(';'):
                    if kv.startswith('ID='):
                        gid = kv.split('=')[1]
                    elif kv.startswith('locus_tag='):
                        locus = kv.split('=')[1]
                    elif kv.startswith('Name='):
                        locus = kv.split('=')[1]
                if gid:
                    gene_id_to_locus[gid] = locus or gid
            elif ftype in ('CDS', 'transcript', 'mRNA'):
                # Extract Parent and product
                parent = None
                product = None
                for kv in attrs.split(';'):
                    if kv.startswith('Parent='):
                        parent = kv.split('=')[1]
                    elif kv.startswith('product='):
                        product = kv.split('=')[1]
                if parent and product:
                    # parent may be a gene-ID; store the product
                    # (keep the first product found, or overwrite —
                    #  for EBV-2 the CDS product is the gene name)
                    if parent not in parent_to_product:
                        parent_to_product[parent] = product

    # Second pass: emit GTF gene lines
    gene_re = re.compile(r'\tgene\t')
    out_lines = []
    with open(gff3_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            if not gene_re.search(line):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9:
                continue
            attrs = fields[8]

            # Extract gene ID
            gid = None
            for kv in attrs.split(';'):
                if kv.startswith('ID='):
                    gid = kv.split('=')[1]
                    break

            # Try to get a meaningful name:
            # 1. product from CDS/transcript (EBV-2 case)
            # 2. Name= or gene= attribute (EBV-1 case)
            # 3. locus_tag fallback
            name = None
            if gid and gid in parent_to_product:
                name = parent_to_product[gid]
            if name is None:
                for kv in attrs.split(';'):
                    if kv.startswith('Name='):
                        name = kv.split('=')[1]
                        break
                    if kv.startswith('gene='):
                        name = kv.split('=')[1]
                        break
            if name is None:
                for kv in attrs.split(';'):
                    if kv.startswith('locus_tag='):
                        name = kv.split('=')[1]
                        break
            if name is None:
                continue

            name = normalise_gene_name(name)
            start, end, strand = fields[3], fields[4], fields[6]
            gtf_attrs = f'gene_id "{name}"; gene_name "{name}";'
            out_lines.append(
                f'{seqid}\tRefSeq\tgene\t{start}\t{end}\t.\t{strand}\t.\t{gtf_attrs}'
            )
    with open(gtf_path, 'w') as out:
        out.write('\n'.join(out_lines) + '\n')
    print(f"[gff3_to_gtf] {len(out_lines)} genes written to {gtf_path}")
    return out_lines


# ---------------------------------------------------------------------------
# 2. Whole-genome alignment of EBV-1 vs EBV-2 with minimap2 and SNP extraction
# ---------------------------------------------------------------------------

def read_fasta(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                continue
            seq.append(line.strip())
    return ''.join(seq).upper()


def find_snp_positions_minimap(seq1_path, seq2_path):
    """Align EBV-1 (query) to EBV-2 (target) with minimap2 and extract SNPs.

    Uses preset 'asm10' which tolerates up to ~10% sequence divergence —
    necessary because the EBNA2 locus differs by ~30% between EBV-1 and
    EBV-2, and asm5 (>5% threshold) breaks alignment in that region.
    asm10 recovers ~250 SNPs in EBNA2 and ~770 in the EBNA3 cluster,
    consistent with the known intertype divergence.

    Returns list of (pos1_1based, pos2_1based, base1, base2).
    """
    # Index the target (EBV-2)
    seq2 = read_fasta(seq2_path)

    idx = mappy.Aligner(seq2_path, preset='asm10')
    if idx is None:
        raise RuntimeError("Failed to build minimap2 index for EBV-2")

    seq1 = read_fasta(seq1_path)
    snps = []

    for hit in idx.map(seq1):
        # hit.r_st, hit.r_en : target (EBV-2) coordinates
        # hit.q_st, hit.q_en : query (EBV-1) coordinates
        # hit.cigar : CIGAR list of (length, op)
        if hit.is_primary:
            # Walk the CIGAR to extract mismatches
            qpos = hit.q_st  # 0-based
            rpos = hit.r_st  # 0-based
            for length, op in hit.cigar:
                if op == 0 or op == 7 or op == 8:  # M, =, X
                    for _ in range(length):
                        if qpos < len(seq1) and rpos < len(seq2):
                            b1 = seq1[qpos]
                            b2 = seq2[rpos]
                            if b1 != b2:
                                snps.append((qpos + 1, rpos + 1, b1, b2))
                        qpos += 1
                        rpos += 1
                elif op == 1:  # insertion in query (gap in target)
                    qpos += length
                elif op == 2:  # deletion in query (gap in target)
                    rpos += length
                elif op == 3:  # N (intron skip) — shouldn't happen for DNA virus
                    rpos += length
                # ignore other ops
    return snps


def write_snp_table(snps, out_path, seqid1, seqid2):
    """Write SNP table: pos1, pos2, base1, base2."""
    with open(out_path, 'w') as out:
        out.write(f'pos_{seqid1}\tpos_{seqid2}\tbase_{seqid1}\tbase_{seqid2}\n')
        for pos1, pos2, b1, b2 in snps:
            out.write(f'{pos1}\t{pos2}\t{b1}\t{b2}\n')
    print(f"[find_snp_positions] {len(snps)} SNPs written to {out_path}")


# ---------------------------------------------------------------------------
# 3. Annotate type-specific SNPs by gene region
# ---------------------------------------------------------------------------

def load_gene_regions(gtf_path):
    """Return list of (gene_name, start, end, strand)."""
    regions = []
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 9 or fields[2] != 'gene':
                continue
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            m = re.search(r'gene_name "([^"]+)"', fields[8])
            name = m.group(1) if m else 'unknown'
            regions.append((name, start, end, strand))
    return regions


def annotate_snps_by_gene(snps, gene_regions, out_path, seqid):
    """Annotate each SNP with all genes it falls in.

    EBV genes are heavily nested (EBNA2 spans 11305-37739 but contains
    BWRF1 repeats and EBNALP).  For typing we need to know whether a SNP
    falls in the EBNA2 locus regardless of nested annotations, so we emit
    a comma-separated gene list and also a separate column flagging whether
    the SNP is in any type-discriminatory gene.
    """
    # Build interval lookup: for each position, find all overlapping genes
    with open(out_path, 'w') as out:
        out.write(f'pos_{seqid}\tbase_{seqid}\tbase_other\tgene\tgenes_all\n')
        for pos1, pos2, b1, b2 in snps:
            genes = []
            for gname, gstart, gend, gstrand in gene_regions:
                if gstart <= pos1 <= gend:
                    genes.append(gname)
            if not genes:
                primary_gene = 'intergenic'
                all_genes = 'intergenic'
            else:
                # Primary: prefer type-discriminatory gene, then latent gene,
                # then the first hit
                primary_gene = genes[0]
                for g in genes:
                    if g in TYPE_DISCRIMINATORY_GENES:
                        primary_gene = g
                        break
                all_genes = ','.join(genes)
            out.write(f'{pos1}\t{b1}\t{b2}\t{primary_gene}\t{all_genes}\n')
    print(f"[annotate_snps_by_gene] annotated SNPs written to {out_path}")


# ---------------------------------------------------------------------------
# 4. EBV gene categories (latent vs lytic)
# ---------------------------------------------------------------------------

# Canonical EBV gene classification based on the literature.
# Keys use normalised names (matching GFF3 after normalisation).
GENE_CATEGORIES = {
    # Latency-associated genes
    'EBNA1': 'latent', 'EBNA2': 'latent', 'EBNA3A': 'latent',
    'EBNA3B': 'latent', 'EBNA3C': 'latent', 'EBNA3BC': 'latent',
    'EBNALP': 'latent',
    'LMP1': 'latent', 'LMP2A': 'latent', 'LMP2B': 'latent',
    'EBER1': 'latent', 'EBER2': 'latent',
    'RPMS1': 'latent', 'A73': 'latent', 'BART': 'latent',
    # Lytic genes
    'BZLF1': 'lytic', 'BRLF1': 'lytic', 'BLLF1': 'lytic',
    'BGLF4': 'lytic', 'BGLF5': 'lytic', 'BALF5': 'lytic',
    'BALF2': 'lytic', 'BALF4': 'lytic', 'BMRF1': 'lytic',
    'BMRF2': 'lytic', 'BFRF1': 'lytic', 'BFRF2': 'lytic',
    'BFRF3': 'lytic', 'BPLF1': 'lytic', 'BOLF1': 'lytic',
    'BORF1': 'lytic', 'BORF2': 'lytic', 'BSLF1': 'lytic',
    'BSLF2': 'lytic', 'BMLF1': 'lytic', 'BKRF4': 'lytic',
    'BILF2': 'lytic', 'BILF1': 'lytic', 'BZLF2': 'lytic',
    'BNLF1': 'lytic', 'BHRF1': 'lytic', 'BARF1': 'lytic',
    'BCRF1': 'lytic', 'BDLF1': 'lytic', 'BDLF2': 'lytic',
    'BDLF3': 'lytic', 'BDLF4': 'lytic', 'BXLF1': 'lytic',
    'BXLF2': 'lytic', 'BKRF1': 'lytic', 'BKRF2': 'lytic',
    'BKRF3': 'lytic', 'BLRF1': 'lytic', 'BLRF2': 'lytic',
    'BBRF1': 'lytic', 'BBRF2': 'lytic', 'BBRF3': 'lytic',
    'BTRF1': 'lytic', 'BTRF2': 'lytic', 'BTRF3': 'lytic',
    'BVLF1': 'lytic', 'BNRF1': 'lytic', 'BFLF1': 'lytic',
    'BFLF2': 'lytic', 'BFRF1A': 'lytic', 'BWRF1': 'lytic',
    'BHLF1': 'lytic',
}

# Genes critical for latency-type determination
LATENCY_MARKER_GENES = {
    'EBNA2': 'Latency III marker',
    'EBNA3A': 'Latency III marker',
    'EBNA3B': 'Latency III marker',
    'EBNA3C': 'Latency III marker',
    'EBNA3BC': 'Latency III marker',
    'EBNALP': 'Latency III marker',
    'LMP1': 'Latency II/III marker',
    'LMP2A': 'Latency II/III marker',
    'LMP2B': 'Latency II/III marker',
    'EBNA1': 'Latency I/II/III (always expressed)',
    'EBER1': 'Latency I/II/III (always expressed)',
    'EBER2': 'Latency I/II/III (always expressed)',
}

# Type-discriminatory genes (where EBV-1 and EBV-2 differ most)
TYPE_DISCRIMINATORY_GENES = [
    'EBNA2', 'EBNA3A', 'EBNA3B', 'EBNA3C', 'EBNA3BC',
    'BZLF1', 'BZLF2', 'BLLF1',
]


def write_gene_categories(gtf_path, out_path):
    """Write gene categories table: gene, category, role."""
    regions = load_gene_regions(gtf_path)
    seen = set()
    with open(out_path, 'w') as out:
        out.write('gene\tcategory\trole\n')
        for gname, start, end, strand in regions:
            if gname in seen:
                continue
            seen.add(gname)
            category = GENE_CATEGORIES.get(gname, 'other')
            role = LATENCY_MARKER_GENES.get(gname, '')
            if gname in TYPE_DISCRIMINATORY_GENES:
                role = role + ('; ' if role else '') + 'type-discriminatory'
            out.write(f'{gname}\t{category}\t{role}\n')
    print(f"[write_gene_categories] gene categories written to {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Prepare EBV reference resources')
    parser.add_argument('--ebv1-fa', required=True)
    parser.add_argument('--ebv2-fa', required=True)
    parser.add_argument('--ebv1-gff3', required=True)
    parser.add_argument('--ebv2-gff3', required=True)
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    seqid1 = 'NC_007605.1'
    seqid2 = 'NC_009334.1'

    # 1. GFF3 -> GTF
    gtf1 = os.path.join(args.outdir, 'ebv1_genes.gtf')
    gtf2 = os.path.join(args.outdir, 'ebv2_genes.gtf')
    gff3_to_gtf(args.ebv1_gff3, gtf1, seqid1)
    gff3_to_gtf(args.ebv2_gff3, gtf2, seqid2)

    # 2. Extract SNPs between EBV-1 and EBV-2 using minimap2
    print("[main] Aligning EBV-1 vs EBV-2 with minimap2 (asm5 preset)...")
    snps = find_snp_positions_minimap(args.ebv1_fa, args.ebv2_fa)
    snp_table = os.path.join(args.outdir, 'ebv1_vs_ebv2_snps.tsv')
    write_snp_table(snps, snp_table, seqid1, seqid2)

    # 3. Annotate SNPs by gene region
    gene_regions1 = load_gene_regions(gtf1)
    type_snps_path = os.path.join(args.outdir, 'type_specific_snps.tsv')
    annotate_snps_by_gene(snps, gene_regions1, type_snps_path, seqid1)

    # 4. Gene categories
    cat_path = os.path.join(args.outdir, 'ebv_gene_categories.tsv')
    write_gene_categories(gtf1, cat_path)

    # Summary
    print("\n=== Summary ===")
    print(f"EBV-1 genes: {len(gene_regions1)}")
    print(f"Total SNPs between EBV-1 and EBV-2: {len(snps)}")
    # Count SNPs in type-discriminatory genes (using genes_all column)
    disc_count = 0
    disc_by_gene = defaultdict(int)
    with open(type_snps_path) as fh:
        next(fh)
        for line in fh:
            fields = line.strip().split('\t')
            all_genes = fields[4] if len(fields) > 4 else fields[3]
            for g in all_genes.split(','):
                if g in TYPE_DISCRIMINATORY_GENES:
                    disc_count += 1
                    disc_by_gene[g] += 1
    print(f"SNPs in type-discriminatory genes: {disc_count}")
    for g, c in sorted(disc_by_gene.items(), key=lambda x: -x[1]):
        print(f"  {g}: {c}")
    print(f"\nResources written to {args.outdir}")


if __name__ == '__main__':
    main()
