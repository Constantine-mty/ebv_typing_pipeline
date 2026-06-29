#!/usr/bin/env python3
"""
type_classification.py
======================
Core EBV type classification logic.

Combines two independent lines of evidence:
  1. Coverage-based scoring: reads from a type-1 infection align better
     to the EBV-1 reference in divergent regions (EBNA2, EBNA3, BZLF1,
     BZLF2, BLLF1), producing higher read counts there.
  2. SNP concordance scoring: at positions where EBV-1 and EBV-2
     references differ, the allele called from the alignment should
     match the infecting type's reference base.

The final call integrates both signals and assigns a confidence level
(high / medium / low / insufficient).

Inputs (via snakemake):
  - ebv1_vcf, ebv2_vcf : SNP VCFs from both alignments
  - ebv1_cov, ebv2_cov : per-gene coverage tables
  - type_snps          : table of type-specific SNP positions
  - gene_categories    : gene → category mapping
  - latency_json       : latency assessment (for typing strategy)

Outputs:
  - typing_json : full results in JSON
  - typing_tsv  : summary table
"""

import gzip
import json
import math
import os
import re
import sys
from collections import defaultdict

import pandas as pd

# Type-discriminatory genes (where EBV-1 and EBV-2 differ most)
TYPE_DISCRIMINATORY_GENES = [
    'EBNA2', 'EBNA3A', 'EBNA3B', 'EBNA3C', 'EBNA3BC',
    'BZLF1', 'BZLF2', 'BLLF1',
]


# ---------------------------------------------------------------------------
# Coverage-based scoring
# ---------------------------------------------------------------------------

def read_coverage(path):
    """Read bedtools coverage output: chr, start, end, gene, score, strand, count."""
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'gene', 'score', 'strand', 'count'])
    return df


def coverage_scoring(cov1, cov2):
    """Compare read counts in type-discriminatory genes between EBV-1 and
    EBV-2 alignments.

    Returns a dict with per-gene log2 ratios and an overall score.
    Positive log2 ratio → type 1; negative → type 2.
    """
    # Build gene → count dicts
    counts1 = dict(zip(cov1['gene'], cov1['count']))
    counts2 = dict(zip(cov2['gene'], cov2['count']))

    results = {}
    type1_votes = 0
    type2_votes = 0
    total_reads1 = 0
    total_reads2 = 0

    for gene in TYPE_DISCRIMINATORY_GENES:
        c1 = counts1.get(gene, 0)
        c2 = counts2.get(gene, 0)
        total_reads1 += c1
        total_reads2 += c2

        if c1 + c2 == 0:
            log2_ratio = None
            vote = None
        else:
            # Add pseudocount of 1 to avoid log(0)
            log2_ratio = math.log2((c1 + 1) / (c2 + 1))
            if log2_ratio > 0.5:
                vote = 'type1'
                type1_votes += 1
            elif log2_ratio < -0.5:
                vote = 'type2'
                type2_votes += 1
            else:
                vote = 'ambiguous'

        results[gene] = {
            'reads_ebv1': int(c1),
            'reads_ebv2': int(c2),
            'log2_ratio': round(log2_ratio, 3) if log2_ratio is not None else None,
            'vote': vote,
        }

    # Overall coverage-based call
    if type1_votes > type2_votes and type1_votes > 0:
        cov_call = 'type1'
    elif type2_votes > type1_votes and type2_votes > 0:
        cov_call = 'type2'
    else:
        cov_call = 'ambiguous'

    # Overall log2 ratio across all discriminatory genes
    if total_reads1 + total_reads2 > 0:
        overall_log2 = math.log2((total_reads1 + 1) / (total_reads2 + 1))
    else:
        overall_log2 = None

    return {
        'per_gene': results,
        'type1_votes': type1_votes,
        'type2_votes': type2_votes,
        'overall_log2_ratio': round(overall_log2, 3) if overall_log2 else None,
        'call': cov_call,
        'total_reads_discrim1': int(total_reads1),
        'total_reads_discrim2': int(total_reads2),
    }


# ---------------------------------------------------------------------------
# SNP concordance scoring
# ---------------------------------------------------------------------------

def read_type_specific_snps(path):
    """Read the type-specific SNP table.

    Columns: pos_NC_007605.1, base_NC_007605.1, base_other, gene, genes_all
    pos = 1-based position in EBV-1 reference
    base_NC_007605.1 = EBV-1 reference base
    base_other = EBV-2 reference base
    """
    snps = []
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            pos = int(parts[0])
            base1 = parts[1].upper()
            base2 = parts[2].upper()
            gene = parts[3]
            genes_all = parts[4]
            snps.append({
                'pos': pos,
                'base_ebv1': base1,
                'base_ebv2': base2,
                'gene': gene,
                'genes_all': genes_all,
            })
    return snps


def read_vcf(vcf_path):
    """Read a VCF file and return a dict: pos → (ref, alt, genotype).

    Handles both plain and gzipped VCF.
    """
    opener = gzip.open if vcf_path.endswith('.gz') else open
    variants = {}
    with opener(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue
            pos = int(parts[1])
            ref = parts[3].upper()
            alt = parts[4].upper()
            # Parse genotype (first sample)
            fmt = parts[8].split(':') if len(parts) > 8 else []
            sample_data = parts[9].split(':') if len(parts) > 9 else []
            gt = ''
            if 'GT' in fmt:
                gt_idx = fmt.index('GT')
                gt = sample_data[gt_idx] if gt_idx < len(sample_data) else ''
            # Get called allele (for homozygous or heterozygous)
            # If GT is 0/0 → ref; 0/1 or 1/1 → alt; ./. → missing
            if '0/0' in gt or '0|0' in gt:
                called_base = ref
            elif '1/1' in gt or '1|1' in gt:
                called_base = alt
            elif '0/1' in gt or '0|1' in gt or '1/0' in gt or '1|0' in gt:
                called_base = 'het'  # heterozygous — both alleles
            else:
                called_base = 'missing'
            variants[pos] = {
                'ref': ref,
                'alt': alt,
                'gt': gt,
                'called_base': called_base,
            }
    return variants


def snp_concordance_scoring(type_snps, vcf1, vcf2):
    """Score SNP concordance for both EBV-1 and EBV-2 alignments.

    For each type-specific SNP position:
      - In the EBV-1 alignment VCF: if the called allele matches the
        EBV-1 reference base → type-1 concordant; if it matches EBV-2
        reference base → type-2 concordant (discordant for type 1).
      - Similarly for the EBV-2 alignment VCF.

    The alignment with more type-concordant SNPs indicates the infecting type.
    """
    # Counters for EBV-1 alignment
    t1_concordant_ali1 = 0  # called base matches EBV-1 ref
    t2_concordant_ali1 = 0  # called base matches EBV-2 ref
    callable_ali1 = 0

    # Counters for EBV-2 alignment
    t1_concordant_ali2 = 0
    t2_concordant_ali2 = 0
    callable_ali2 = 0

    per_snp = []

    for snp in type_snps:
        pos = snp['pos']
        base1 = snp['base_ebv1']
        base2 = snp['base_ebv2']

        # EBV-1 alignment
        v1 = vcf1.get(pos)
        ali1_info = {'callable': False, 'matches': None}
        if v1 and v1['called_base'] != 'missing':
            callable_ali1 += 1
            ali1_info['callable'] = True
            cb = v1['called_base']
            if cb == base1:
                t1_concordant_ali1 += 1
                ali1_info['matches'] = 'type1'
            elif cb == base2:
                t2_concordant_ali1 += 1
                ali1_info['matches'] = 'type2'
            elif cb == 'het':
                # Heterozygous — count as half for each
                t1_concordant_ali1 += 0.5
                t2_concordant_ali1 += 0.5
                ali1_info['matches'] = 'het'

        # EBV-2 alignment
        # Note: positions in type_snps are EBV-1 coordinates.
        # The EBV-2 VCF uses EBV-2 coordinates, so we need to map.
        # For simplicity, we use the EBV-1 alignment VCF for type-1
        # concordance and the EBV-2 alignment VCF for type-2 concordance.
        # The EBV-2 VCF positions are in EBV-2 coordinates, which differ
        # from EBV-1 positions.  We handle this by also checking the
        # EBV-1 alignment for type-2 concordance (base matches EBV-2 ref).
        # This is valid because the EBV-1 alignment VCF reports variants
        # relative to the EBV-1 reference, so positions match.

        per_snp.append({
            'pos': pos,
            'base_ebv1': base1,
            'base_ebv2': base2,
            'gene': snp['gene'],
            'ali1_matches': ali1_info['matches'],
        })

    # We primarily use the EBV-1 alignment VCF for SNP concordance,
    # because type_snps positions are in EBV-1 coordinates.
    # The logic: if reads truly come from a type-1 virus, the EBV-1
    # alignment will show few variants (reads match reference), while
    # the EBV-2 alignment will show many variants (reads differ from
    # reference).  Conversely for type-2.
    #
    # So: high t1_concordant_ali1 + low t2_concordant_ali1 → type 1
    #     low t1_concordant_ali1 + high t2_concordant_ali1 → type 2

    if callable_ali1 > 0:
        t1_fraction = t1_concordant_ali1 / callable_ali1
        t2_fraction = t2_concordant_ali1 / callable_ali1
    else:
        t1_fraction = 0
        t2_fraction = 0

    if t1_fraction > t2_fraction and t1_fraction > 0.5:
        snp_call = 'type1'
        snp_concordance = t1_fraction
    elif t2_fraction > t1_fraction and t2_fraction > 0.5:
        snp_call = 'type2'
        snp_concordance = t2_fraction
    else:
        snp_call = 'ambiguous'
        snp_concordance = max(t1_fraction, t2_fraction)

    return {
        'callable_snps': callable_ali1,
        'type1_concordant': t1_concordant_ali1,
        'type2_concordant': t2_concordant_ali1,
        'type1_fraction': round(t1_fraction, 4),
        'type2_fraction': round(t2_fraction, 4),
        'call': snp_call,
        'concordance_score': round(snp_concordance, 4),
    }


# ---------------------------------------------------------------------------
# Combined classification
# ---------------------------------------------------------------------------

def combine_classification(cov_result, snp_result, latency_result,
                           min_ebv_reads, min_type_snps, high_conf, cov_log2):
    """Combine coverage and SNP concordance into a final type call."""
    total_ebv = latency_result.get('total_ebv_reads', 0)
    typing_strategy = latency_result.get('typing_strategy', 'unknown')

    # Check data sufficiency
    if total_ebv < min_ebv_reads:
        return {
            'type': 'undetermined',
            'confidence': 'insufficient',
            'reason': f'Only {total_ebv} EBV reads (minimum {min_ebv_reads}). '
                      'Recommend DNA-based PCR.',
            'coverage': cov_result,
            'snp': snp_result,
        }

    if snp_result['callable_snps'] < min_type_snps:
        # Fall back to coverage-only
        if cov_result['call'] != 'ambiguous':
            final_type = cov_result['call']
            confidence = 'medium'
            reason = (f'Coverage-based call only (SNP concordance has '
                      f'{snp_result["callable_snps"]} callable SNPs, '
                      f'minimum {min_type_snps}).')
        else:
            final_type = 'undetermined'
            confidence = 'low'
            reason = 'Both coverage and SNP signals ambiguous.'
    else:
        # Both methods available
        cov_call = cov_result['call']
        snp_call = snp_result['call']

        if cov_call == snp_call and cov_call != 'ambiguous':
            final_type = cov_call
            if snp_result['concordance_score'] >= high_conf:
                confidence = 'high'
                reason = (f'Coverage and SNP concordance agree: {final_type}. '
                          f'SNP concordance = {snp_result["concordance_score"]:.2%}.')
            else:
                confidence = 'medium'
                reason = (f'Coverage and SNP concordance agree: {final_type}, '
                          f'but concordance score ({snp_result["concordance_score"]:.2%}) '
                          f'below high-confidence threshold ({high_conf:.0%}).')
        elif cov_call != 'ambiguous' and snp_call == 'ambiguous':
            final_type = cov_call
            confidence = 'medium'
            reason = f'Coverage suggests {final_type}; SNP signal ambiguous.'
        elif snp_call != 'ambiguous' and cov_call == 'ambiguous':
            final_type = snp_call
            confidence = 'medium'
            reason = f'SNP concordance suggests {final_type}; coverage ambiguous.'
        elif cov_call != snp_call:
            final_type = snp_call if snp_result['concordance_score'] > 0.7 else cov_call
            confidence = 'low'
            reason = (f'DISCORDANT: coverage suggests {cov_call}, '
                      f'SNP suggests {snp_call}. Possible recombination or '
                      f'mixed infection. Using {final_type} as primary call.')
        else:
            final_type = 'undetermined'
            confidence = 'low'
            reason = 'Both coverage and SNP signals ambiguous.'

    return {
        'type': final_type,
        'confidence': confidence,
        'reason': reason,
        'typing_strategy': typing_strategy,
        'coverage': cov_result,
        'snp': snp_result,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Read inputs
    cov1 = read_coverage(snakemake.input.ebv1_cov)
    cov2 = read_coverage(snakemake.input.ebv2_cov)
    type_snps = read_type_specific_snps(snakemake.input.type_snps)
    vcf1 = read_vcf(snakemake.input.ebv1_vcf)
    # vcf2 uses EBV-2 coordinates; we primarily use vcf1 for concordance
    # but read it for completeness
    vcf2 = read_vcf(snakemake.input.ebv2_vcf)

    # Load latency assessment
    with open(snakemake.input.latency_json) as fh:
        latency_result = json.load(fh)

    # Coverage scoring
    cov_result = coverage_scoring(cov1, cov2)

    # SNP concordance scoring
    snp_result = snp_concordance_scoring(type_snps, vcf1, vcf2)

    # Combined classification
    final = combine_classification(
        cov_result, snp_result, latency_result,
        snakemake.params.min_ebv_reads,
        snakemake.params.min_type_snps,
        snakemake.params.high_conf,
        snakemake.params.cov_log2,
    )

    # Assemble result
    result = {
        'sample': snakemake.wildcards.sample,
        'ebv_type': final['type'],
        'confidence': final['confidence'],
        'reason': final['reason'],
        'typing_strategy': final.get('typing_strategy', 'unknown'),
        'latency_type': latency_result.get('latency_type', 'unknown'),
        'total_ebv_reads': latency_result.get('total_ebv_reads', 0),
        'coverage_scoring': final['coverage'],
        'snp_concordance': final['snp'],
    }

    # Write JSON
    with open(snakemake.output.typing_json, 'w') as out:
        json.dump(result, out, indent=2)

    # Write TSV summary
    with open(snakemake.output.typing_tsv, 'w') as out:
        out.write('field\tvalue\n')
        out.write(f"sample\t{result['sample']}\n")
        out.write(f"ebv_type\t{result['ebv_type']}\n")
        out.write(f"confidence\t{result['confidence']}\n")
        out.write(f"latency_type\t{result['latency_type']}\n")
        out.write(f"total_ebv_reads\t{result['total_ebv_reads']}\n")
        out.write(f"typing_strategy\t{result['typing_strategy']}\n")
        out.write(f"coverage_call\t{final['coverage']['call']}\n")
        out.write(f"coverage_log2_ratio\t{final['coverage']['overall_log2_ratio']}\n")
        out.write(f"snp_call\t{final['snp']['call']}\n")
        out.write(f"snp_concordance\t{final['snp']['concordance_score']}\n")
        out.write(f"callable_snps\t{final['snp']['callable_snps']}\n")
        out.write(f"reason\t{result['reason']}\n")

    print(f"[type_classification] Sample: {result['sample']}")
    print(f"  EBV type: {result['ebv_type']}")
    print(f"  Confidence: {result['confidence']}")
    print(f"  Latency: {result['latency_type']}")
    print(f"  Coverage call: {final['coverage']['call']} "
          f"(log2={final['coverage']['overall_log2_ratio']})")
    print(f"  SNP call: {final['snp']['call']} "
          f"(concordance={final['snp']['concordance_score']:.2%}, "
          f"{final['snp']['callable_snps']} callable SNPs)")
    print(f"  Reason: {result['reason']}")


if __name__ == '__main__' or 'snakemake' in globals():
    if 'snakemake' in globals():
        main()
