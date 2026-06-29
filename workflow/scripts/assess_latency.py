#!/usr/bin/env python3
"""
assess_latency.py
=================
Assess EBV latency type from gene expression patterns.

Reads featureCounts output for both EBV-1 and EBV-2 alignments, computes
TPM for each EBV gene, and classifies the latency type (I, II, III, or
lytic-active) based on which latency genes are expressed.

Latency classification logic:
  Latency III : EBNA2 + EBNA3A/3B/3C + LMP1 + LMP2A all expressed (>10 TPM)
  Latency II  : LMP1 and/or LMP2A expressed, EBNA2 absent (<1 TPM)
  Latency I   : only EBNA1 + EBER1/2 expressed
  Lytic       : lytic genes > 50% of total EBV expression

Also determines which typing strategy is viable:
  - EBNA2 reads > 100  → EBNA2-based typing (highest confidence)
  - EBNA3 reads > 100  → EBNA3-based typing
  - Both < 100         → genome-wide SNP concordance (lower confidence)
  - Total EBV < 1000   → insufficient data

Snakemake passes inputs/outputs/params as snakecall globals.
"""

import json
import os
import re
import sys
from collections import defaultdict

import pandas as pd

# Snakemake globals (injected by Snakemake's script: directive)
# snakemake.input.ebv1_counts, snakemake.input.ebv2_counts, etc.


def read_featurecounts(path):
    """Read featureCounts output into a DataFrame.

    featureCounts format: comment lines starting with '#', then a header
    row, then one row per gene.  Columns:
      Geneid, Chr, Start, End, Strand, Length, <sample_bam>
    """
    # Read, skipping comment lines
    df = pd.read_csv(path, sep='\t', comment='#')
    # The last column is the count
    count_col = df.columns[-1]
    df = df.rename(columns={count_col: 'count', 'Geneid': 'gene',
                            'Length': 'length'})
    # Keep only relevant columns
    cols = ['gene', 'length', 'count']
    for c in cols:
        if c not in df.columns:
            # featureCounts may use different column names
            pass
    df = df[['gene', 'length', 'count']].copy()
    return df


def compute_tpm(df):
    """Add a TPM column to the featureCounts DataFrame."""
    # TPM = (count / length) * (1e6 / sum(count/length))
    df = df.copy()
    df['rpk'] = df['count'] / df['length'].replace(0, 1)
    total_rpk = df['rpk'].sum()
    if total_rpk == 0:
        df['tpm'] = 0.0
    else:
        df['tpm'] = df['rpk'] / total_rpk * 1e6
    return df


def load_gene_categories(path):
    """Load gene → category mapping."""
    cat = {}
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                cat[parts[0]] = parts[1]
    return cat


def classify_latency(df1, df2, gene_categories, expressed_tpm, absent_tpm):
    """Classify latency type from expression patterns.

    Uses the alignment (EBV-1 or EBV-2) that has more total EBV reads,
    since that is the matched-type alignment.
    """
    # Pick the alignment with more total reads
    total1 = df1['count'].sum()
    total2 = df2['count'].sum()
    primary_df = df1 if total1 >= total2 else df2
    primary_type = 'EBV-1' if total1 >= total2 else 'EBV-2'
    secondary_df = df2 if total1 >= total2 else df1

    # Build gene → TPM dict
    tpm = dict(zip(primary_df['gene'], primary_df['tpm']))
    counts = dict(zip(primary_df['gene'], primary_df['count']))

    # Check expression of key latency genes
    # Normalise gene names (featureCounts uses GTF gene_name)
    def get_tpm(gene):
        return tpm.get(gene, 0.0)

    def get_count(gene):
        return counts.get(gene, 0)

    def is_expressed(gene):
        return get_tpm(gene) > expressed_tpm

    def is_absent(gene):
        return get_tpm(gene) < absent_tpm

    # Latency markers
    ebna2_expr = is_expressed('EBNA2')
    ebna3a_expr = is_expressed('EBNA3A')
    ebna3b_expr = is_expressed('EBNA3B') or is_expressed('EBNA3BC')
    ebna3c_expr = is_expressed('EBNA3C') or is_expressed('EBNA3BC')
    lmp1_expr = is_expressed('LMP1')
    lmp2a_expr = is_expressed('LMP2A')
    ebna1_expr = is_expressed('EBNA1')
    eber1_expr = is_expressed('EBER1')
    eber2_expr = is_expressed('EBER2')

    # Lytic gene expression
    lytic_total = 0.0
    latent_total = 0.0
    for _, row in primary_df.iterrows():
        cat = gene_categories.get(row['gene'], 'other')
        if cat == 'lytic':
            lytic_total += row['tpm']
        elif cat == 'latent':
            latent_total += row['tpm']
    total_expr = lytic_total + latent_total
    lytic_fraction = lytic_total / total_expr if total_expr > 0 else 0.0

    # Classify
    if ebna2_expr and (ebna3a_expr or ebna3b_expr or ebna3c_expr) and lmp1_expr:
        latency_type = 'Latency III'
    elif lmp1_expr or lmp2a_expr:
        if is_absent('EBNA2'):
            latency_type = 'Latency II'
        else:
            latency_type = 'Latency II/III (ambiguous)'
    elif ebna1_expr and is_absent('EBNA2') and is_absent('LMP1'):
        latency_type = 'Latency I'
    else:
        latency_type = 'Unclassified'

    # Override if lytic genes dominate
    if lytic_fraction > 0.5:
        latency_type = latency_type + ' (lytic-active)'

    # Determine typing viability
    ebna2_reads = get_count('EBNA2')
    ebna3_reads = get_count('EBNA3A') + get_count('EBNA3B') + \
                  get_count('EBNA3BC') + get_count('EBNA3C')
    total_ebv_reads = int(primary_df['count'].sum())

    if total_ebv_reads < 1000:
        typing_strategy = 'insufficient_data'
        typing_note = (f'Only {total_ebv_reads} EBV reads detected. '
                       'Recommend DNA-based PCR for reliable typing.')
    elif ebna2_reads > 100:
        typing_strategy = 'EBNA2-based'
        typing_note = (f'EBNA2 has {ebna2_reads} reads — sufficient for '
                       'high-confidence EBNA2-based typing.')
    elif ebna3_reads > 100:
        typing_strategy = 'EBNA3-based'
        typing_note = (f'EBNA2 has only {ebna2_reads} reads but EBNA3 has '
                       f'{ebna3_reads} — EBNA3-based typing recommended.')
    else:
        typing_strategy = 'genome-wide SNP'
        typing_note = (f'EBNA2 ({ebna2_reads} reads) and EBNA3 ({ebna3_reads} '
                       'reads) both have low coverage — falling back to '
                       'genome-wide SNP concordance (lower confidence).')

    # Per-gene expression table
    gene_table = []
    for _, row in primary_df.iterrows():
        gene_table.append({
            'gene': row['gene'],
            'count': int(row['count']),
            'tpm': round(row['tpm'], 2),
            'category': gene_categories.get(row['gene'], 'other'),
            'expressed': row['tpm'] > expressed_tpm,
        })
    gene_table.sort(key=lambda x: -x['tpm'])

    result = {
        'sample': snakemake.wildcards.sample,
        'latency_type': latency_type,
        'primary_alignment': primary_type,
        'total_ebv_reads': total_ebv_reads,
        'ebv1_total_reads': int(total1),
        'ebv2_total_reads': int(total2),
        'lytic_fraction': round(lytic_fraction, 4),
        'typing_strategy': typing_strategy,
        'typing_note': typing_note,
        'marker_expression': {
            'EBNA2': {'tpm': round(get_tpm('EBNA2'), 2),
                      'reads': int(get_count('EBNA2')),
                      'expressed': ebna2_expr},
            'EBNA3A': {'tpm': round(get_tpm('EBNA3A'), 2),
                       'reads': int(get_count('EBNA3A')),
                       'expressed': ebna3a_expr},
            'EBNA3B': {'tpm': round(get_tpm('EBNA3B'), 2),
                       'reads': int(get_count('EBNA3B')),
                       'expressed': ebna3b_expr},
            'EBNA3C': {'tpm': round(get_tpm('EBNA3C'), 2),
                       'reads': int(get_count('EBNA3C')),
                       'expressed': ebna3c_expr},
            'LMP1': {'tpm': round(get_tpm('LMP1'), 2),
                     'reads': int(get_count('LMP1')),
                     'expressed': lmp1_expr},
            'LMP2A': {'tpm': round(get_tpm('LMP2A'), 2),
                      'reads': int(get_count('LMP2A')),
                      'expressed': lmp2a_expr},
            'EBNA1': {'tpm': round(get_tpm('EBNA1'), 2),
                      'reads': int(get_count('EBNA1')),
                      'expressed': ebna1_expr},
            'EBER1': {'tpm': round(get_tpm('EBER1'), 2),
                      'reads': int(get_count('EBER1')),
                      'expressed': eber1_expr},
            'EBER2': {'tpm': round(get_tpm('EBER2'), 2),
                      'reads': int(get_count('EBER2')),
                      'expressed': eber2_expr},
        },
        'gene_table': gene_table,
    }
    return result


def main():
    # Read inputs
    df1 = read_featurecounts(snakemake.input.ebv1_counts)
    df2 = read_featurecounts(snakemake.input.ebv2_counts)

    # Compute TPM
    df1 = compute_tpm(df1)
    df2 = compute_tpm(df2)

    # Load gene categories
    gene_categories = load_gene_categories(snakemake.input.gene_categories)

    # Classify
    expressed_tpm = snakemake.params.expressed_tpm
    absent_tpm = snakemake.params.absent_tpm
    result = classify_latency(df1, df2, gene_categories, expressed_tpm, absent_tpm)

    # Write JSON
    with open(snakemake.output.latency_json, 'w') as out:
        json.dump(result, out, indent=2)

    # Write TSV (per-gene expression)
    with open(snakemake.output.latency_tsv, 'w') as out:
        out.write('gene\tcount\ttpm\tcategory\texpressed\n')
        for g in result['gene_table']:
            out.write(f"{g['gene']}\t{g['count']}\t{g['tpm']}\t{g['category']}\t{g['expressed']}\n")

    print(f"[assess_latency] Sample: {result['sample']}")
    print(f"  Latency type: {result['latency_type']}")
    print(f"  Primary alignment: {result['primary_alignment']}")
    print(f"  Total EBV reads: {result['total_ebv_reads']}")
    print(f"  Typing strategy: {result['typing_strategy']}")
    print(f"  Typing note: {result['typing_note']}")


if __name__ == '__main__' or 'snakemake' in globals():
    if 'snakemake' in globals():
        main()
