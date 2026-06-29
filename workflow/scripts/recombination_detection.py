#!/usr/bin/env python3
"""
recombination_detection.py
===========================
Detect intertypic recombination or mixed infection.

Strategy:
  1. Independently type the EBNA2 region (using only SNPs/coverage in EBNA2).
  2. Independently type the EBNA3 region (using only SNPs/coverage in EBNA3A/3B/3C).
  3. Check BZLF1 region for the Zp-V3 promoter variant (type-2 hallmark).
  4. If EBNA2 type != EBNA3 type → flag as intertypic recombinant.
  5. If both types detected in the same region → flag as mixed infection.

Inputs (via snakemake):
  - ebv1_vcf, ebv2_vcf : SNP VCFs
  - ebv1_cov, ebv2_cov : per-gene coverage tables
  - type_snps          : type-specific SNP table
  - typing_json        : overall typing result

Outputs:
  - recomb_json : recombination analysis in JSON
  - recomb_tsv  : summary table
"""

import gzip
import json
import math
import os
import sys
from collections import defaultdict

import pandas as pd


# Region definitions (gene groups for independent typing)
REGIONS = {
    'EBNA2': ['EBNA2'],
    'EBNA3': ['EBNA3A', 'EBNA3B', 'EBNA3C', 'EBNA3BC'],
    'BZLF1': ['BZLF1'],
    'BZLF2': ['BZLF2'],
    'BLLF1': ['BLLF1'],
}


def read_coverage(path):
    """Read bedtools coverage output."""
    df = pd.read_csv(path, sep='\t', header=None,
                     names=['chr', 'start', 'end', 'gene', 'score', 'strand', 'count'])
    return df


def read_type_specific_snps(path):
    """Read type-specific SNP table."""
    snps = []
    with open(path) as fh:
        next(fh)
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
    """Read VCF and return pos → variant dict."""
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
            fmt = parts[8].split(':') if len(parts) > 8 else []
            sample_data = parts[9].split(':') if len(parts) > 9 else []
            gt = ''
            if 'GT' in fmt:
                gt_idx = fmt.index('GT')
                gt = sample_data[gt_idx] if gt_idx < len(sample_data) else ''
            if '0/0' in gt or '0|0' in gt:
                called_base = ref
            elif '1/1' in gt or '1|1' in gt:
                called_base = alt
            elif '0/1' in gt or '0|1' in gt or '1/0' in gt or '1|0' in gt:
                called_base = 'het'
            else:
                called_base = 'missing'
            variants[pos] = {'ref': ref, 'alt': alt, 'called_base': called_base}
    return variants


def type_region(snps_in_region, vcf1):
    """Type a single region using SNP concordance.

    Returns: type ('type1', 'type2', 'ambiguous'), details dict.
    """
    t1_conc = 0
    t2_conc = 0
    callable_n = 0

    for snp in snps_in_region:
        v = vcf1.get(snp['pos'])
        if v and v['called_base'] != 'missing':
            callable_n += 1
            cb = v['called_base']
            if cb == snp['base_ebv1']:
                t1_conc += 1
            elif cb == snp['base_ebv2']:
                t2_conc += 1
            elif cb == 'het':
                t1_conc += 0.5
                t2_conc += 0.5

    if callable_n == 0:
        return 'unknown', {'callable_snps': 0, 't1_concordant': 0, 't2_concordant': 0}

    t1_frac = t1_conc / callable_n
    t2_frac = t2_conc / callable_n

    if t1_frac > t2_frac and t1_frac > 0.6:
        call = 'type1'
    elif t2_frac > t1_frac and t2_frac > 0.6:
        call = 'type2'
    else:
        call = 'ambiguous'

    return call, {
        'callable_snps': callable_n,
        't1_concordant': t1_conc,
        't2_concordant': t2_conc,
        't1_fraction': round(t1_frac, 4),
        't2_fraction': round(t2_frac, 4),
    }


def coverage_region(cov1, cov2, genes):
    """Compare coverage for a set of genes between EBV-1 and EBV-2 alignments."""
    c1 = sum(int(cov1[cov1['gene'] == g]['count'].iloc[0]) for g in genes
             if g in cov1['gene'].values)
    c2 = sum(int(cov2[cov2['gene'] == g]['count'].iloc[0]) for g in genes
             if g in cov2['gene'].values)
    if c1 + c2 == 0:
        return 'unknown', 0, 0, None
    log2_ratio = math.log2((c1 + 1) / (c2 + 1))
    if log2_ratio > 0.5:
        call = 'type1'
    elif log2_ratio < -0.5:
        call = 'type2'
    else:
        call = 'ambiguous'
    return call, c1, c2, round(log2_ratio, 3)


def main():
    # Read inputs
    cov1 = read_coverage(snakemake.input.ebv1_cov)
    cov2 = read_coverage(snakemake.input.ebv2_cov)
    type_snps = read_type_specific_snps(snakemake.input.type_snps)
    vcf1 = read_vcf(snakemake.input.ebv1_vcf)

    # Load overall typing result
    with open(snakemake.input.typing_json) as fh:
        overall_typing = json.load(fh)

    # Group SNPs by region (using genes_all for overlap-aware assignment)
    snps_by_region = defaultdict(list)
    for snp in type_snps:
        genes_all = snp['genes_all'].split(',')
        for region_name, region_genes in REGIONS.items():
            if any(g in genes_all for g in region_genes):
                snps_by_region[region_name].append(snp)

    # Type each region independently
    region_results = {}
    for region_name in REGIONS:
        snp_call, snp_details = type_region(snps_by_region.get(region_name, []), vcf1)
        cov_call, c1, c2, log2r = coverage_region(cov1, cov2, REGIONS[region_name])
        region_results[region_name] = {
            'snp_call': snp_call,
            'snp_details': snp_details,
            'coverage_call': cov_call,
            'reads_ebv1': c1,
            'reads_ebv2': c2,
            'coverage_log2_ratio': log2r,
            'combined_call': snp_call if snp_call != 'unknown' else cov_call,
        }

    # Recombination detection
    ebna2_type = region_results.get('EBNA2', {}).get('combined_call', 'unknown')
    ebna3_type = region_results.get('EBNA3', {}).get('combined_call', 'unknown')
    bzlf1_type = region_results.get('BZLF1', {}).get('combined_call', 'unknown')

    recombination_status = 'no_recombination'
    recombination_note = ''

    if ebna2_type != 'unknown' and ebna3_type != 'unknown':
        if ebna2_type != ebna3_type and ebna2_type != 'ambiguous' and ebna3_type != 'ambiguous':
            recombination_status = 'intertypic_recombinant'
            recombination_note = (
                f'EBNA2 = {ebna2_type}, EBNA3 = {ebna3_type} — '
                f'discordant region types suggest intertypic recombination. '
                f'Breakpoint likely between EBNA2 and EBNA3 loci.'
            )
        elif ebna2_type == 'ambiguous' or ebna3_type == 'ambiguous':
            recombination_status = 'possible_recombination'
            recombination_note = (
                f'EBNA2 = {ebna2_type}, EBNA3 = {ebna3_type} — '
                f'one region ambiguous, cannot rule out recombination.'
            )

    # Check for mixed infection (both types in same region)
    mixed_signals = []
    for region_name, rr in region_results.items():
        sd = rr.get('snp_details', {})
        t1f = sd.get('t1_fraction', 0)
        t2f = sd.get('t2_fraction', 0)
        if sd.get('callable_snps', 0) >= 10:
            if 0.2 < t1f < 0.8 and 0.2 < t2f < 0.8:
                mixed_signals.append(
                    f'{region_name} (t1={t1f:.0%}, t2={t2f:.0%})'
                )
    if mixed_signals:
        if recombination_status == 'no_recombination':
            recombination_status = 'possible_mixed_infection'
        recombination_note += (
            f' Mixed type signals in: {", ".join(mixed_signals)}. '
            f'This may indicate co-infection with both EBV types.'
        )

    # BZLF1 Zp-V3 check
    bzlf1_note = ''
    if bzlf1_type == 'type2':
        bzlf1_note = 'BZLF1 region shows type-2 pattern (Zp-V3 variant likely present).'
    elif bzlf1_type == 'type1':
        bzlf1_note = 'BZLF1 region shows type-1 pattern (Zp-P variant likely).'
    else:
        bzlf1_note = 'BZLF1 region ambiguous or insufficient data.'

    result = {
        'sample': snakemake.wildcards.sample,
        'overall_type': overall_typing.get('ebv_type', 'unknown'),
        'recombination_status': recombination_status,
        'recombination_note': recombination_note,
        'bzlf1_note': bzlf1_note,
        'region_typing': region_results,
    }

    # Write JSON
    with open(snakemake.output.recomb_json, 'w') as out:
        json.dump(result, out, indent=2)

    # Write TSV
    with open(snakemake.output.recomb_tsv, 'w') as out:
        out.write('field\tvalue\n')
        out.write(f"sample\t{result['sample']}\n")
        out.write(f"overall_type\t{result['overall_type']}\n")
        out.write(f"recombination_status\t{result['recombination_status']}\n")
        out.write(f"ebna2_type\t{ebna2_type}\n")
        out.write(f"ebna3_type\t{ebna3_type}\n")
        out.write(f"bzlf1_type\t{bzlf1_type}\n")
        out.write(f"note\t{result['recombination_note']}\n")
        for region_name, rr in region_results.items():
            out.write(f"{region_name}_snp_call\t{rr['snp_call']}\n")
            out.write(f"{region_name}_cov_call\t{rr['coverage_call']}\n")
            out.write(f"{region_name}_callable_snps\t{rr['snp_details'].get('callable_snps', 0)}\n")

    print(f"[recombination_detection] Sample: {result['sample']}")
    print(f"  Recombination status: {recombination_status}")
    print(f"  EBNA2 type: {ebna2_type}")
    print(f"  EBNA3 type: {ebna3_type}")
    print(f"  BZLF1 type: {bzlf1_type}")
    if recombination_note:
        print(f"  Note: {recombination_note}")


if __name__ == '__main__' or 'snakemake' in globals():
    if 'snakemake' in globals():
        main()
