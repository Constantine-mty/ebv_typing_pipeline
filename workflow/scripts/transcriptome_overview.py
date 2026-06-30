#!/usr/bin/env python3
"""
transcriptome_overview.py
==========================
Generate EBV transcriptome overview plots and tables.

Outputs:
  - gene_expression.tsv : per-gene count, TPM, category
  - gene_expression_barplot.png : bar plot of EBV gene TPM (sorted)
  - latent_lytic_pie.png : pie chart of latent vs lytic expression
  - genome_coverage.png : read depth across EBV genome
  - region_coverage_comparison.png : EBV-1 vs EBV-2 coverage in
    type-discriminatory regions

Inputs (via snakemake):
  - ebv1_counts, ebv2_counts : featureCounts output
  - ebv1_bam, ebv2_bam : sorted BAM files
  - ebv1_gtf, ebv2_gtf : gene annotations
  - gene_categories : gene → category mapping
  - latency_json, typing_json : prior results
"""

import json
import os
import re
import sys
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import pysam

# Phylo color palette
PHYLO_COLORS = {
    'black': '#000000',
    'cream': '#ECE9E2',
    'offwhite': '#FAF9F3',
    'yellow': '#E9ED4C',
    'orange': '#FF9400',
    'green': '#75A025',
    'pink': '#FD9BED',
    'blue': '#0279EE',
}

# Font settings
matplotlib.rcParams['font.family'] = ['Liberation Sans', 'Arimo', 'DejaVu Sans']
matplotlib.rcParams['svg.fonttype'] = 'none'


def read_featurecounts(path):
    """Read featureCounts output."""
    df = pd.read_csv(path, sep='\t', comment='#')
    count_col = df.columns[-1]
    df = df.rename(columns={count_col: 'count', 'Geneid': 'gene', 'Length': 'length'})
    return df[['gene', 'length', 'count']]


def compute_tpm(df):
    """Add TPM column."""
    df = df.copy()
    df['rpk'] = df['count'] / df['length'].replace(0, 1)
    total_rpk = df['rpk'].sum()
    df['tpm'] = df['rpk'] / total_rpk * 1e6 if total_rpk > 0 else 0.0
    return df


def load_gene_categories(path):
    """Load gene → category mapping."""
    cat = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                cat[parts[0]] = parts[1]
    return cat


def plot_gene_expression_barplot(df, gene_categories, out_path, sample_name):
    """Bar plot of EBV gene TPM, sorted by expression, colored by category."""
    df = df[df['tpm'] > 0].sort_values('tpm', ascending=True)
    if df.empty:
        # Create empty placeholder
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'No EBV gene expression detected',
                ha='center', va='center', fontsize=14)
        ax.set_title(f'{sample_name}: EBV Gene Expression')
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close()
        return

    colors = []
    for gene in df['gene']:
        cat = gene_categories.get(gene, 'other')
        if cat == 'latent':
            colors.append(PHYLO_COLORS['blue'])
        elif cat == 'lytic':
            colors.append(PHYLO_COLORS['orange'])
        else:
            colors.append(PHYLO_COLORS['green'])

    fig, ax = plt.subplots(figsize=(10, max(6, len(df) * 0.3)))
    ax.barh(df['gene'], df['tpm'], color=colors, edgecolor='none')
    ax.set_xlabel('TPM', fontsize=12)
    ax.set_title(f'{sample_name}: EBV Gene Expression (TPM)', fontsize=13)
    ax.set_xscale('symlog', linthresh=1)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=PHYLO_COLORS['blue'], label='Latent'),
        Patch(facecolor=PHYLO_COLORS['orange'], label='Lytic'),
        Patch(facecolor=PHYLO_COLORS['green'], label='Other'),
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_latent_lytic_pie(df, gene_categories, out_path, sample_name):
    """Pie chart of latent vs lytic expression proportion."""
    latent_tpm = 0
    lytic_tpm = 0
    other_tpm = 0
    for _, row in df.iterrows():
        cat = gene_categories.get(row['gene'], 'other')
        if cat == 'latent':
            latent_tpm += row['tpm']
        elif cat == 'lytic':
            lytic_tpm += row['tpm']
        else:
            other_tpm += row['tpm']

    total = latent_tpm + lytic_tpm + other_tpm
    if total == 0:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.text(0.5, 0.5, 'No EBV expression detected',
                ha='center', va='center', fontsize=14)
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close()
        return

    sizes = [latent_tpm, lytic_tpm, other_tpm]
    labels = ['Latent', 'Lytic', 'Other']
    colors = [PHYLO_COLORS['blue'], PHYLO_COLORS['orange'], PHYLO_COLORS['green']]
    # Filter out zero entries
    sizes_f = [s for s in sizes if s > 0]
    labels_f = [l for s, l in zip(sizes, labels) if s > 0]
    colors_f = [c for s, c in zip(sizes, colors) if s > 0]

    fig, ax = plt.subplots(figsize=(6, 6))
    wedges, texts, autotexts = ax.pie(
        sizes_f, labels=labels_f, colors=colors_f,
        autopct='%1.1f%%', startangle=90, textprops={'fontsize': 11}
    )
    ax.set_title(f'{sample_name}: Latent vs Lytic Expression', fontsize=13)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


def compute_genome_coverage_profile(bam_path):
    """Return windowed read depth across one EBV genome alignment."""
    bam = pysam.AlignmentFile(bam_path, 'rb')
    ref_name = bam.references[0]
    ref_length = bam.get_reference_length(ref_name)

    # Compute per-base depth (sample every 100 bp for speed)
    window = 100
    positions = []
    depths = []
    for pos in range(1, ref_length + 1, window):
        end = min(pos + window - 1, ref_length)
        depth = bam.count_coverage(ref_name, start=pos-1, end=end)
        # Sum across all 4 bases
        total = sum(sum(base_depths) for base_depths in depth)
        avg_depth = total / max(1, end - pos + 1)
        positions.append(pos)
        depths.append(avg_depth)
    bam.close()
    return positions, depths, ref_length


def draw_genome_coverage(ax, positions, depths, ref_length, sample_name, ref_name, color):
    """Draw one EBV genome coverage profile on an axis."""
    ax.fill_between(positions, depths, color=color, alpha=0.7)
    ax.set_xlabel(f'Position on {ref_name}', fontsize=11)
    ax.set_ylabel('Average read depth', fontsize=11)
    ax.set_title(f'{sample_name}: EBV Genome Coverage ({ref_name})', fontsize=13)
    ax.set_xlim(0, ref_length)


def plot_genome_coverage(bam_path, out_path, sample_name, ref_name, color):
    """Plot read depth across one EBV genome."""
    positions, depths, ref_length = compute_genome_coverage_profile(bam_path)

    fig, ax = plt.subplots(figsize=(12, 4))
    draw_genome_coverage(ax, positions, depths, ref_length, sample_name, ref_name, color)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


def plot_combined_genome_coverage(ebv1_bam, ebv2_bam, out_path, sample_name):
    """Plot EBV-1 and EBV-2 genome coverage in stacked panels."""
    p1, d1, len1 = compute_genome_coverage_profile(ebv1_bam)
    p2, d2, len2 = compute_genome_coverage_profile(ebv2_bam)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=False)
    draw_genome_coverage(
        axes[0], p1, d1, len1, sample_name,
        'EBV-1 (NC_007605.1)', PHYLO_COLORS['blue']
    )
    draw_genome_coverage(
        axes[1], p2, d2, len2, sample_name,
        'EBV-2 (NC_009334.1)', PHYLO_COLORS['orange']
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


def plot_region_coverage_comparison(cov1_df, cov2_df, out_path, sample_name):
    """Compare EBV-1 vs EBV-2 coverage in type-discriminatory regions."""
    type_genes = ['EBNA2', 'EBNA3A', 'EBNA3B', 'EBNA3C', 'EBNA3BC',
                  'BZLF1', 'BZLF2', 'BLLF1']

    counts1 = dict(zip(cov1_df['gene'], cov1_df['count']))
    counts2 = dict(zip(cov2_df['gene'], cov2_df['count']))

    genes_present = [g for g in type_genes
                     if g in counts1 or g in counts2]
    if not genes_present:
        genes_present = type_genes  # show all even if zero

    c1 = [counts1.get(g, 0) for g in genes_present]
    c2 = [counts2.get(g, 0) for g in genes_present]

    x = np.arange(len(genes_present))
    width = 0.35

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(x - width/2, c1, width, label='EBV-1 alignment',
           color=PHYLO_COLORS['blue'], edgecolor='none')
    ax.bar(x + width/2, c2, width, label='EBV-2 alignment',
           color=PHYLO_COLORS['orange'], edgecolor='none')
    ax.set_ylabel('Read count', fontsize=11)
    ax.set_title(f'{sample_name}: Coverage in Type-Discriminatory Regions',
                 fontsize=13)
    ax.set_xticks(x)
    ax.set_xticklabels(genes_present, rotation=45, ha='right')
    ax.legend(fontsize=10)
    ax.set_yscale('symlog', linthresh=1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved: {out_path}")


def main():
    sample = snakemake.wildcards.sample

    # Read expression data
    df1 = read_featurecounts(snakemake.input.ebv1_counts)
    df2 = read_featurecounts(snakemake.input.ebv2_counts)
    df1 = compute_tpm(df1)
    df2 = compute_tpm(df2)

    # Load gene categories
    gene_categories = load_gene_categories(snakemake.input.gene_categories)

    # Load latency and typing results
    with open(snakemake.input.latency_json) as fh:
        latency = json.load(fh)
    with open(snakemake.input.typing_json) as fh:
        typing = json.load(fh)

    # Use the primary alignment (more total reads)
    total1 = df1['count'].sum()
    total2 = df2['count'].sum()
    primary_df = df1 if total1 >= total2 else df2
    primary_ref = 'EBV-1 (NC_007605.1)' if total1 >= total2 else 'EBV-2 (NC_009334.1)'
    primary_bam = snakemake.input.ebv1_bam if total1 >= total2 else snakemake.input.ebv2_bam

    # Write gene expression table
    expr_table = primary_df[['gene', 'count', 'tpm']].copy()
    expr_table['category'] = expr_table['gene'].map(
        lambda g: gene_categories.get(g, 'other')
    )
    expr_table = expr_table.sort_values('tpm', ascending=False)
    expr_table.to_csv(snakemake.output.gene_expression_tsv, sep='\t', index=False)
    print(f"  Saved: {snakemake.output.gene_expression_tsv}")

    # Generate plots
    os.makedirs(os.path.dirname(snakemake.output.barplot), exist_ok=True)

    plot_gene_expression_barplot(
        primary_df, gene_categories, snakemake.output.barplot, sample
    )
    plot_latent_lytic_pie(
        primary_df, gene_categories, snakemake.output.piechart, sample
    )
    plot_genome_coverage(
        snakemake.input.ebv1_bam,
        snakemake.output.ebv1_coverage_plot,
        sample,
        'EBV-1 (NC_007605.1)',
        PHYLO_COLORS['blue'],
    )
    plot_genome_coverage(
        snakemake.input.ebv2_bam,
        snakemake.output.ebv2_coverage_plot,
        sample,
        'EBV-2 (NC_009334.1)',
        PHYLO_COLORS['orange'],
    )
    plot_combined_genome_coverage(
        snakemake.input.ebv1_bam,
        snakemake.input.ebv2_bam,
        snakemake.output.coverage_plot,
        sample,
    )

    # Region coverage comparison (need coverage tables — read from typing)
    # We re-read the coverage from the BAMs via bedtools-like approach
    # Actually, we can read the coverage TSVs from typing output
    cov1_path = f"results/typing/{sample}_ebv1_coverage.tsv"
    cov2_path = f"results/typing/{sample}_ebv2_coverage.tsv"
    if os.path.exists(cov1_path) and os.path.exists(cov2_path):
        cov1 = pd.read_csv(cov1_path, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'gene', 'score', 'strand', 'count'])
        cov2 = pd.read_csv(cov2_path, sep='\t', header=None,
                           names=['chr', 'start', 'end', 'gene', 'score', 'strand', 'count'])
        plot_region_coverage_comparison(cov1, cov2, snakemake.output.region_comparison, sample)
    else:
        # Create placeholder
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'Coverage data not available',
                ha='center', va='center', fontsize=14)
        plt.tight_layout()
        plt.savefig(snakemake.output.region_comparison, dpi=150)
        plt.close()

    print(f"[transcriptome_overview] Sample: {sample}")
    print(f"  Primary reference: {primary_ref}")
    print(f"  Total EBV genes detected: {(primary_df['count'] > 0).sum()}")


if __name__ == '__main__' or 'snakemake' in globals():
    if 'snakemake' in globals():
        main()
