# report.smk — Final report generation
# -------------------------------------------------------------------
# Assemble all results into a single Markdown report and a
# machine-readable JSON summary.

rule generate_report:
    """Generate the final EBV typing report."""
    input:
        latency_json = "results/expression/{sample}_latency.json",
        typing_json = "results/typing/{sample}_typing.json",
        recomb_json = "results/recombination/{sample}_recombination.json",
        gene_expression_tsv = "results/transcriptome/{sample}_gene_expression.tsv",
        barplot = "results/figures/{sample}_gene_expression_barplot.png",
        piechart = "results/figures/{sample}_latent_lytic_pie.png",
        coverage_plot = "results/figures/{sample}_genome_coverage.png",
        region_comparison = "results/figures/{sample}_region_coverage_comparison.png",
        ebv1_stats = "results/alignment/{sample}_ebv1_stats.txt",
        ebv2_stats = "results/alignment/{sample}_ebv2_stats.txt",
    output:
        report = "results/{sample}_EBV_typing_report.md",
        summary_json = "results/{sample}_EBV_typing_summary.json",
    script:
        "../scripts/generate_report.py"
