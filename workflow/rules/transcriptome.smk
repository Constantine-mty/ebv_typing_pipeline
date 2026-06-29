# transcriptome.smk — EBV transcriptome overview
# -------------------------------------------------------------------
# Generate visualisations and tables summarising EBV gene expression,
# latent vs lytic proportions, and genome-wide coverage.
#
# Outputs:
#   - Bar plot of EBV gene TPM (sorted by expression)
#   - Pie chart of latent vs lytic expression proportion
#   - Genome-wide coverage plot
#   - Type-discriminatory region coverage comparison (EBV-1 vs EBV-2)
#   - Per-gene expression table

rule transcriptome_overview:
    """Generate EBV transcriptome overview plots and tables."""
    input:
        ebv1_counts = "results/expression/{sample}_ebv1_counts.txt",
        ebv2_counts = "results/expression/{sample}_ebv2_counts.txt",
        ebv1_bam = "results/alignment/{sample}_ebv1.bam",
        ebv2_bam = "results/alignment/{sample}_ebv2.bam",
        ebv1_gtf = config["refs"]["ebv1_gtf"],
        ebv2_gtf = config["refs"]["ebv2_gtf"],
        gene_categories = config["refs"]["gene_categories"],
        latency_json = "results/expression/{sample}_latency.json",
        typing_json = "results/typing/{sample}_typing.json",
    output:
        gene_expression_tsv = "results/transcriptome/{sample}_gene_expression.tsv",
        barplot = "results/figures/{sample}_gene_expression_barplot.png",
        piechart = "results/figures/{sample}_latent_lytic_pie.png",
        coverage_plot = "results/figures/{sample}_genome_coverage.png",
        region_comparison = "results/figures/{sample}_region_coverage_comparison.png",
    script:
        "../scripts/transcriptome_overview.py"
