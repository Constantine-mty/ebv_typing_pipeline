# expression.smk — EBV gene quantification + latency assessment
# -------------------------------------------------------------------
# Step 1: Count reads per EBV gene using featureCounts on both EBV-1
#         and EBV-2 alignments.
# Step 2: Assess latency type from expression patterns (Python script).
#
# Rationale:
#   - Gene-level counts tell us which EBV genes are expressed, which
#     determines the latency type (I, II, or III).
#   - Latency type determines which typing strategy is viable:
#     Latency III → EBNA2-based typing (highest confidence)
#     Latency II  → EBNA3-based typing
#     Latency I   → genome-wide SNP concordance (lower confidence)

rule count_ebv1_genes:
    """Count reads per EBV gene on the EBV-1 alignment."""
    input:
        bam = "results/alignment/{sample}_ebv1.bam",
        gtf = config["refs"]["ebv1_gtf"],
    output:
        counts = "results/expression/{sample}_ebv1_counts.txt",
        summary = "results/expression/{sample}_ebv1_counts.txt.summary",
    params:
        outprefix = "results/expression/{sample}_ebv1_counts",
    threads: 4
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        featureCounts \
            -a {input.gtf} \
            -o {output.counts} \
            -T {threads} \
            -p --countReadPairs \
            -t gene \
            -g gene_id \
            -s 0 \
            {input.bam} 2> {output.counts}.log
        """


rule count_ebv2_genes:
    """Count reads per EBV gene on the EBV-2 alignment."""
    input:
        bam = "results/alignment/{sample}_ebv2.bam",
        gtf = config["refs"]["ebv2_gtf"],
    output:
        counts = "results/expression/{sample}_ebv2_counts.txt",
        summary = "results/expression/{sample}_ebv2_counts.txt.summary",
    threads: 4
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        featureCounts \
            -a {input.gtf} \
            -o {output.counts} \
            -T {threads} \
            -p --countReadPairs \
            -t gene \
            -g gene_id \
            -s 0 \
            {input.bam} 2> {output.counts}.log
        """


rule assess_latency:
    """Assess EBV latency type from gene expression patterns."""
    input:
        ebv1_counts = "results/expression/{sample}_ebv1_counts.txt",
        ebv2_counts = "results/expression/{sample}_ebv2_counts.txt",
        ebv1_stats = "results/alignment/{sample}_ebv1_stats.txt",
        ebv2_stats = "results/alignment/{sample}_ebv2_stats.txt",
        gene_categories = config["refs"]["gene_categories"],
    output:
        latency_json = "results/expression/{sample}_latency.json",
        latency_tsv = "results/expression/{sample}_latency.tsv",
    params:
        expressed_tpm = config["latency"]["expressed_tpm"],
        absent_tpm = config["latency"]["absent_tpm"],
    script:
        "../scripts/assess_latency.py"
