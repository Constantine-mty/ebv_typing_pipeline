# typing.smk — EBV type classification (coverage + SNP concordance)
# -------------------------------------------------------------------
# Step 1: Call SNPs on EBV-1 and EBV-2 alignments (bcftools).
# Step 2: Compute coverage on type-discriminatory regions (bedtools).
# Step 3: Classify type using combined coverage + SNP concordance scoring.
#
# Rationale:
#   - Coverage-based: reads from a type-1 infection align better to the
#     EBV-1 reference (especially in divergent EBNA2/EBNA3 regions),
#     producing higher read counts there.
#   - SNP concordance: at positions where EBV-1 and EBV-2 references
#     differ, the called allele should match the infecting type's
#     reference.  This is the most direct evidence of type identity.
#   - Combining both methods provides robustness and allows confidence
#     assessment.

rule call_snps_ebv1:
    """Call SNPs on the EBV-1 alignment."""
    input:
        bam = "results/alignment/{sample}_ebv1.bam",
        fa = config["refs"]["ebv1_fa"],
    output:
        vcf = "results/typing/{sample}_ebv1_snps.vcf.gz",
    params:
        min_bq = config["bcftools"]["min_base_quality"],
        min_mq = config["bcftools"]["min_mapping_quality"],
        min_dp = config["bcftools"]["min_depth"],
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        bcftools mpileup \
            -f {input.fa} \
            -q {params.min_mq} \
            -Q {params.min_bq} \
            -a FORMAT/DP \
            {input.bam} \
            | bcftools call -mv -Oz \
            | bcftools view -i 'DP>={params.min_dp}' -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule call_snps_ebv2:
    """Call SNPs on the EBV-2 alignment."""
    input:
        bam = "results/alignment/{sample}_ebv2.bam",
        fa = config["refs"]["ebv2_fa"],
    output:
        vcf = "results/typing/{sample}_ebv2_snps.vcf.gz",
    params:
        min_bq = config["bcftools"]["min_base_quality"],
        min_mq = config["bcftools"]["min_mapping_quality"],
        min_dp = config["bcftools"]["min_depth"],
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        bcftools mpileup \
            -f {input.fa} \
            -q {params.min_mq} \
            -Q {params.min_bq} \
            -a FORMAT/DP \
            {input.bam} \
            | bcftools call -mv -Oz \
            | bcftools view -i 'DP>={params.min_dp}' -Oz -o {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule compute_coverage:
    """Compute per-base coverage on both EBV alignments."""
    input:
        ebv1_bam = "results/alignment/{sample}_ebv1.bam",
        ebv2_bam = "results/alignment/{sample}_ebv2.bam",
        ebv1_gtf = config["refs"]["ebv1_gtf"],
        ebv2_gtf = config["refs"]["ebv2_gtf"],
    output:
        ebv1_cov = "results/typing/{sample}_ebv1_coverage.tsv",
        ebv2_cov = "results/typing/{sample}_ebv2_coverage.tsv",
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        # Per-gene read counts via bedtools coverage
        # Convert GTF to BED (gene-level)
        awk -F'\t' '$3=="gene" {{split($9,a,"\""); print $1"\t"$4"\t"$5"\t"a[2]"\t.\t"$7}}' \
            {input.ebv1_gtf} > results/typing/{wildcards.sample}_ebv1_genes.bed
        awk -F'\t' '$3=="gene" {{split($9,a,"\""); print $1"\t"$4"\t"$5"\t"a[2]"\t.\t"$7}}' \
            {input.ebv2_gtf} > results/typing/{wildcards.sample}_ebv2_genes.bed

        bedtools coverage \
            -a results/typing/{wildcards.sample}_ebv1_genes.bed \
            -b {input.ebv1_bam} \
            -counts > {output.ebv1_cov}

        bedtools coverage \
            -a results/typing/{wildcards.sample}_ebv2_genes.bed \
            -b {input.ebv2_bam} \
            -counts > {output.ebv2_cov}
        """


rule classify_type:
    """Classify EBV type using combined coverage + SNP concordance scoring."""
    input:
        ebv1_vcf = "results/typing/{sample}_ebv1_snps.vcf.gz",
        ebv2_vcf = "results/typing/{sample}_ebv2_snps.vcf.gz",
        ebv1_cov = "results/typing/{sample}_ebv1_coverage.tsv",
        ebv2_cov = "results/typing/{sample}_ebv2_coverage.tsv",
        type_snps = config["refs"]["type_specific_snps"],
        gene_categories = config["refs"]["gene_categories"],
        latency_json = "results/expression/{sample}_latency.json",
    output:
        typing_json = "results/typing/{sample}_typing.json",
        typing_tsv = "results/typing/{sample}_typing.tsv",
    params:
        min_ebv_reads = config["typing"]["min_ebv_reads"],
        min_ebna2_reads = config["typing"]["min_ebna2_reads"],
        min_ebna3_reads = config["typing"]["min_ebna3_reads"],
        min_type_snps = config["typing"]["min_type_snps"],
        high_conf = config["typing"]["high_confidence_threshold"],
        cov_log2 = config["typing"]["coverage_log2_threshold"],
    script:
        "../scripts/type_classification.py"
