# alignment.smk — Human read filtering + dual EBV alignment
# -------------------------------------------------------------------
# Step 1: Align trimmed reads to GRCh38; retain unmapped reads (EBV + other).
# Step 2: Align human-unmapped reads separately to EBV-1 and EBV-2 references.
#
# Rationale:
#   - Human reads dominate RNA-seq from infected cells (>99%).  Filtering
#     them first reduces noise and speeds EBV alignment.
#   - Dual-reference alignment is the core typing strategy: reads from a
#     type-1 infection will align better to the EBV-1 reference, and vice
#     versa, because EBNA2/EBNA3 differ by up to 30% between types.
#   - HISAT2 is chosen over STAR for memory efficiency (~8 GB vs ~30 GB
#     for the human genome index).

# Helper: get R1 trimmed FASTQ path for a sample
def get_r1(wc):
    sample = wc.sample
    if config["samples"][sample].get("r2"):
        return f"results/qc/{sample}_R1.trimmed.fastq.gz"
    return f"results/qc/{sample}_S.trimmed.fastq.gz"

# Helper: get R2 trimmed FASTQ path for a sample (empty list for single-end)
def get_r2(wc):
    sample = wc.sample
    if config["samples"][sample].get("r2"):
        return f"results/qc/{sample}_R2.trimmed.fastq.gz"
    return []


rule align_to_human:
    """Align reads to GRCh38; extract unmapped reads for EBV alignment."""
    input:
        r1 = get_r1,
        r2 = get_r2,
        human_idx = "resources/grch38.fa.hisat2_idx",
    output:
        r1 = "results/alignment/{sample}_ebv_input_R1.fastq.gz",
        r2 = "results/alignment/{sample}_ebv_input_R2.fastq.gz",
    params:
        human_fa = config["refs"]["human_fa"],
        extra = config["hisat2"]["human_filter_extra"],
    threads: 8
    conda:
        "../envs/alignment.yaml"
    log:
        "results/alignment/logs/{sample}_human_align.log",
    shell:
        r"""
        # Determine paired-end vs single-end
        if [ -f "{input.r2}" ]; then
            echo "Paired-end alignment to human genome..."
            hisat2 {params.extra} \
                -x {params.human_fa} \
                -p {threads} \
                -1 {input.r1} -2 {input.r2} \
                --un-conc-gz results/alignment/{wildcards.sample}_ebv_input_R%.fastq.gz \
                --summary-file {log} \
                2> {log}.err \
                > /dev/null

            # HISAT2 normally creates both files, but keep empty valid gzip
            # placeholders so downstream rules fail only on real aligner errors.
            for fq in {output.r1} {output.r2}; do
                if [ ! -e "$fq" ]; then
                    printf '' | gzip -c > "$fq"
                fi
            done
        else
            echo "Single-end alignment to human genome..."
            hisat2 {params.extra} \
                -x {params.human_fa} \
                -p {threads} \
                -U {input.r1} \
                --un-gz {output.r1} \
                --summary-file {log} \
                2> {log}.err \
                > /dev/null
            if [ ! -e "{output.r1}" ]; then
                printf '' | gzip -c > {output.r1}
            fi
            touch {output.r2}
        fi
        """


rule align_to_ebv1:
    """Align human-unmapped reads to EBV-1 reference (NC_007605.1)."""
    input:
        r1 = "results/alignment/{sample}_ebv_input_R1.fastq.gz",
        r2 = "results/alignment/{sample}_ebv_input_R2.fastq.gz",
        ebv1_idx = "resources/ebv1_reference.fa.hisat2_idx",
    output:
        bam = "results/alignment/{sample}_ebv1.bam",
        stats = "results/alignment/{sample}_ebv1_stats.txt",
    params:
        ebv1_fa = config["refs"]["ebv1_fa"],
        extra = config["hisat2"]["ebv_extra"],
    threads: 4
    conda:
        "../envs/alignment.yaml"
    log:
        "results/alignment/logs/{sample}_ebv1_align.log",
    shell:
        r"""
        if [ -s "{input.r2}" ]; then
            hisat2 {params.extra} \
                -x {params.ebv1_fa} \
                -p {threads} \
                -1 {input.r1} -2 {input.r2} \
                --summary-file {output.stats} \
                2> {log} \
                | samtools sort -@ {threads} -o {output.bam}
        else
            hisat2 {params.extra} \
                -x {params.ebv1_fa} \
                -p {threads} \
                -U {input.r1} \
                --summary-file {output.stats} \
                2> {log} \
                | samtools sort -@ {threads} -o {output.bam}
        fi
        samtools index {output.bam}
        """


rule align_to_ebv2:
    """Align human-unmapped reads to EBV-2 reference (NC_009334.1)."""
    input:
        r1 = "results/alignment/{sample}_ebv_input_R1.fastq.gz",
        r2 = "results/alignment/{sample}_ebv_input_R2.fastq.gz",
        ebv2_idx = "resources/ebv2_reference.fa.hisat2_idx",
    output:
        bam = "results/alignment/{sample}_ebv2.bam",
        stats = "results/alignment/{sample}_ebv2_stats.txt",
    params:
        ebv2_fa = config["refs"]["ebv2_fa"],
        extra = config["hisat2"]["ebv_extra"],
    threads: 4
    conda:
        "../envs/alignment.yaml"
    log:
        "results/alignment/logs/{sample}_ebv2_align.log",
    shell:
        r"""
        if [ -s "{input.r2}" ]; then
            hisat2 {params.extra} \
                -x {params.ebv2_fa} \
                -p {threads} \
                -1 {input.r1} -2 {input.r2} \
                --summary-file {output.stats} \
                2> {log} \
                | samtools sort -@ {threads} -o {output.bam}
        else
            hisat2 {params.extra} \
                -x {params.ebv2_fa} \
                -p {threads} \
                -U {input.r1} \
                --summary-file {output.stats} \
                2> {log} \
                | samtools sort -@ {threads} -o {output.bam}
        fi
        samtools index {output.bam}
        """
