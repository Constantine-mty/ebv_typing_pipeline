# qc.smk — Quality control and adapter trimming with fastp
# -------------------------------------------------------------------
# Input  : raw FASTQ files (from config[samples])
# Output : trimmed FASTQ + fastp JSON/HTML report
#
# Rationale: RNA-seq reads may contain adapter sequences and low-quality
# bases that would interfere with alignment. fastp performs both QC and
# trimming in a single pass and is the standard tool for Illumina data.

rule all_qc:
    input:
        expand("results/qc/{sample}_R1.trimmed.fastq.gz", sample=SAMPLES),
        expand("results/qc/{sample}_R2.trimmed.fastq.gz", sample=PAIRED_SAMPLES),
        expand("results/qc/{sample}_S.trimmed.fastq.gz", sample=SINGLE_SAMPLES),
        expand("results/qc/{sample}_fastp.json", sample=SAMPLES),


rule fastp_paired:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["r1"],
        r2 = lambda wc: config["samples"][wc.sample]["r2"],
    output:
        r1 = "results/qc/{sample}_R1.trimmed.fastq.gz",
        r2 = "results/qc/{sample}_R2.trimmed.fastq.gz",
        json = "results/qc/{sample}_fastp.json",
        html = "results/qc/{sample}_fastp.html",
    params:
        qual = lambda w: config["fastp"]["qualified_quality_phred"],
        lenreq = lambda w: config["fastp"]["length_required"],
    log:
        "results/qc/logs/{sample}_fastp.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -q {params.qual} \
            --length_required {params.lenreq} \
            --detect_adapter_for_pe \
            --thread {threads} \
            --json {output.json} \
            --html {output.html} \
            2> {log}
        """


rule fastp_single:
    input:
        r1 = lambda wc: config["samples"][wc.sample]["r1"],
    output:
        s = "results/qc/{sample}_S.trimmed.fastq.gz",
        json = "results/qc/{sample}_fastp.json",
        html = "results/qc/{sample}_fastp.html",
    params:
        qual = lambda w: config["fastp"]["qualified_quality_phred"],
        lenreq = lambda w: config["fastp"]["length_required"],
    log:
        "results/qc/logs/{sample}_fastp.log",
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        fastp \
            -i {input.r1} \
            -o {output.s} \
            -q {params.qual} \
            --length_required {params.lenreq} \
            --thread {threads} \
            --json {output.json} \
            --html {output.html} \
            2> {log}
        """
