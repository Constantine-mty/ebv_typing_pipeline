# references.smk — Reference genome index preparation
# -------------------------------------------------------------------
# This rule builds HISAT2 indices for the human, EBV-1, and EBV-2
# reference genomes.  The EBV references and their GTF/SNP tables are
# pre-built by prepare_references.py (run once during setup).
#
# The human reference (GRCh38) must be provided by the user.  If the
# file does not exist, the pipeline will error with a clear message.

rule build_hisat2_indices:
    """Build HISAT2 indices for human, EBV-1, and EBV-2 references."""
    input:
        human_fa = config["refs"]["human_fa"],
        ebv1_fa  = config["refs"]["ebv1_fa"],
        ebv2_fa  = config["refs"]["ebv2_fa"],
    output:
        human_idx = touch("resources/grch38.fa.hisat2_idx"),
        ebv1_idx  = touch("resources/ebv1_reference.fa.hisat2_idx"),
        ebv2_idx  = touch("resources/ebv2_reference.fa.hisat2_idx"),
    conda:
        "../envs/alignment.yaml"
    shell:
        r"""
        # Check human reference exists
        if [ ! -f "{input.human_fa}" ]; then
            echo "ERROR: Human reference genome not found at {input.human_fa}"
            echo "Please download GRCh38 and place it there, or update config.yaml."
            exit 1
        fi

        # Build HISAT2 indices (EBV genomes are tiny, human takes ~30 min)
        echo "Building EBV-1 HISAT2 index..."
        hisat2-build -p {threads} {input.ebv1_fa} {input.ebv1_fa}
        echo "Building EBV-2 HISAT2 index..."
        hisat2-build -p {threads} {input.ebv2_fa} {input.ebv2_fa}
        echo "Building human HISAT2 index (this may take ~30 min)..."
        hisat2-build -p {threads} {input.human_fa} {input.human_fa}
        echo "All HISAT2 indices built."
        """
