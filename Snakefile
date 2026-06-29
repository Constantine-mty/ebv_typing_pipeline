# =============================================================================
# EBV Type 1/2 Identification Pipeline — Main Snakefile
# =============================================================================
# Identifies EBV type (1 or 2) from RNA-seq data of EBV-infected human cells.
#
# Usage:
#   snakemake --use-conda -j 4 --configfile config/config.yaml
#
# For a single sample:
#   snakemake --use-conda -j 4 --configfile config/config.yaml results/sample1_EBV_typing_report.md
#
# Requirements:
#   - Edit config/config.yaml to point at your FASTQ files and human reference
#   - EBV references and resource files are pre-built in resources/
# =============================================================================

import os
import sys
from pathlib import Path

# -----------------------------------------------------------------------------
# Load configuration
# -----------------------------------------------------------------------------
configfile: "config/config.yaml"

# -----------------------------------------------------------------------------
# Resolve sample lists
# -----------------------------------------------------------------------------
SAMPLES = list(config["samples"].keys())

# Determine paired-end vs single-end samples
PAIRED_SAMPLES = [s for s in SAMPLES if "r2" in config["samples"][s]]
SINGLE_SAMPLES = [s for s in SAMPLES if "r2" not in config["samples"][s]]

# -----------------------------------------------------------------------------
# Output directories
# -----------------------------------------------------------------------------
OUTDIRS = config.get("output_dirs", {})
QC_DIR = OUTDIRS.get("qc", "results/qc")
ALIGN_DIR = OUTDIRS.get("alignment", "results/alignment")
EXPR_DIR = OUTDIRS.get("expression", "results/expression")
TYPING_DIR = OUTDIRS.get("typing", "results/typing")
RECOMB_DIR = OUTDIRS.get("recombination", "results/recombination")
TRANS_DIR = OUTDIRS.get("transcriptome", "results/transcriptome")
FIG_DIR = OUTDIRS.get("figures", "results/figures")

# Ensure directories exist
for d in [QC_DIR, ALIGN_DIR, EXPR_DIR, TYPING_DIR, RECOMB_DIR, TRANS_DIR, FIG_DIR]:
    Path(d).mkdir(parents=True, exist_ok=True)
    Path(f"{d}/logs").mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# Import all rule modules
# -----------------------------------------------------------------------------
include: "workflow/rules/references.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/expression.smk"
include: "workflow/rules/typing.smk"
include: "workflow/rules/recombination.smk"
include: "workflow/rules/transcriptome.smk"
include: "workflow/rules/report.smk"

# -----------------------------------------------------------------------------
# Default target: generate reports for all samples
# -----------------------------------------------------------------------------
rule all:
    input:
        expand("results/{sample}_EBV_typing_report.md", sample=SAMPLES),
        expand("results/{sample}_EBV_typing_summary.json", sample=SAMPLES),
