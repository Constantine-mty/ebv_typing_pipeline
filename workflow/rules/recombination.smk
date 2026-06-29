# recombination.smk — Recombination / mixed infection detection
# -------------------------------------------------------------------
# Independently type the EBNA2 region and the EBNA3 region.
# If EBNA2 type != EBNA3 type → flag as potential intertypic recombinant.
# Also check BZLF1 promoter variant (Zp-V3 = type 2 hallmark).
#
# Rationale:
#   - Intertypic recombinants are well-documented in EBV [Midgley et al.
#     2000; Yao et al. 1996; Blazquez et al. 2025].  They have, e.g.,
#     type-1 EBNA2 with type-2 EBNA3, or vice versa.
#   - Detecting recombination requires typing each discriminatory region
#     independently and checking for discordance.
#   - The BZLF1 Zp-V3 promoter variant is strictly associated with EBV-2
#     [Blazquez et al. 2025; Correia et al. 2017] and provides an
#     independent marker.

rule detect_recombination:
    """Detect intertypic recombination or mixed infection."""
    input:
        ebv1_vcf = "results/typing/{sample}_ebv1_snps.vcf.gz",
        ebv2_vcf = "results/typing/{sample}_ebv2_snps.vcf.gz",
        ebv1_cov = "results/typing/{sample}_ebv1_coverage.tsv",
        ebv2_cov = "results/typing/{sample}_ebv2_coverage.tsv",
        type_snps = config["refs"]["type_specific_snps"],
        typing_json = "results/typing/{sample}_typing.json",
    output:
        recomb_json = "results/recombination/{sample}_recombination.json",
        recomb_tsv = "results/recombination/{sample}_recombination.tsv",
    script:
        "../scripts/recombination_detection.py"
