# workflow/rules/variability.smk
# ------------------------------------------------------------
# Compute per-site Shannon entropy + smoothed plot
# ------------------------------------------------------------
# Globals provided by Snakefile:
#   MSA_DIR  – path to multiple-sequence alignment folder
#   VAR_DIR  – folder to write variability outputs

rule compute_variability:
    input:
        msa = f"{MSA_DIR}/aligned_scaffolds.fasta"
    output:
        txt      = f"{VAR_DIR}/variability.txt",
        windowed = f"{VAR_DIR}/variability_windowed.txt",
        plot     = f"{VAR_DIR}/variability_plot.png",
    conda: "../envs/variability.yaml"
    threads: 1
    script: "../scripts/compute_variability.py"

