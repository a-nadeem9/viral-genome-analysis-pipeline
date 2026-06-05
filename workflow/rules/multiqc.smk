# workflow/rules/multiqc.smk
# ------------------------------------------------------------
# Aggregate per-sample Fastp JSONs into a single MultiQC report
# ------------------------------------------------------------

# `SAMPLES`   – list of sample IDs          (from Snakefile)
# `MULTIQC_DIR` – fixed results directory  (from Snakefile)

rule multiqc:
    input:
        # depend on every Fastp JSON
        jsons = expand("logs/fastp/{sample}_fastp.json", sample=SAMPLES)
    output:
        report = f"{MULTIQC_DIR}/fastp_report.html"
    conda: "../envs/multiqc.yaml"
    log:   "logs/fastp/multiqc.log"
    shell: r"""
        multiqc -f logs/fastp \
               -o {MULTIQC_DIR} \
               -n fastp_report.html \
               2> {log}
    """
