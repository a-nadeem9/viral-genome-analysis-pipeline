# workflow/rules/kraken.smk
# ------------------------------------------------------------
# Taxonomic screening with Kraken 2  ➜  summary via MultiQC
# ------------------------------------------------------------

# `KRAKEN_DIR`, `QC_TRIM_DIR`, and `SAMPLES` come from Snakefile globals
KRAKEN_DB  = config["kraken2_db"]
LOG_DIR    = "logs/kraken2"          # single fixed log folder

# ------------------------------------------------------------
# When a DB is supplied ­→ run Kraken + MultiQC
# ------------------------------------------------------------
if KRAKEN_DB:

    rule kraken2_classify:
        """Classify trimmed paired reads with Kraken 2."""
        threads: 8
        input:
            fq1 = lambda wc: f"{QC_TRIM_DIR}/{wc.sample}_1.trimmed.fastq.gz",
            fq2 = lambda wc: f"{QC_TRIM_DIR}/{wc.sample}_2.trimmed.fastq.gz",
        output:
            report = f"{KRAKEN_DIR}/{{sample}}.report",
            kraken = f"{KRAKEN_DIR}/{{sample}}.kraken",
        params:
            db = KRAKEN_DB
        conda: "../envs/kraken2.yaml"
        log:   f"{LOG_DIR}/{{sample}}.log"
        shell: r"""
            kraken2 --db {params.db} --threads {threads} --paired \
                    --report {output.report} --output {output.kraken} \
                    {input.fq1} {input.fq2} \
                    2> {log}
        """

    rule kraken_multiqc:
        """Aggregate Kraken reports with MultiQC."""
        input:  expand(f"{KRAKEN_DIR}/{{sample}}.report", sample=SAMPLES)
        output: html = f"{KRAKEN_DIR}/multiqc/kraken_report.html"
        conda:  "../envs/multiqc.yaml"
        log:    f"{LOG_DIR}/multiqc.log"
        shell: r"""
            multiqc {KRAKEN_DIR} \
                    -o {KRAKEN_DIR}/multiqc \
                    -n kraken_report.html \
                    2> {log}
        """

    # convenience target so you can run:  snakemake screen
    rule screen:
        input:  rules.kraken_multiqc.output.html

# ------------------------------------------------------------
# No DB given ­→ define a no-op `screen` target so pipeline still runs
# ------------------------------------------------------------
else:

    rule screen:
        shell: "echo 'kraken2_db not set – skipping taxonomic screening.'"
