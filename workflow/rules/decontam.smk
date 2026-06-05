# workflow/rules/decontam.smk
# ------------------------------------------------------------
# Optional contaminant removal via Bowtie2.
# Enabled when config["contaminants_fasta"] is not empty.
# Downstream rules use get_clean_reads() from common.smk
# to select trimmed or decontaminated reads as needed.
# ------------------------------------------------------------

if config["contaminants_fasta"]:

    rule build_contaminant_index:
        input:
            fasta = config["contaminants_fasta"]
        output:
            multiext("resources/contaminants/bt2_idx/contaminants",
                     ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params:
            prefix = "resources/contaminants/bt2_idx/contaminants"
        conda:  "../envs/bowtie2.yaml"
        log:    "logs/decont/bowtie2_build.log"
        shell:  r"""
            bowtie2-build {input.fasta} {params.prefix} > /dev/null 2> {log}
        """

    rule remove_contaminants:
        threads: 8
        input:
            fq1 = f"{QC_TRIM_DIR}/{{sample}}_1.trimmed.fastq.gz",
            fq2 = f"{QC_TRIM_DIR}/{{sample}}_2.trimmed.fastq.gz",
            idx = rules.build_contaminant_index.output,
        output:
            clean1 = f"{QC_DECONT_DIR}/{{sample}}_1.clean.fastq.gz",
            clean2 = f"{QC_DECONT_DIR}/{{sample}}_2.clean.fastq.gz",
        params:
            prefix = "resources/contaminants/bt2_idx/contaminants"
        conda: "../envs/bowtie2.yaml"
        log:   "logs/decont/{sample}.log"
        shell: r"""
            bowtie2 --very-sensitive --threads {threads} \
                    -x {params.prefix} -1 {input.fq1} -2 {input.fq2} \
                    > /dev/null 2> {log} \
            | samtools sort -n -@ {threads} \
            | samtools fastq -@ {threads} -f 12 -F 256 \
                  -1 {output.clean1} -2 {output.clean2} -
        """
