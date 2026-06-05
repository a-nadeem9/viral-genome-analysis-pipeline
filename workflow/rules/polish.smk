# workflow/rules/polish.smk
# ------------------------------------------------------------
# Optional short-read polishing with Pilon.
# Uses get_clean_reads() and get_final_assembly() from common.smk.
# finalize_assembly removed — downstream rules call get_final_assembly() directly.
# ------------------------------------------------------------

ENABLE_POLISH = str(config.get("enable_polishing", True)).lower() in {"true", "yes", "1"}

if ENABLE_POLISH:

    rule build_bowtie2_index_pilon:
        input:
            fasta = f"{ASSEMBLY_DIR}/{{sample}}/merged/merged.fasta"
        output:
            multiext(f"{ASSEMBLY_DIR}/{{sample}}/pilon/idx",
                     ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params:
            prefix = f"{ASSEMBLY_DIR}/{{sample}}/pilon/idx"
        conda: "../envs/bowtie2.yaml"
        log:   "logs/pilon/{sample}_index.log"
        shell: r"""
            bowtie2-build {input.fasta} {params.prefix} > /dev/null 2> {log}
        """

    rule map_reads_for_pilon:
        threads: 8
        input:
            r1  = lambda wc: get_clean_reads(wc, 1),
            r2  = lambda wc: get_clean_reads(wc, 2),
            idx = multiext(f"{ASSEMBLY_DIR}/{{sample}}/pilon/idx",
                           ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        output:
            bam = temp(f"{ASSEMBLY_DIR}/{{sample}}/pilon/aln.sorted.bam")
        params:
            prefix = f"{ASSEMBLY_DIR}/{{sample}}/pilon/idx"
        conda: "../envs/bowtie2.yaml"
        log:   "logs/pilon/{sample}_map.log"
        shell: r"""
            bowtie2 --very-sensitive --threads {threads} \
                    -x {params.prefix} \
                    -1 {input.r1} -2 {input.r2} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam}
            samtools index {output.bam}
        """

    rule run_pilon:
        threads: config["parameters"]["pilon"]["threads"]
        input:
            bam = rules.map_reads_for_pilon.output.bam,
            ref = f"{ASSEMBLY_DIR}/{{sample}}/merged/merged.fasta"
        output:
            polished = f"{ASSEMBLY_DIR}/{{sample}}/pilon/polished.fasta"
        params:
            outdir = f"{ASSEMBLY_DIR}/{{sample}}/pilon"
        conda: "../envs/pilon.yaml"
        log:   "logs/pilon/{sample}.log"
        shell: r"""
            pilon --genome {input.ref} --frags {input.bam} \
                  --output polished --outdir {params.outdir} 2> {log}
        """
