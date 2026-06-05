# workflow/rules/reference_pipeline.smk
# ------------------------------------------------------------
# Reference selection → patch → mapping → samtools consensus
# Enabled when enable_reference_pipeline: true AND candidates_fasta is set.
# ------------------------------------------------------------

PIPE_ON = (
    str(config.get("enable_reference_pipeline", True)).lower() in {"true", "1", "yes"}
    and config.get("candidates_fasta", "") != ""
)

if PIPE_ON:

    CAND_FASTA = config["candidates_fasta"]

    # ---------- Mash sketch of candidate references ----------
    rule mash_sketch_candidates:
        input:
            fasta = CAND_FASTA
        output:
            sketch = "results/mash/candidates.msh"
        conda: "../envs/mash.yaml"
        log:   "logs/mash/candidates.log"
        shell: r"""
            mash sketch -o results/mash/candidates {input.fasta} > {log} 2>&1
        """

    # ---------- Mash sketch of each assembly ----------
    rule mash_sketch_sample:
        input:
            fasta = lambda wc: get_final_assembly(wc)
        output:
            sketch = temp(f"{ASSEMBLY_DIR}/{{sample}}/mash.msh")
        conda: "../envs/mash.yaml"
        log:   "logs/mash/{sample}.log"
        shell: "mash sketch -o {output.sketch} {input.fasta} > {log} 2>&1"

    # ---------- Pick best reference ----------
    rule select_best_reference:
        input:
            sample_sk = rules.mash_sketch_sample.output.sketch,
            cand_sk   = rules.mash_sketch_candidates.output.sketch,
        output:
            best_ref = f"{REFSEL_DIR}/{{sample}}_best_ref.txt"
        conda: "../envs/mash.yaml"
        log:   "logs/mash/{sample}_select.log"
        shell: r"""
            mash dist {input.sample_sk} {input.cand_sk} \
              | sort -k3,3n | head -n1 | cut -f2 | cut -d':' -f1 > {output.best_ref}
        """

    # ---------- RagTag patch: target (assembly) first, query (ref) second ----------
    rule ragtag_patch:
        threads: config["parameters"]["ragtag"]["threads"]
        input:
            contigs = lambda wc: get_final_assembly(wc),
            ref_txt = rules.select_best_reference.output.best_ref,
        output:
            patched = f"{PATCH_DIR}/{{sample}}/patched.fasta"
        params:
            outdir = f"{PATCH_DIR}/{{sample}}"
        conda: "../envs/ragtag.yaml"
        log:   "logs/ragtag/{sample}_patch.log"
        shell: r"""
            ragtag.py patch {input.contigs} $(cat {input.ref_txt}) \
                      -o {params.outdir} 2> {log} || true
            if [[ -f "{params.outdir}/ragtag.patch.fasta" ]]; then
                cp {params.outdir}/ragtag.patch.fasta {output.patched}
            else
                echo "RagTag found no alignments – using assembly as-is." >&2
                cp {input.contigs} {output.patched}
            fi
        """

    # ---------- Bowtie2 index of patched assembly ----------
    rule bowtie2_index_patched:
        input:
            fasta = rules.ragtag_patch.output.patched
        output:
            multiext(f"{PATCH_DIR}/{{sample}}/bt2_idx",
                     ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
        params:
            prefix = f"{PATCH_DIR}/{{sample}}/bt2_idx"
        conda: "../envs/bowtie2.yaml"
        log:   "logs/bwa/{sample}_index.log"
        shell: "bowtie2-build {input.fasta} {params.prefix} > /dev/null 2> {log}"

    # ---------- Map reads to patched assembly ----------
    READ_DIR = QC_DECONT_DIR if config["contaminants_fasta"] else QC_TRIM_DIR
    READ_TAG = "clean"       if config["contaminants_fasta"] else "trimmed"

    rule map_reads_for_consensus:
        threads: 8
        input:
            r1  = lambda wc: get_clean_reads(wc, 1),
            r2  = lambda wc: get_clean_reads(wc, 2),
            ref = rules.ragtag_patch.output.patched,
            idx = multiext(f"{PATCH_DIR}/{{sample}}/bt2_idx",
                           ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"),
        output:
            bam = temp(f"{CONS_DIR}/{{sample}}/aln.sorted.bam"),
        params:
            prefix = f"{PATCH_DIR}/{{sample}}/bt2_idx"
        conda: "../envs/bowtie2.yaml"
        log:   "logs/consensus/{sample}_map.log"
        shell: r"""
            bowtie2 --very-sensitive --threads {threads} \
                    -x {params.prefix} \
                    -1 {input.r1} -2 {input.r2} 2> {log} \
            | samtools sort -@ {threads} -o {output.bam}
            samtools index {output.bam}
        """

    # ---------- consensus via mpileup + python script (no REF BIAS) ----------
    rule call_consensus:
        input:
            bam = rules.map_reads_for_consensus.output.bam,
            ref = rules.ragtag_patch.output.patched,
        output:
            cons = f"{CONS_DIR}/{{sample}}/consensus.fasta",
        params:
            min_depth = config.get("consensus_min_depth", 5),
            script    = "workflow/scripts/mpileup_consensus.py",
        conda: "../envs/bowtie2.yaml"
        log:   "logs/consensus/{sample}.log"
        shell: r"""
            samtools faidx {input.ref}
            samtools mpileup -aa -q 0 -Q 0 -f {input.ref} {input.bam} \
              | python {params.script} \
                  --min-depth {params.min_depth} \
                  --sample {wildcards.sample} \
              > {output.cons} 2> {log}
        """

else:
    rule call_consensus:
        shell: "echo 'Reference pipeline disabled – consensus step skipped.'"
