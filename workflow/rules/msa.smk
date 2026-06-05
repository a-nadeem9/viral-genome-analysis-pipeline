# workflow/rules/msa.smk
# ------------------------------------------------------------
# Multi-sequence alignment of per-sample consensus sequences.
# Uses get_msa_input() from common.smk:
#   - consensus FASTA if reference pipeline is enabled
#   - polished/merged assembly otherwise
# ------------------------------------------------------------

rule concat_scaffolds:
    input:
        expand("{cons_dir}/{sample}/consensus.fasta" if (
            str(config.get("enable_reference_pipeline", True)).lower() in {"true", "yes", "1"}
            and config.get("candidates_fasta", "")
        ) else "{assembly_dir}/{sample}/merged/merged.fasta",
               sample=SAMPLES,
               cons_dir=CONS_DIR,
               assembly_dir=ASSEMBLY_DIR)
    output:
        combined = f"{MSA_DIR}/scaffolds_all.fasta"
    shell:
        "cat {input} > {output.combined}"

rule run_mafft:
    threads: config["parameters"]["mafft"]["threads"]
    input:
        combined = rules.concat_scaffolds.output.combined
    output:
        aligned = f"{MSA_DIR}/aligned_scaffolds.fasta"
    log:
        "logs/mafft/mafft.log"
    conda: "../envs/mafft.yaml"
    shell:
        "mafft --thread {threads} {input.combined} > {output.aligned} 2> {log}"
