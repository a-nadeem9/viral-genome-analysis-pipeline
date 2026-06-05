# workflow/rules/assembly_stats.smk
# ------------------------------------------------------------
# Per-sample assembly statistics using seqkit stats
# Runs on the final polished (or merged) assembly per sample
# ------------------------------------------------------------

rule compute_assembly_stats:
    input:
        fasta = lambda wc: get_final_assembly(wc)
    output:
        tsv = "results/assembly_stats/{sample}_stats.tsv"
    conda: "../envs/seqkit.yaml"
    log:   "logs/assembly_stats/{sample}.log"
    shell: r"""
        seqkit stats -a -T {input.fasta} > {output.tsv} 2> {log}
    """

rule aggregate_assembly_stats:
    input:
        tsvs = expand("results/assembly_stats/{sample}_stats.tsv", sample=SAMPLES)
    output:
        combined = "results/assembly_stats/all_samples_stats.tsv"
    conda: "../envs/plots.yaml"
    script: "../scripts/aggregate_assembly_stats.py"


rule plot_assembly_stats:
    input:
        tsv = rules.aggregate_assembly_stats.output.combined
    output:
        plot = "results/plots/assembly_stats.png"
    conda: "../envs/plots.yaml"
    script: "../scripts/plot_assembly_stats.py"
