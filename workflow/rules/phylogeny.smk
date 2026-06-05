# workflow/rules/phylogeny.smk
# ------------------------------------------------------------
# IQ-TREE + Toytree plotting, using the fixed folder layout
# ------------------------------------------------------------
# Globals:  MSA_DIR, PLOTS_DIR

rule run_iqtree:
    threads: config["parameters"]["iqtree"]["threads"]
    input:
        aln = f"{MSA_DIR}/aligned_scaffolds.fasta"
    output:
        tree = f"{MSA_DIR}/tree.nwk"
    params:
        out_prefix = f"{MSA_DIR}/tree"
    log:
        "logs/iqtree/iqtree.log"
    conda: "../envs/iqtree.yaml"
    shell:
        """
        iqtree2 -s {input.aln} -nt {threads} -m GTR+G -bb 1000 \
                -pre {params.out_prefix} -redo 2> {log}
        mv {params.out_prefix}.treefile {output.tree}
        """

rule plot_pretty:
    input:
        windowed = f"{VAR_DIR}/variability_windowed.txt",
        tree     = f"{MSA_DIR}/tree.nwk",
    output:
        variability_plot = f"{PLOTS_DIR}/variability_clean.png",
        tree_plot        = f"{PLOTS_DIR}/tree_clean.png",
    params:
        n_strains = lambda wildcards: len(SAMPLES),
    conda:
        "../envs/plots.yaml"
    script:
        "../scripts/make_plots.py"

