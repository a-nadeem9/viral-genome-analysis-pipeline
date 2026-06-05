# workflow/rules/assembly.smk
# ------------------------------------------------------------
# SPAdes assembly + RagTag scaffolding
# Uses get_clean_reads() from common.smk
# ------------------------------------------------------------

rule spades_assemble:
    threads: config["parameters"]["spades"]["threads"]
    input:
        r1 = lambda wc: get_clean_reads(wc, 1),
        r2 = lambda wc: get_clean_reads(wc, 2),
    output:
        contigs = f"{ASSEMBLY_DIR}/{{sample}}/contigs/contigs.fasta",
    params:
        mem    = config["parameters"]["spades"]["memory"],
        outdir = f"{ASSEMBLY_DIR}/{{sample}}/contigs",
    conda: "../envs/spades.yaml"
    log:   "logs/spades/{sample}.log"
    shell: r"""
        spades.py \
          --phred-offset 33 \
          -1 {input.r1} -2 {input.r2} \
          -o {params.outdir} \
          -t {threads} -m {params.mem} \
          > {log} 2>&1
    """

rule ragtag_scaffold:
    threads: config["parameters"]["ragtag"]["threads"]
    input:
        contigs = rules.spades_assemble.output.contigs,
        ref     = config["reference"],
    output:
        scaffolds = f"{ASSEMBLY_DIR}/{{sample}}/scaffolds.fasta",
    params:
        outdir = f"{ASSEMBLY_DIR}/{{sample}}/ragtag_out",
    conda: "../envs/ragtag.yaml"
    log:   "logs/ragtag/{sample}.log"
    shell: r"""
        ragtag.py scaffold {input.ref} {input.contigs} \
               -o {params.outdir} 2> {log} || true
        if [[ -f "{params.outdir}/ragtag.scaffold.fasta" ]]; then
            mv {params.outdir}/ragtag.scaffold.fasta {output.scaffolds}
        else
            echo "WARNING: RagTag scaffold failed (no useful alignments). Using raw contigs." >> {log}
            cp {input.contigs} {output.scaffolds}
        fi
    """

rule extract_ragtag_scaffold:
    """Extract the RagTag-scaffolded sequence (suffix _RagTag) and rename header.
    Falls back to longest contig if no _RagTag scaffold exists (divergent samples)."""
    input:
        scaff = rules.ragtag_scaffold.output.scaffolds
    output:
        merged = f"{ASSEMBLY_DIR}/{{sample}}/merged/merged.fasta"
    conda:  "../envs/seqkit.yaml"
    shell:  r"""
        MATCH_COUNT=$(grep -c "_RagTag" {input.scaff} || true)
        if [ "$MATCH_COUNT" -gt 0 ]; then
            seqkit grep -r -p "_RagTag$" {input.scaff} \
            | sed 's/^>.*/>{wildcards.sample}/' > {output.merged}
        else
            # No RagTag scaffold — sample too divergent; use longest contig
            seqkit sort --by-length --reverse {input.scaff} -o {output.merged}.tmp
            awk 'NR==1{{print ">{wildcards.sample}"; next}} /^>/{{exit}} {{print}}' {output.merged}.tmp > {output.merged}
            rm -f {output.merged}.tmp
        fi
    """
