# -------- fastp trimming --------
rule fastp_trim:
    threads: 4
    input:
        fq1 = lambda wc: samples_df.set_index("sample").loc[wc.sample, "fq1"],
        fq2 = lambda wc: samples_df.set_index("sample").loc[wc.sample, "fq2"],
    output:
        trimmed1 = f"{QC_TRIM_DIR}/{{sample}}_1.trimmed.fastq.gz",
        trimmed2 = f"{QC_TRIM_DIR}/{{sample}}_2.trimmed.fastq.gz",
        html     = f"logs/fastp/{{sample}}_fastp.html",
        json     = f"logs/fastp/{{sample}}_fastp.json",
    params:
        adapter = config["adapter_path"],
        q  = config["parameters"]["fastp"]["qualified_quality_phred"],
        u  = config["parameters"]["fastp"]["unqualified_percent_limit"],
        n  = config["parameters"]["fastp"]["n_base_limit"],
        l  = config["parameters"]["fastp"]["length_required"],
        detect = "--detect_adapter_for_pe" if config["parameters"]["fastp"]["detect_adapter_for_pe"] else "",
    conda: "../envs/fastp.yaml"
    log:   "logs/fastp/{sample}.log"
    shell: r"""
        fastp \
          --in1  {input.fq1}  --in2  {input.fq2} \
          --out1 {output.trimmed1} --out2 {output.trimmed2} \
          --adapter_fasta {params.adapter} \
          --qualified_quality_phred {params.q} \
          --unqualified_percent_limit {params.u} \
          --n_base_limit {params.n} \
          --length_required {params.l} \
          {params.detect} \
          --thread {threads} \
          --html {output.html} --json {output.json} \
          2> {log}
    """
