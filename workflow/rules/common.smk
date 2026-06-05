# workflow/rules/common.smk
# ------------------------------------------------------------
# Shared input functions used across multiple rules.
# Included first in the Snakefile so all rules can access them.
# ------------------------------------------------------------


def get_clean_reads(wc, read_num):
    """Return the correct input reads — decontaminated if enabled, else trimmed."""
    if config["contaminants_fasta"]:
        return f"{QC_DECONT_DIR}/{wc.sample}_{read_num}.clean.fastq.gz"
    return f"{QC_TRIM_DIR}/{wc.sample}_{read_num}.trimmed.fastq.gz"


def get_final_assembly(wc):
    """Return the polished assembly if polishing is enabled, else the merged scaffold."""
    if str(config.get("enable_polishing", True)).lower() in {"true", "yes", "1"}:
        return f"{ASSEMBLY_DIR}/{wc.sample}/pilon/polished.fasta"
    return f"{ASSEMBLY_DIR}/{wc.sample}/merged/merged.fasta"


def get_msa_input(wc):
    """Return consensus FASTA if reference pipeline is on, else polished/merged assembly."""
    if (
        str(config.get("enable_reference_pipeline", True)).lower() in {"true", "yes", "1"}
        and config.get("candidates_fasta", "")
    ):
        return f"{CONS_DIR}/{wc.sample}/consensus.fasta"
    return get_final_assembly(wc)
