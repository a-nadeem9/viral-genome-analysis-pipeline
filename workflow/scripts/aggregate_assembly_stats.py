"""
workflow/scripts/aggregate_assembly_stats.py
---------------------------------------------
Combine per-sample seqkit stats TSV files into one table.
"""

import pandas as pd

dfs = []
for f in snakemake.input.tsvs:
    df = pd.read_csv(f, sep="\t")
    sample = f.split("/")[-1].replace("_stats.tsv", "")
    df["sample"] = sample
    dfs.append(df)

combined = pd.concat(dfs, ignore_index=True)
combined.to_csv(snakemake.output.combined, sep="\t", index=False)
