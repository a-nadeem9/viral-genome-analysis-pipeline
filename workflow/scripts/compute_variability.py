# scripts/compute_variability.py

from Bio import AlignIO
import math
import numpy as np
import matplotlib.pyplot as plt

# Snakemake-provided paths
alignment_file  = snakemake.input.msa
output_raw      = snakemake.output.txt
output_windowed = snakemake.output.windowed
output_plot     = snakemake.output.plot

# Configurable window size
window_size = snakemake.config["parameters"]["variability"]["window_size"]

def shannon_entropy(column):
    counts = {}
    for base in column:
        if base != "-":
            counts[base] = counts.get(base, 0) + 1
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)
    return entropy

# Load alignment
alignment = AlignIO.read(alignment_file, "fasta")
length    = alignment.get_alignment_length()

# Raw per-position entropy
entropies = [shannon_entropy(alignment[:, i]) for i in range(length)]
with open(output_raw, "w") as f:
    for val in entropies:
        f.write(f"{val:.5f}\n")

# Smoothed (windowed) entropy
smoothed = [
    np.mean(entropies[i : i + window_size])
    for i in range(0, length - window_size + 1, window_size)
]
with open(output_windowed, "w") as f:
    for idx, val in enumerate(smoothed):
        pos = idx * window_size
        f.write(f"{pos}\t{val:.5f}\n")

# Plot
plt.figure(figsize=(14, 5))
plt.plot(
    range(0, length - window_size + 1, window_size),
    smoothed,
    lw=1.2
)
plt.title("Sequence Variability Across Genome")
plt.xlabel("Alignment Position")
plt.ylabel("Shannon Entropy (Smoothed)")
plt.tight_layout()
plt.savefig(output_plot)
