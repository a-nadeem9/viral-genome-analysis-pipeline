"""
workflow/scripts/plot_assembly_stats.py
---------------------------------------
Snakemake script: plots assembly statistics per sample.
Style matches the variability and tree plots (light theme).

Inputs:
    snakemake.input.tsv  — aggregated seqkit stats TSV
Outputs:
    snakemake.output.plot — PNG figure
"""

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np

# ── font ───────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["DejaVu Sans", "Arial", "Helvetica"],
    "axes.spines.top":  False,
    "axes.spines.right":False,
})

# ── load data ──────────────────────────────────────────────────────────────
df = pd.read_csv(snakemake.input.tsv, sep="\t")
df["sample"] = df["sample"].str.replace("sample_", "Strain ", regex=False)
df = df.sort_values("N50", ascending=False).reset_index(drop=True)

samples  = df["sample"].tolist()
n50      = df["N50"].tolist()
sum_len  = df["sum_len"].tolist()
num_seqs = df["num_seqs"].tolist()
x        = np.arange(len(samples))

# ── colours (matching the variability plot palette) ────────────────────────
COL_GREEN  = "#6a994e"   # genome size — conserved green
COL_ORANGE = "#e07b39"   # N50 — moderate orange
COL_RED    = "#8b3a3a"   # contigs — hypervariable dark red
BG         = "#f7f7f5"   # very light panel background
GRID       = "#dddddd"

# ── figure ─────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.patch.set_facecolor("#ffffff")

def draw_panel(ax, values, colour, title, ylabel, fmt_top=None):
    ax.set_facecolor(BG)
    bars = ax.bar(x, values, width=0.6, color=colour,
                  edgecolor="white", linewidth=1.0, zorder=3)
    # value labels on top
    for bar, val in zip(bars, values):
        label = fmt_top(val) if fmt_top else str(int(val))
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(values) * 0.015,
            label, ha="center", va="bottom",
            fontsize=8.5, color="#333333", fontweight="bold",
        )
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=35, ha="right", fontsize=9)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=8, color="#111111")
    ax.set_ylabel(ylabel, fontsize=9.5, color="#444444")
    ax.tick_params(colors="#444444")
    ax.grid(axis="y", linestyle="--", linewidth=0.7, color=GRID, zorder=0)
    ax.set_axisbelow(True)
    ax.spines["left"].set_color("#bbbbbb")
    ax.spines["bottom"].set_color("#bbbbbb")

def kbp(v, _=None):
    return f"{v/1000:.0f}k" if v >= 1000 else str(int(v))

draw_panel(axes[0], sum_len,  COL_GREEN,
           "Assembled Genome Size", "Total Length (bp)",
           fmt_top=lambda v: f"{v/1000:.1f}k")
axes[0].yaxis.set_major_formatter(mticker.FuncFormatter(kbp))

draw_panel(axes[1], n50, COL_ORANGE,
           "Assembly N50", "N50 (bp)",
           fmt_top=lambda v: f"{v/1000:.1f}k")
axes[1].yaxis.set_major_formatter(mticker.FuncFormatter(kbp))

draw_panel(axes[2], num_seqs, COL_RED,
           "Number of Contigs", "Count")
axes[2].yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: str(int(v))))

# ── suptitle (same style as variability/tree plots) ────────────────────────
n = len(samples)
fig.suptitle(
    "Human Adenovirus  ·  Polished Assembly Statistics",
    fontsize=14, fontweight="bold", color="#111111", y=1.02,
)
fig.text(
    0.5, 0.985,
    f"{n} assembled strains  ·  seqkit stats  ·  polished assemblies",
    ha="center", fontsize=9.5, color="#777777",
)

plt.tight_layout()
fig.savefig(
    snakemake.output.plot,
    dpi=180, bbox_inches="tight",
    facecolor="#ffffff",
)
plt.close(fig)
