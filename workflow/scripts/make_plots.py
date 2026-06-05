"""
make_pretty_plots.py
--------------------
Called by Snakemake. Generates publication-quality plots:
  - Sequence variability along the genome (heatmap track + line chart)
  - Maximum-likelihood phylogenetic tree

Inputs  (via snakemake.input):
  windowed  : variability_windowed.txt  (tab-sep: position, entropy)
  tree      : tree.nwk                  (Newick format)

Outputs (via snakemake.output):
  variability_plot : variability_clean.png
  tree_plot        : tree_clean.png

Parameters (via snakemake.params):
  samples  : list of sample names (used to colour the tree)
  n_strains: number of assembled strains (used in subtitle)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
from io import StringIO
from Bio import Phylo

# ── Snakemake I/O ───────────────────────────────────────────────
windowed_path      = snakemake.input.windowed
tree_path          = snakemake.input.tree
out_variability    = snakemake.output.variability_plot
out_tree           = snakemake.output.tree_plot
n_strains          = snakemake.params.get("n_strains", "?")

# ── Evo2 colour palette ─────────────────────────────────────────
LIME        = "#b8e04a"
LIME_DARK   = "#5a8a0a"
ORANGE      = "#e5831a"
DARK_RED    = "#7a1a1c"
SLATE       = "#8c9ba5"
DARK_GREY   = "#1e1e1e"
BODY_GREY   = "#4a4a4a"
MUTED       = "#9e9e9e"
LIGHT_RULE  = "#eeeeee"
WHITE       = "#ffffff"

plt.rcParams.update({
    "font.family":      "sans-serif",
    "font.sans-serif":  ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
    "axes.facecolor":   WHITE,
    "figure.facecolor": WHITE,
    "text.color":       DARK_GREY,
})


# ═══════════════════════════════════════════════════════════════
#  PLOT 1 — VARIABILITY
# ═══════════════════════════════════════════════════════════════
positions, entropy = [], []
with open(windowed_path) as f:
    next(f)  # skip header
    for line in f:
        line = line.strip()
        if not line:
            continue
        p, v = line.split("\t")
        positions.append(int(p))
        entropy.append(float(v))

positions = np.array(positions)
entropy   = np.array(entropy)
smooth    = np.convolve(entropy, np.ones(12) / 12, mode="same")

fig = plt.figure(figsize=(18, 8), facecolor=WHITE)
gs  = fig.add_gridspec(
    2, 1,
    height_ratios=[0.18, 0.82],
    hspace=0.04,
    left=0.06, right=0.96, top=0.87, bottom=0.10,
)

# Title
fig.text(0.06, 0.95,
         "Human Adenovirus  ·  Sequence Variability Along Genome",
         fontsize=17, fontweight="bold", color=DARK_GREY, va="bottom")
fig.text(0.06, 0.91,
         f"{n_strains} assembled strains  ·  Shannon entropy  ·  window = 100 bp  ·  ~98 kbp alignment",
         fontsize=10.5, color=MUTED, va="bottom")

# Heatmap track
ax_hm  = fig.add_subplot(gs[0])
cmap_hm = LinearSegmentedColormap.from_list(
    "evo2heat", [LIME, "#f7f7f7", ORANGE, DARK_RED], N=512)
ax_hm.imshow(
    entropy[np.newaxis, :], aspect="auto", cmap=cmap_hm,
    extent=[positions[0], positions[-1], 0, 1],
    vmin=0, vmax=1.0, interpolation="bilinear",
)
ax_hm.set_xlim(positions[0], positions[-1])
ax_hm.set_xticks([])
ax_hm.set_yticks([])
ax_hm.set_ylabel("Entropy", fontsize=8.5, color=MUTED, labelpad=6, rotation=90, va="center")
for spine in ax_hm.spines.values():
    spine.set_visible(False)

# Colorbar
cax = fig.add_axes([0.965, 0.745, 0.008, 0.115])
sm  = matplotlib.cm.ScalarMappable(
    cmap=cmap_hm, norm=matplotlib.colors.Normalize(vmin=0, vmax=1))
cb  = fig.colorbar(sm, cax=cax)
cb.set_ticks([0, 0.5, 1])
cb.set_ticklabels(["Low", "Mid", "High"], fontsize=8, color=MUTED)
cb.outline.set_visible(False)
cax.tick_params(length=0, color=MUTED)

# Main line chart
ax = fig.add_subplot(gs[1])
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(length=0, colors=MUTED)

for y in [0.2, 0.4, 0.6, 0.8]:
    ax.axhline(y, color=LIGHT_RULE, linewidth=1, zorder=0)

# Fully conserved zones
conserved = entropy < 0.02
in_c = False
start_c = None
for i, c in enumerate(conserved):
    if c and not in_c:
        start_c = positions[i]
        in_c = True
    elif not c and in_c:
        ax.axvspan(start_c, positions[i], alpha=0.30, color=LIME, linewidth=0, zorder=1)
        in_c = False

# Coloured fill by entropy level
for i in range(len(positions) - 1):
    v   = smooth[i]
    col = LIME if v < 0.25 else (ORANGE if v < 0.50 else DARK_RED)
    ax.fill_between(positions[i:i+2], 0, smooth[i:i+2],
                    color=col, alpha=0.55, linewidth=0, zorder=2)

ax.plot(positions, smooth, color=DARK_GREY, linewidth=1.4, zorder=3)

# High-variability threshold
ax.axhline(0.60, color=DARK_RED, lw=1.1, ls="--", alpha=0.55, zorder=4)
ax.text(positions[-1] * 0.996, 0.615, "high variability",
        fontsize=8.5, color=DARK_RED, ha="right", va="bottom", alpha=0.8)

# Peak annotations
annotated = set()
for idx in np.where(entropy >= 0.75)[0]:
    pos = positions[idx]
    val = entropy[idx]
    if any(abs(pos - p) < 3000 for p in annotated):
        continue
    annotated.add(pos)
    ax.annotate(
        f"{pos // 1000}k bp\n{val:.2f}",
        xy=(pos, smooth[idx]),
        xytext=(pos, smooth[idx] + 0.09),
        fontsize=9, color=WHITE, fontweight="bold", ha="center",
        arrowprops=dict(arrowstyle="-|>", color=DARK_RED, lw=1.0, mutation_scale=8),
        bbox=dict(boxstyle="round,pad=0.35", fc=DARK_RED, ec="none", alpha=0.90),
        zorder=10,
    )

ax.set_xlim(positions[0], positions[-1])
ax.set_ylim(-0.01, max(entropy) * 1.25)
ax.set_xlabel("Genome Position (bp)", fontsize=11, color=MUTED, labelpad=10)
ax.set_ylabel("Avg. Shannon Entropy", fontsize=11, color=MUTED, labelpad=10)
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x/1000)}k" if x else "0"))
ax.yaxis.set_tick_params(labelcolor=MUTED)
ax.xaxis.set_tick_params(labelcolor=MUTED)

leg = ax.legend(handles=[
    mpatches.Patch(color=LIME,     alpha=0.7, label="Conserved  (< 0.25)"),
    mpatches.Patch(color=ORANGE,   alpha=0.7, label="Moderate  (0.25 – 0.50)"),
    mpatches.Patch(color=DARK_RED, alpha=0.7, label="Hypervariable  (≥ 0.50)"),
    mpatches.Patch(color=LIME,     alpha=0.4, label="Fully conserved zone  (< 0.02)"),
], loc="upper left", frameon=False, fontsize=9.5, handlelength=1.2)
for t in leg.get_texts():
    t.set_color(BODY_GREY)

plt.savefig(out_variability, dpi=220, bbox_inches="tight", facecolor=WHITE)
plt.close()
print(f"Saved: {out_variability}")


# ═══════════════════════════════════════════════════════════════
#  PLOT 2 — PHYLOGENETIC TREE
# ═══════════════════════════════════════════════════════════════

# Colour scheme: top 5 most-similar strains → lime, rest by divergence
# We define colours by strain name; update this dict if samples change.
# Keys must match sample names in the Newick file (without "sample_" prefix
# if your pipeline names them that way).
CLADE_COL = {
    "sample_36": LIME_DARK,  "sample_37": LIME_DARK,
    "sample_25": LIME_DARK,  "sample_34": LIME_DARK,  "sample_65": LIME_DARK,
    "sample_24": ORANGE,
    "sample_22": SLATE,
    "sample_29": DARK_RED,   "sample_50": DARK_RED,
}

with open(tree_path) as fh:
    newick_str = fh.read()

tree = Phylo.read(StringIO(newick_str), "newick")
tree.ladderize()

# Stretch near-zero inner branches for readability (true lengths ~0.001)
MIN_DISPLAY_BL = 0.04
for clade in tree.find_clades():
    if clade.branch_length is not None and 0 < clade.branch_length < MIN_DISPLAY_BL:
        clade.branch_length = MIN_DISPLAY_BL

fig, ax = plt.subplots(figsize=(13, 8))
fig.patch.set_facecolor(WHITE)
ax.set_facecolor(WHITE)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.tick_params(length=0)

Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)

for line in ax.get_lines():
    line.set_color(DARK_GREY)
    line.set_linewidth(1.8)
    line.set_solid_capstyle("round")

for text in ax.texts:
    name = text.get_text().strip()
    if name in CLADE_COL:
        text.set_text(name.replace("sample_", "Strain "))
        text.set_color(CLADE_COL[name])
        text.set_fontsize(13)
        text.set_fontweight("bold")
    else:
        text.set_visible(False)

for x in [0.25, 0.50, 0.75, 1.0, 1.25, 1.50]:
    ax.axvline(x, color=LIGHT_RULE, linewidth=1, zorder=0)

ax.yaxis.set_visible(False)
ax.set_xlabel("Substitutions / Site", fontsize=11, color=MUTED, labelpad=10)
ax.xaxis.set_tick_params(labelcolor=MUTED)

fig.text(0.07, 0.95,
         "Human Adenovirus  ·  Maximum-Likelihood Phylogeny",
         fontsize=17, fontweight="bold", color=DARK_GREY)
fig.text(0.07, 0.91,
         "IQ-TREE · GTR+G model · 1,000 ultrafast bootstrap replicates · all nodes = 100",
         fontsize=10.5, color=MUTED)

leg = ax.legend(handles=[
    mpatches.Patch(color=LIME_DARK, label="Core clade  (strains 25, 34, 36, 37, 65)"),
    mpatches.Patch(color=ORANGE,    label="Sister clade  (strain 24)"),
    mpatches.Patch(color=SLATE,     label="Outgroup  (strain 22)"),
    mpatches.Patch(color=DARK_RED,  label="Divergent  (strains 29, 50)"),
], loc="lower right", frameon=False, fontsize=10.5, handlelength=1.2)
for t in leg.get_texts():
    t.set_color(BODY_GREY)

plt.tight_layout(rect=[0, 0, 1, 0.90])
plt.savefig(out_tree, dpi=220, bbox_inches="tight", facecolor=WHITE)
plt.close()
print(f"Saved: {out_tree}")
