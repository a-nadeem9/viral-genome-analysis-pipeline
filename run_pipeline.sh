#!/bin/bash
#SBATCH --job-name=adv-pipeline
#SBATCH --output=logs/slurm/pipeline_%j.log
#SBATCH --error=logs/slurm/pipeline_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00

# ── Activate Snakemake environment ───────────────────────────────
eval "$(micromamba shell hook --shell bash)"
micromamba activate snakemake

# ── Create log directory ─────────────────────────────────────────
mkdir -p logs/slurm

# ── Dry run first (comment out once confirmed) ───────────────────
# snakemake -n --use-conda

# ── Run the full pipeline ────────────────────────────────────────
snakemake \
    --use-conda \
    --cores 32 \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds
