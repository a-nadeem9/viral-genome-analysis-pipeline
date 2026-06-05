# Adenovirus Genome Analysis Pipeline

A complete Snakemake pipeline for Human Adenovirus genome assembly, comparative genomics, and phylogenetic analysis.

---

## What does this pipeline do?

You start with raw paired-end sequencing reads (FASTQ files) from one or more adenovirus samples. The pipeline takes them all the way from raw reads to publication-quality figures showing how different your virus strains are from each other.

### Step-by-step overview

```
Raw reads (FASTQ)
     │
     ▼
① Quality control (fastp)
     │  Trim adapters, remove low-quality bases and short reads
     │
     ▼
② [Optional] Host decontamination (Bowtie2)
     │  Remove any human/host DNA that ended up in your sample
     │
     ├──► [Optional] Taxonomic screening (Kraken2)
     │         "What organisms are in this sample?"
     │
     ▼
③ De novo genome assembly (SPAdes)
     │  Puzzle the cleaned reads together into genome fragments (contigs)
     │
     ▼
④ Scaffolding (RagTag scaffold)
     │  Order and orient contigs using a reference genome
     │
     ▼
⑤ [Optional] Assembly polishing (Pilon)
     │  Map reads back to the assembly and fix remaining errors
     │
     ▼
⑥ Reference selection (Mash)
     │  For each sample, find the most similar known adenovirus strain
     │
     ▼
⑦ Gap filling / patching (RagTag patch)
     │  Fill assembly gaps using the best-matching reference
     │
     ▼
⑧ Consensus calling (BWA-MEM2 + samtools consensus)
     │  Map reads to the patched assembly → build a final consensus sequence
     │  Positions with low read coverage get 'N' (unknown) — no reference bias
     │
     ▼
⑨ Multiple sequence alignment (MAFFT)
     │  Align all consensus sequences side by side for comparison
     │
     ├──► ⑩ Variability analysis
     │         How variable is each genomic position across all samples?
     │         Uses Shannon entropy. Outputs a windowed plot.
     │
     └──► ⑪ Phylogenetic tree (IQ-TREE, GTR+G model)
               Which samples are most closely related?
               Bootstrapped (1000 replicates) for confidence.
                    │
                    ▼
              ⑫ Publication-quality plots
                   • Genome variability plot (variability_clean.png)
                   • Phylogenetic tree (tree_clean.png)

⑬ QC summary report (MultiQC)
     Aggregates all per-sample fastp statistics into one HTML report
```

---

## Requirements

- [Micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html) or Conda
- Snakemake ≥ 9

Create the Snakemake environment:
```bash
micromamba create -n snakemake -c conda-forge snakemake
micromamba activate snakemake
```

All other tools (SPAdes, RagTag, Pilon, MAFFT, IQ-TREE, etc.) are installed automatically by Snakemake into isolated conda environments when you first run the pipeline.

---

## Project structure

```
adenovirus-genome-analysis-pipeline/
├── config/
│   ├── config.yaml        # All settings — edit this before running
│   └── samples.tsv        # Your sample names and read file paths
├── data/
│   ├── reads/             # Place your FASTQ files here
│   └── reference/
│       └── ref.fasta      # Scaffolding reference genome
├── resources/
│   └── adapters/          # Sequencing adapter sequences (for fastp)
├── run_pipeline.sh        # SLURM submission script
└── workflow/
    ├── Snakefile           # Main pipeline definition
    ├── envs/               # Conda environment files (one per tool)
    ├── rules/              # Modular rule files (one per analysis step)
    └── scripts/            # Custom Python scripts
        ├── make_plots.py
        └── compute_variability.py
```

---

## Setup

### 1. Add your samples

Put your FASTQ files in `data/reads/` and edit `config/samples.tsv`:

```
sample      fq1                               fq2
sample_22   data/reads/sample_22_1.fastq.gz   data/reads/sample_22_2.fastq.gz
sample_24   data/reads/sample_24_1.fastq.gz   data/reads/sample_24_2.fastq.gz
```

### 2. Configure the pipeline

Edit `config/config.yaml`. The most important settings:

```yaml
# --- Required ---
reference:    "data/reference/ref.fasta"      # reference genome for scaffolding
samples_tsv:  "config/samples.tsv"

# --- Optional steps (set to "" to disable) ---
contaminants_fasta: ""     # path to host genome FASTA, or "" to skip decontamination
kraken2_db:         ""     # path to Kraken2 database, or "" to skip

# --- Reference-guided consensus (recommended: true) ---
enable_polishing:           true
enable_reference_pipeline:  true
candidates_fasta: "resources/reference/candidates.fasta"

# --- Consensus quality threshold ---
consensus_min_depth: 5     # positions with fewer reads than this get 'N'
```

### 3. Add candidate reference genomes (for consensus calling)

Place a FASTA file containing representative adenovirus reference genomes at `resources/reference/candidates.fasta`. The pipeline will automatically select the best-matching reference for each sample using Mash distance.
---

## Outputs

All results are written to the `results/` directory.

| Path | Description |
|---|---|
| `results/qc/trimmed/` | Quality-trimmed FASTQ files |
| `results/assembly/{sample}/merged/merged.fasta` | Scaffolded assembly |
| `results/assembly/{sample}/pilon/polished.fasta` | Polished assembly |
| `results/consensus/{sample}/consensus.fasta` | Final consensus sequence |
| `results/msa/aligned_scaffolds.fasta` | Multiple sequence alignment |
| `results/msa/tree.nwk` | Phylogenetic tree (Newick format) |
| `results/variability/variability_windowed.txt` | Per-window variability scores |
| `results/plots/variability_clean.png` | 🎨 Genome variability figure |
| `results/plots/tree_clean.png` | 🎨 Phylogenetic tree figure |
| `results/multiqc/fastp_report.html` | 📊 QC summary report |
| `results/kraken/multiqc/kraken_report.html` | 📊 Taxonomic screening report (if enabled) |

---

## Reference Genomes used for candidates.fasta

The `candidates.fasta` file contains full genomes from 8 Human Adenovirus strains downloaded from NCBI RefSeq, representing diverse adenovirus species:

- Human Adenovirus A
- Human Adenovirus B
- Human Adenovirus C
- Human Adenovirus D
- Human Adenovirus E
- Human Adenovirus F
- Human Adenovirus 1
- Human Adenovirus 2

---