#!/bin/bash
# ============================================================
# config.sh — Pipeline configuration
# Source this file in every SLURM script:  source ./config.sh
# Edit paths below to match your environment.
# ============================================================

# ── Project root ─────────────────────────────────────────────
export PROJECT_ROOT="${PROJECT_ROOT:-/path/to/project}"

# ── Reference databases ──────────────────────────────────────
export HOST_REF="${PROJECT_ROOT}/refs/Gigantidas_platifrons.mmi"   # minimap2 index of GCA_026163415.1
export GENOMAD_DB="${PROJECT_ROOT}/db/genomad_db"
export CHECKV_DB="${PROJECT_ROOT}/db/checkv-db-v1.5"
export VS2_DB="${PROJECT_ROOT}/db/vs2_db"
export IPHOP_DB="${PROJECT_ROOT}/db/iphop_db/Jun_2025_pub_rw"
export GTDBTK_DATA_PATH="${PROJECT_ROOT}/db/gtdbtk_r220"

# ── Working directories ─────────────────────────────────────
export READS_DIR="${PROJECT_ROOT}/data/clean_reads"        # *.filtered_1.fastq, *.filtered_2.fastq per sample
export ASSEMBLY_DIR="${PROJECT_ROOT}/work/assemblies"
export BINNING_DIR="${PROJECT_ROOT}/work/binning"
export VIRUS_DIR="${PROJECT_ROOT}/work/viruses"
export GLOBAL_DIR="${PROJECT_ROOT}/work/global"
export RESULTS_DIR="${PROJECT_ROOT}/results"

# ── Samples ─────────────────────────────────────────────────
export SAMPLES_LIST="${PROJECT_ROOT}/samples.txt"          # one sample ID per line, 48 lines

# ── Conda environments (adjust to your installation) ────────
export CONDA_BASE="${CONDA_BASE:-$HOME/miniconda3}"

# Auto-create work directories
mkdir -p "$ASSEMBLY_DIR" "$BINNING_DIR" "$VIRUS_DIR" "$GLOBAL_DIR" "$RESULTS_DIR" logs
