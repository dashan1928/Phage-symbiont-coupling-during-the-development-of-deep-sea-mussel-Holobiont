#!/bin/bash
# ============================================================
# 02_assembly.slurm
# Per-sample metagenomic assembly with MEGAHIT
# Input :  $READS_DIR/<SAMPLE>.filtered_{1,2}.fastq
# Output:  $ASSEMBLY_DIR/Assembly_<SAMPLE>/final.contigs.fa
# ============================================================
#SBATCH --job-name=MEGAHIT
#SBATCH --output=logs/assembly_%A_%a.out
#SBATCH --error=logs/assembly_%A_%a.err
#SBATCH --partition=normal,high_memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --array=1-48

set -euo pipefail
source ./config.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate megahit

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")
OUT_DIR="${ASSEMBLY_DIR}/Assembly_${SAMPLE}"

# MEGAHIT requires the output dir not to pre-exist
if [ -d "$OUT_DIR" ]; then
    echo "[WARN] $(date) Assembly already exists for $SAMPLE — skipping"
    exit 0
fi

echo "[INFO] $(date) Assembling $SAMPLE with MEGAHIT meta-large preset"
megahit \
    -1 "${READS_DIR}/${SAMPLE}.filtered_1.fastq" \
    -2 "${READS_DIR}/${SAMPLE}.filtered_2.fastq" \
    -o "$OUT_DIR" \
    --presets meta-large \
    -t 36 \
    --min-contig-len 1000 \
    --memory 0.85

echo "[INFO] $(date) $SAMPLE assembly complete"
conda deactivate
