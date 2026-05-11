#!/bin/bash
# ============================================================
# 08_abundance_quantification.slurm
# Per-sample competitive mapping of reads against
# (vOTU catalogue + dereplicated bacterial MAGs) with CoverM
#
# Output: $RESULTS_DIR/CoverM/<SAMPLE>_abundance.tsv
#         (TPM, raw count, covered_fraction columns)
# ============================================================
#SBATCH --job-name=CoverM
#SBATCH --output=logs/coverm_%A_%a.out
#SBATCH --error=logs/coverm_%A_%a.err
#SBATCH --partition=normal,high_memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-48

set -euo pipefail
source ./config.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate coverm

VIRUS_FA="${RESULTS_DIR}/iPHoP_Output_Network/Global_Dereplicated_vOTUs_clean.fna"
HOST_FA="${GLOBAL_DIR}/Global_Refined_Bacterial_MAGs.fa"
COMBINED_FA="${GLOBAL_DIR}/Holobiont_Combined_Reference.fa"
COVERM_OUT="${RESULTS_DIR}/CoverM"
mkdir -p "$COVERM_OUT"

# Build combined reference once (only the first array task creates it)
# Use a lockfile to prevent race conditions across array tasks
LOCKFILE="${COMBINED_FA}.lock"
if [ ! -f "$COMBINED_FA" ]; then
    (
        flock -n 9 || { echo "[INFO] another task is building reference; waiting..."; flock 9; }
        if [ ! -f "$COMBINED_FA" ]; then
            echo "[INFO] $(date) Building competitive reference (vOTUs + MAGs)"
            cat "$VIRUS_FA" "$HOST_FA" > "${COMBINED_FA}.tmp"
            mv "${COMBINED_FA}.tmp" "$COMBINED_FA"
        fi
    ) 9>"$LOCKFILE"
fi

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")
READS1="${READS_DIR}/${SAMPLE}.filtered_1.fastq"
READS2="${READS_DIR}/${SAMPLE}.filtered_2.fastq"

echo "[INFO] $(date) CoverM for $SAMPLE"
coverm contig \
    -1 "$READS1" -2 "$READS2" \
    --reference "$COMBINED_FA" \
    --methods tpm count covered_fraction \
    --min-read-percent-identity 0.95 \
    --min-read-aligned-percent 0.75 \
    --output-file "${COVERM_OUT}/${SAMPLE}_abundance.tsv" \
    -t 16

echo "[INFO] $(date) $SAMPLE abundance complete"
conda deactivate
