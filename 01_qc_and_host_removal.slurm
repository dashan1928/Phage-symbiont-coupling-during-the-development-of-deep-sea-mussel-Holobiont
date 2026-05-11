#!/bin/bash
# ============================================================
# 01_qc_and_host_removal.slurm
# QC trimming with fastp + host read removal with minimap2
# Input :  raw_reads/<SAMPLE>_1.fastq.gz, raw_reads/<SAMPLE>_2.fastq.gz
# Output:  $READS_DIR/<SAMPLE>.filtered_{1,2}.fastq
# ============================================================
#SBATCH --job-name=QC_HostRm
#SBATCH --output=logs/qc_%A_%a.out
#SBATCH --error=logs/qc_%A_%a.err
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=64G
#SBATCH --array=1-48

set -euo pipefail
source ./config.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")
RAW_R1="${PROJECT_ROOT}/data/raw_reads/${SAMPLE}_1.fastq.gz"
RAW_R2="${PROJECT_ROOT}/data/raw_reads/${SAMPLE}_2.fastq.gz"
QC_DIR="${PROJECT_ROOT}/work/qc"
mkdir -p "$QC_DIR" "$READS_DIR"

echo "[INFO] $(date) Sample: $SAMPLE — fastp QC trimming"
conda activate fastp
fastp \
    -i "$RAW_R1" -I "$RAW_R2" \
    -o "${QC_DIR}/${SAMPLE}.qc_1.fastq.gz" \
    -O "${QC_DIR}/${SAMPLE}.qc_2.fastq.gz" \
    --json "${QC_DIR}/${SAMPLE}.fastp.json" \
    --html "${QC_DIR}/${SAMPLE}.fastp.html" \
    --length_required 50 \
    --thread 16
conda deactivate

echo "[INFO] $(date) Sample: $SAMPLE — host read removal (minimap2 vs $HOST_REF)"
conda activate minimap2

# Strategy: align reads to host index. Reads with edit distance NM < 15
# (i.e. high identity to host) are flagged as host-derived and removed.
# Reads with NM >= 15 OR no NM tag (e.g. unmapped reads with secondary/none) are kept.
# Note: this script extracts the IDs of host-derived reads, then filters them out.

minimap2 -t 32 -N 50 -ax sr "$HOST_REF" \
    "${QC_DIR}/${SAMPLE}.qc_1.fastq.gz" \
    "${QC_DIR}/${SAMPLE}.qc_2.fastq.gz" \
  | samtools view -F 4 - \
  | awk 'BEGIN{OFS="\t"} {
        nm = -1
        for (i=12; i<=NF; i++) if ($i ~ /^NM:i:/) { split($i,a,":"); nm=a[3]; break }
        if (nm >= 0 && nm < 15) print $1
    }' \
  | sort -u > "${QC_DIR}/${SAMPLE}.host_ids.list"

# Use seqkit to remove host-derived read pairs
seqkit grep -v -f "${QC_DIR}/${SAMPLE}.host_ids.list" \
    "${QC_DIR}/${SAMPLE}.qc_1.fastq.gz" \
    -o "${READS_DIR}/${SAMPLE}.filtered_1.fastq"
seqkit grep -v -f "${QC_DIR}/${SAMPLE}.host_ids.list" \
    "${QC_DIR}/${SAMPLE}.qc_2.fastq.gz" \
    -o "${READS_DIR}/${SAMPLE}.filtered_2.fastq"

conda deactivate
echo "[INFO] $(date) Sample: $SAMPLE — completed"
