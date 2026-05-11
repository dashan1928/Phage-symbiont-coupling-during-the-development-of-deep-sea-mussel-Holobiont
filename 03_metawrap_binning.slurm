#!/bin/bash
# ============================================================
# 03_metawrap_binning.slurm
# Three-tool binning + bin refinement with MetaWRAP
# Input :  $ASSEMBLY_DIR/Assembly_<SAMPLE>/final.contigs.fa
#          $READS_DIR/<SAMPLE>.filtered_{1,2}.fastq
# Output:  $BINNING_DIR/BinRefine_<SAMPLE>/metawrap_40_10_bins/
# ============================================================
#SBATCH --job-name=MetaWrap_Bin
#SBATCH --output=logs/binning_%A_%a.out
#SBATCH --error=logs/binning_%A_%a.err
#SBATCH --partition=normal,high_memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=1-48

set -euo pipefail
source ./config.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate metawrap

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")
CONTIGS="${ASSEMBLY_DIR}/Assembly_${SAMPLE}/final.contigs.fa"
READS1="${READS_DIR}/${SAMPLE}.filtered_1.fastq"
READS2="${READS_DIR}/${SAMPLE}.filtered_2.fastq"

BIN_OUT="${BINNING_DIR}/Binning_${SAMPLE}"
REFINE_OUT="${BINNING_DIR}/BinRefine_${SAMPLE}"

echo "[INFO] $(date) Binning $SAMPLE with MetaBAT2 + MaxBin2 + CONCOCT"
metawrap binning \
    -o "$BIN_OUT" \
    -t 36 \
    -a "$CONTIGS" \
    --metabat2 --maxbin2 --concoct \
    "$READS1" "$READS2"

echo "[INFO] $(date) Refining bins (completeness >= 40%, contamination <= 10%)"
metawrap bin_refinement \
    -o "$REFINE_OUT" \
    -t 32 \
    -A "${BIN_OUT}/metabat2_bins/" \
    -B "${BIN_OUT}/maxbin2_bins/" \
    -C "${BIN_OUT}/concoct_bins/" \
    -c 40 -x 10

echo "[INFO] $(date) Validating refined bins with CheckM"
checkm lineage_wf -t 32 -x fa \
    "${REFINE_OUT}/metawrap_40_10_bins/" \
    "${BINNING_DIR}/CheckM_${SAMPLE}"

echo "[INFO] $(date) $SAMPLE binning + refinement complete"
conda deactivate
