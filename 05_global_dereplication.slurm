#!/bin/bash
# ============================================================
# 05_global_dereplication.slurm
# Pool & dereplicate MAGs with dRep, vOTUs with cd-hit-est
# Input :  per-sample BinRefine_<SAMPLE>/metawrap_40_10_bins/*.fa
#          per-sample verified_vOTUs_<SAMPLE>.fa
# Output:  $GLOBAL_DIR/Global_Refined_Bacterial_MAGs.fa
#          $GLOBAL_DIR/Global_Dereplicated_vOTUs.fasta
# ============================================================
#SBATCH --job-name=Global_Derep
#SBATCH --output=logs/derep_%j.out
#SBATCH --error=logs/derep_%j.err
#SBATCH --partition=normal,high_memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=256G
#SBATCH --time=48:00:00

set -euo pipefail
source ./config.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ── 1. Pool & rename MAGs from all samples ────────────────────
echo "[INFO] $(date) Stage 1: pooling MAGs from 48 samples"
mkdir -p "${GLOBAL_DIR}/Global_MAGs_Pool"

while read -r SAMPLE; do
    BIN_DIR="${BINNING_DIR}/BinRefine_${SAMPLE}/metawrap_40_10_bins"
    [ -d "$BIN_DIR" ] || { echo "[WARN] $SAMPLE: no refined bins"; continue; }

    for bin in "$BIN_DIR"/*.fa; do
        [ -f "$bin" ] || continue
        bin_name=$(basename "$bin" .fa)
        # Prefix every contig name with sample_bin to avoid collisions
        awk -v p="${SAMPLE}_${bin_name}_" \
            '/^>/{print ">" p substr($0, 2); next} {print}' \
            "$bin" > "${GLOBAL_DIR}/Global_MAGs_Pool/${SAMPLE}_${bin_name}.fa"
    done
done < "$SAMPLES_LIST"

# ── 2. dRep dereplication of MAGs ─────────────────────────────
echo "[INFO] $(date) Stage 2: dRep dereplication (completeness >= 50, contam <= 10)"
conda activate drep
dRep dereplicate "${GLOBAL_DIR}/dRep_out" \
    -g "${GLOBAL_DIR}/Global_MAGs_Pool/"*.fa \
    -p 48 \
    -comp 50 -con 10
conda deactivate

# Concatenate representative MAGs
cat "${GLOBAL_DIR}/dRep_out/dereplicated_genomes/"*.fa \
    > "${GLOBAL_DIR}/Global_Refined_Bacterial_MAGs.fa"
echo "[INFO] $(date) Bacterial MAGs catalogue: ${GLOBAL_DIR}/Global_Refined_Bacterial_MAGs.fa"

# ── 3. Pool & dereplicate vOTUs ───────────────────────────────
echo "[INFO] $(date) Stage 3: pooling and clustering vOTUs"
cat "${VIRUS_DIR}/verified_vOTUs_"*.fa > "${GLOBAL_DIR}/all_samples_verified_viruses.fa"

conda activate cdhit
cd-hit-est \
    -i "${GLOBAL_DIR}/all_samples_verified_viruses.fa" \
    -o "${GLOBAL_DIR}/Global_Dereplicated_vOTUs.fasta" \
    -c 0.95 -aS 0.85 \
    -T 48 -M 0 -d 0
conda deactivate

echo "[INFO] $(date) vOTU catalogue: ${GLOBAL_DIR}/Global_Dereplicated_vOTUs.fasta"
echo "[INFO] $(date) All dereplication complete"
