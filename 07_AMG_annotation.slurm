#!/bin/bash
# ============================================================
# 07_AMG_annotation.slurm
# Two complementary AMG callers:
#   (a) VIBRANT  — virome-optimised lifecycle + AMG annotation
#   (b) DRAM-v   — high-resolution KEGG/Pfam/UniRef AMG annotation
# Input:  $RESULTS_DIR/iPHoP_Output_Network/Global_Dereplicated_vOTUs_clean.fna
# ============================================================
#SBATCH --job-name=AMG
#SBATCH --output=logs/amg_%j.out
#SBATCH --error=logs/amg_%j.err
#SBATCH --partition=normal,high_memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=128G
#SBATCH --time=48:00:00

set -euo pipefail
source ./config.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"

VIRUS_FA="${RESULTS_DIR}/iPHoP_Output_Network/Global_Dereplicated_vOTUs_clean.fna"
AMG_DIR="${RESULTS_DIR}/AMG"
mkdir -p "$AMG_DIR"
cd "$AMG_DIR"

# ── (a) VIBRANT ───────────────────────────────────────────────
echo "[INFO] $(date) Running VIBRANT (-virome mode)"
conda activate vibrant
VIBRANT_run.py \
    -i "$VIRUS_FA" \
    -t 48 \
    -virome \
    -folder "${AMG_DIR}/VIBRANT"
conda deactivate
echo "[INFO] $(date) VIBRANT outputs in ${AMG_DIR}/VIBRANT"

# ── (b) DRAM-v ────────────────────────────────────────────────
# Note: DRAM-v ideally uses a VirSorter2 affi-contigs.tab file when available.
# If you ran VirSorter2 on the dereplicated vOTUs, supply that here.
# Otherwise DRAM-v can still annotate but flags will be less informative.
echo "[INFO] $(date) Running DRAM-v annotation"
conda activate dramv

# Generate VirSorter2 input for DRAM-v (recommended by DRAM-v docs)
virsorter run --prep-for-dramv \
    -i "$VIRUS_FA" \
    -w "${AMG_DIR}/vs2_for_dramv" \
    --include-groups dsDNAphage,ssDNA \
    --min-score 0.5 \
    -j 48 \
    --db-dir "$VS2_DB"

DRAM-v.py annotate \
    -i "${AMG_DIR}/vs2_for_dramv/for-dramv/final-viral-combined-for-dramv.fa" \
    -v "${AMG_DIR}/vs2_for_dramv/for-dramv/viral-affi-contigs-for-dramv.tab" \
    -o "${AMG_DIR}/DRAMv_annotation" \
    --threads 48 \
    --min_contig_size 1000

DRAM-v.py distill \
    -i "${AMG_DIR}/DRAMv_annotation/annotations.tsv" \
    -o "${AMG_DIR}/DRAMv_distill"
conda deactivate

echo "[INFO] $(date) AMG annotation complete"
echo "      VIBRANT AMGs: ${AMG_DIR}/VIBRANT/VIBRANT_*/VIBRANT_results_*/VIBRANT_AMG_individuals_*.tsv"
echo "      DRAM-v AMGs:  ${AMG_DIR}/DRAMv_distill/amg_summary.tsv"
