#!/bin/bash
# ============================================================
# 06_taxonomy_and_host_prediction.slurm
# (a) GTDB-Tk classification of dereplicated bacterial MAGs
# (b) iPHoP host prediction for the dereplicated vOTU catalogue
# ============================================================
#SBATCH --job-name=Tax_Host
#SBATCH --output=logs/tax_host_%j.out
#SBATCH --error=logs/tax_host_%j.err
#SBATCH --partition=normal,high_memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=256G
#SBATCH --time=48:00:00

set -euo pipefail
source ./config.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"

# ── (a) GTDB-Tk taxonomic classification ─────────────────────
echo "[INFO] $(date) GTDB-Tk classify_wf on dereplicated MAGs"
conda activate gtdbtk
gtdbtk classify_wf \
    --genome_dir "${GLOBAL_DIR}/dRep_out/dereplicated_genomes/" \
    --extension fa \
    --out_dir "${RESULTS_DIR}/gtdbtk_out" \
    --cpus 48 \
    --skip_ani_screen
conda deactivate
echo "[INFO] $(date) GTDB-Tk output: ${RESULTS_DIR}/gtdbtk_out/gtdbtk.bac120.summary.tsv"

# ── (b) iPHoP host prediction ────────────────────────────────
VIRUS_FA="${GLOBAL_DIR}/Global_Dereplicated_vOTUs.fasta"
IPHOP_OUT="${RESULTS_DIR}/iPHoP_Output_Network"

echo "[INFO] $(date) iPHoP host prediction (composite score >= 90)"
conda activate iphop
iphop predict \
    -f "$VIRUS_FA" \
    --db_dir "$IPHOP_DB" \
    --out_dir "$IPHOP_OUT" \
    --num_threads 48
conda deactivate

# iPHoP outputs a sequence-cleaned vOTU file used by VIBRANT downstream
CLEAN_VOTUS="${IPHOP_OUT}/Global_Dereplicated_vOTUs_clean.fna"
[ -f "$CLEAN_VOTUS" ] || echo "[WARN] iPHoP clean fasta not found"

echo "[INFO] $(date) Host prediction table: ${IPHOP_OUT}/Host_prediction_to_genus_m90.csv"
echo "[INFO] $(date) Done"
