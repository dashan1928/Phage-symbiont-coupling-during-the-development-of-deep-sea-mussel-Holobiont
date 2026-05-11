#!/bin/bash
# ============================================================
# 04_viral_discovery.slurm
# Per-sample viral mining: geNomad + VirSorter2 + CheckV
# Input :  $ASSEMBLY_DIR/Assembly_<SAMPLE>/final.contigs.fa
# Output:  $VIRUS_DIR/verified_vOTUs_<SAMPLE>.fa
# ============================================================
#SBATCH --job-name=ViralDisc
#SBATCH --output=logs/virome_%A_%a.out
#SBATCH --error=logs/virome_%A_%a.err
#SBATCH --partition=normal,high_memory
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --array=1-48

set -euo pipefail
source ./config.sh
source "${CONDA_BASE}/etc/profile.d/conda.sh"

export CHECKVDB="$CHECKV_DB"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLES_LIST")
CONTIGS="${ASSEMBLY_DIR}/Assembly_${SAMPLE}/final.contigs.fa"
WORK="${VIRUS_DIR}/${SAMPLE}"
mkdir -p "$WORK"
cd "$WORK"

# 1. geNomad end-to-end
echo "[INFO] $(date) $SAMPLE — geNomad"
conda activate genomad
genomad end-to-end \
    --cleanup --splits 8 \
    "$CONTIGS" "genomad_${SAMPLE}" "$GENOMAD_DB"
conda deactivate

# 2. VirSorter2
echo "[INFO] $(date) $SAMPLE — VirSorter2"
conda activate vs2
virsorter run \
    -i "$CONTIGS" \
    -w "vs2_${SAMPLE}" \
    -j 48 \
    --include-groups dsDNAphage,ssDNA \
    --min-score 0.7 \
    --min-length 1500 \
    --provirus-off \
    --db-dir "$VS2_DB"
conda deactivate

# 3. Merge predictions
cat "genomad_${SAMPLE}/${SAMPLE}_summary/${SAMPLE}_virus.fna" \
    "vs2_${SAMPLE}/final-viral-combined.fa" \
    > "all_putative_viruses_${SAMPLE}.fa"

# 4. CheckV quality assessment
echo "[INFO] $(date) $SAMPLE — CheckV"
conda activate checkv
checkv end_to_end \
    "all_putative_viruses_${SAMPLE}.fa" \
    "checkv_${SAMPLE}" \
    -t 48 -d "$CHECKV_DB"
conda deactivate

# 5. Retain CheckV-verified viruses + proviruses
cat "checkv_${SAMPLE}/viruses.fna" "checkv_${SAMPLE}/proviruses.fna" \
    > "${VIRUS_DIR}/verified_vOTUs_${SAMPLE}.fa"

echo "[INFO] $(date) $SAMPLE viral discovery complete"
