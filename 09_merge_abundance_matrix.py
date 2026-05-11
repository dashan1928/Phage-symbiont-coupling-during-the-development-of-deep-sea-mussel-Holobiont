#!/usr/bin/env python3
"""
09_merge_abundance_matrix.py
Merge per-sample CoverM outputs into wide TPM and count matrices.

Output:
    results/virus_only_tpm_matrix.csv          (3324 vOTU x 48 sample)
    results/bacteria_host_MAG_TPM.tsv          (23 MAG x 48 sample)
    results/Final_Holobiont_TPM_Matrix.csv     (combined, 19457 contig x 48)
    results/holobiont_count_matrix.csv         (raw counts, for DESeq2)
"""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--coverm_dir", required=True,
                   help="Directory containing <SAMPLE>_abundance.tsv files")
    p.add_argument("--samples", required=True, help="samples.txt (one ID per line)")
    p.add_argument("--virus_fa", required=True, help="vOTU catalogue fasta")
    p.add_argument("--mag_fa",   required=True, help="MAG catalogue fasta")
    p.add_argument("--out_dir",  required=True)
    return p.parse_args()


def fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def main():
    args = parse_args()
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    samples = [s.strip() for s in open(args.samples) if s.strip()]

    # Identify which contigs are viruses vs MAGs
    virus_ids = set(fasta_ids(args.virus_fa))
    mag_ids   = set(fasta_ids(args.mag_fa))

    tpm_frames, count_frames = [], []
    for s in samples:
        f = Path(args.coverm_dir) / f"{s}_abundance.tsv"
        if not f.exists():
            sys.exit(f"[ERROR] missing CoverM output for {s}: {f}")

        df = pd.read_csv(f, sep="\t")
        # CoverM column names: 'Contig', '<sample> TPM', '<sample> Read Count', '<sample> Covered Fraction'
        df = df.rename(columns={df.columns[0]: "Contig"})
        tpm_col   = [c for c in df.columns if "TPM" in c][0]
        count_col = [c for c in df.columns if "Read Count" in c][0]

        tpm_frames.append(df.set_index("Contig")[[tpm_col]].rename(columns={tpm_col: s}))
        count_frames.append(df.set_index("Contig")[[count_col]].rename(columns={count_col: s}))

    tpm_mat   = pd.concat(tpm_frames,   axis=1).fillna(0.0)
    count_mat = pd.concat(count_frames, axis=1).fillna(0).astype(int)

    # Re-order contigs: viruses first, then MAGs (preserve fasta order)
    ordered = [c for c in fasta_ids(args.virus_fa) + fasta_ids(args.mag_fa)
               if c in tpm_mat.index]
    tpm_mat   = tpm_mat.loc[ordered, samples]
    count_mat = count_mat.loc[ordered, samples]

    # Split outputs
    virus_idx = [i for i in tpm_mat.index if i in virus_ids]
    mag_idx   = [i for i in tpm_mat.index if i in mag_ids]

    tpm_mat.loc[virus_idx].to_csv(out / "virus_only_tpm_matrix.csv",
                                  index_label="Contig_norm")
    tpm_mat.loc[mag_idx].to_csv(out / "bacteria_host_MAG_TPM.tsv",
                                index_label="MAG_ID", sep="\t")
    tpm_mat.to_csv(out / "Final_Holobiont_TPM_Matrix.csv",
                   index_label="Contig")
    count_mat.to_csv(out / "holobiont_count_matrix.csv",
                     index_label="Contig")

    print(f"[OK] Wrote 4 matrices to {out}")
    print(f"     vOTU TPM   : {tpm_mat.loc[virus_idx].shape}")
    print(f"     MAG TPM    : {tpm_mat.loc[mag_idx].shape}")
    print(f"     Combined   : {tpm_mat.shape}")
    print(f"     Counts     : {count_mat.shape}")


if __name__ == "__main__":
    main()
