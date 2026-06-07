#!/usr/bin/env python3
"""Cut the RHS3 locus from each genome: anchor marker fixes the chromosome, then take the span between two markers.

Usage:
    python 01_extract_region_between_markers.py --genome-dir genomes/ --blast-dir blast_results/ \
        --anchor "ORF4(JIA)-gene" --start "ORF3-gene" --end "ORF9-gene" --out region.fa

Marker coordinates come from the module-01 hit tables.
"""

import argparse
import glob
import os

import pandas as pd
from Bio import SeqIO


def region_bounds(blast_df, anchor, start_marker, end_marker):
    chrom = blast_df[blast_df["marker"] == anchor].iloc[0]["chrome"]
    on_chr = blast_df[blast_df["chrome"] == chrom]
    s = on_chr[on_chr["marker"] == start_marker].iloc[0]["s_start"]
    e = on_chr[on_chr["marker"] == end_marker].iloc[0]["s_end"]
    return chrom, min(s, e), max(s, e)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--blast-dir", required=True)
    ap.add_argument("--anchor", required=True, help="marker used to fix the chromosome")
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--out", default="region.fa")
    args = ap.parse_args()

    with open(args.out, "w") as out:
        for genome in sorted(glob.glob(os.path.join(args.genome_dir, "*.fasta"))):
            label = os.path.splitext(os.path.basename(genome))[0]
            blast_csv = os.path.join(args.blast_dir, f"{label}_blast.csv")
            if not os.path.exists(blast_csv):
                print(f"[SKIP] {label}: no hit table")
                continue

            df = pd.read_csv(blast_csv).sort_values("MatchingDegree", ascending=False)
            try:
                chrom, start, end = region_bounds(df, args.anchor, args.start, args.end)
            except IndexError:
                print(f"[WARN] {label}: anchor/start/end marker missing, skipped")
                continue

            seq = None
            for rec in SeqIO.parse(genome, "fasta"):
                if rec.id in (str(chrom), f"Chr{chrom}"):
                    seq = str(rec.seq[start - 1:end])
                    break
            if seq is None:
                print(f"[WARN] {label}: chromosome {chrom} not found")
                continue

            out.write(f">{label}\n{seq}\n")
            print(f"[OK] {label}: {chrom}:{start}-{end} ({len(seq)} bp)")


if __name__ == "__main__":
    main()
