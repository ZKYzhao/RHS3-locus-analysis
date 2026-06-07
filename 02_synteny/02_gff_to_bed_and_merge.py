#!/usr/bin/env python3
"""Convert per-species marker GFF to BED and merge, for building synteny blocks.

Usage:
    python 02_gff_to_bed_and_merge.py --gff-dir gffs/ --bed-dir beds/ --merged merged_all.bed
"""

import argparse
import glob
import os

import pandas as pd

GFF_COLS = ["chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]


def gff_to_bed(gff_file, bed_file):
    df = pd.read_csv(gff_file, sep="\t", header=None, names=GFF_COLS)
    df["name"] = df["attributes"].str.extract(r"ID=([^;]+)")
    pd.DataFrame({
        "chr": df["chr"],
        "start": df["start"] - 1,   # BED is 0-based, half-open
        "end": df["end"],
        "name": df["name"],
        "score": df["score"],
        "strand": df["strand"],
    }).to_csv(bed_file, sep="\t", header=False, index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff-dir", default="gffs")
    ap.add_argument("--bed-dir", default="beds")
    ap.add_argument("--merged", default="merged_all.bed")
    args = ap.parse_args()

    os.makedirs(args.bed_dir, exist_ok=True)
    frames = []
    for gff_file in sorted(glob.glob(os.path.join(args.gff_dir, "*.gff"))):
        species = os.path.splitext(os.path.basename(gff_file))[0]
        bed_file = os.path.join(args.bed_dir, f"{species}.bed")
        gff_to_bed(gff_file, bed_file)
        frames.append(pd.read_csv(bed_file, sep="\t", header=None))

    pd.concat(frames, ignore_index=True).to_csv(args.merged, sep="\t", header=False, index=False)
    print(f"[DONE] merged BED -> {args.merged}")


if __name__ == "__main__":
    main()
