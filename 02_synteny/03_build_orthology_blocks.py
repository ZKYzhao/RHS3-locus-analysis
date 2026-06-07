#!/usr/bin/env python3
"""Align per-species BED hits into an orthology block matrix (rows=marker, cols=species).

Usage:
    python 03_build_orthology_blocks.py --bed-dir beds/ --out blocks

Missing orthologs are filled with '.'; output feeds jcvi.graphics.synteny.
"""

import argparse
import glob
import os

import pandas as pd


def marker_of(name):
    return name.split("_")[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bed-dir", default="beds")
    ap.add_argument("--out", default="blocks")
    args = ap.parse_args()

    per_species = {}
    for bed_file in sorted(glob.glob(os.path.join(args.bed_dir, "*.bed"))):
        species = os.path.splitext(os.path.basename(bed_file))[0]
        names = pd.read_csv(bed_file, sep="\t", header=None).iloc[:, 3]
        per_species[species] = {marker_of(n): n for n in names}

    all_markers = set()
    for genes in per_species.values():
        all_markers.update(genes)

    matrix = {m: {sp: genes.get(m, ".") for sp, genes in per_species.items()} for m in all_markers}
    df = pd.DataFrame.from_dict(matrix, orient="index").sort_index()
    df.to_csv(args.out, sep="\t", index=False, header=False)
    print(f"[DONE] {df.shape[0]} markers x {df.shape[1]} species -> {args.out}")


if __name__ == "__main__":
    main()
