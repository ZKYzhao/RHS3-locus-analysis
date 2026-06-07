#!/usr/bin/env python3
"""Split the marker hit table (best_matches_results.csv) into per-species GFF for synteny.

Usage:
    python 01_matches_to_gff.py --matches best_matches_results.csv --out-dir gffs/

Each hit becomes one mRNA feature; ID=<species>_<marker> encodes orthology, match_score
stores the coverage-weighted match (not strict sequence identity).
"""

import argparse
import os

import pandas as pd

STRAND_MAP = {"plus": "+", "minus": "-"}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matches", required=True)
    ap.add_argument("--out-dir", default="gffs")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    matches = pd.read_csv(args.matches).sort_values(["line", "s_start"], ascending=False)

    for species in matches["line"].unique():
        sample = matches[matches["line"] == species]
        gff = pd.DataFrame({
            "seqid": sample["line"],
            "source": "marker_blast",
            "type": "mRNA",
            "start": sample["s_start"],
            "end": sample["s_end"],
            "score": sample["MatchingDegree"],
            "strand": sample["sstrand"].map(STRAND_MAP),
            "phase": ".",
            "attributes": (
                "gene_orientation . ; "
                "match_score " + sample["MatchingDegree"].round(2).astype(str) + " ; "
                "ID=" + species + "_" + sample["marker"]
            ),
        })
        gff.to_csv(os.path.join(args.out_dir, f"{species}.gff"), sep="\t", header=False, index=False)
        print(f"[OK] {species}: {len(gff)} features")


if __name__ == "__main__":
    main()
