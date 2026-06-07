#!/usr/bin/env python3
"""Extract predicted CDS (exons only) from exonerate alignments; introns/UTR are lowercase and dropped.

Usage:
    python 02_extract_cds_from_exonerate.py --gff-dir exonerate_result/ --out CDSs.fa

Sequences are named <genome>_<target>.
"""

import argparse
import glob
import os
import re

DROP = str.maketrans("", "", ".-atcg")  # remove lowercase bases and alignment fillers


def extract_from_block(block):
    name = re.findall(r"Target:(.+)\n", block)[0].strip()
    seq = "".join(re.findall(r"\|.+\n\s*\d+ : (.+) :", block)).translate(DROP)
    return name, seq


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gff-dir", required=True)
    ap.add_argument("--out", default="CDSs.fa")
    args = ap.parse_args()

    records = []
    for gff in sorted(glob.glob(os.path.join(args.gff_dir, "*.gff"))):
        label = os.path.basename(gff).split(".")[0]
        for block in open(gff).read().split("C4 Alignment:")[1:]:
            name, seq = extract_from_block(block)
            if seq:
                records.append(f">{label}_{name}\n{seq}")

    with open(args.out, "w") as f:
        f.write("\n".join(records) + "\n")
    print(f"[DONE] {len(records)} CDS -> {args.out}")


if __name__ == "__main__":
    main()
