#!/usr/bin/env python3
"""Write colors.txt for jcvi --colorfile: highlight key genes, others light gray.

Usage:
    python 05_make_color_config.py --blocks blocks --out colors.txt
    python 05_make_color_config.py --blocks blocks --highlight "ORF4(JIA)=#4D4DFF,ORF10(MAO)=#FF0000"

Matching is by substring, so a marker named ORF4(JIA) or ORF4(JIA)-gene both match.
"""

import argparse

import pandas as pd

DEFAULT_HIGHLIGHT = {"ORF4(JIA)": "#4D4DFF", "ORF10(MAO)": "#FF0000", "ORF16(DUN)": "#00FF00"}
DEFAULT_COLOR = "#C0C0C0"


def color_of(marker, highlight):
    for key, color in highlight.items():
        if key in marker:
            return color
    return DEFAULT_COLOR


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--blocks", default="blocks")
    ap.add_argument("--out", default="colors.txt")
    ap.add_argument("--highlight", help='"name=color,name=color"')
    args = ap.parse_args()

    highlight = dict(DEFAULT_HIGHLIGHT)
    if args.highlight:
        highlight = dict(pair.split("=") for pair in args.highlight.split(","))

    blocks = pd.read_csv(args.blocks, sep="\t", header=None)
    with open(args.out, "w") as f:
        for _, row in blocks.iterrows():
            for gene in row:
                if gene != "." and pd.notna(gene):
                    f.write(f"{gene}\t{color_of(gene.split('_')[-1], highlight)}\n")
    print(f"[DONE] -> {args.out}")


if __name__ == "__main__":
    main()
