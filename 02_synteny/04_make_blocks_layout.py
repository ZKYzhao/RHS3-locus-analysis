#!/usr/bin/env python3
"""Write blocks.layout (track positions, colors and edges) for jcvi.graphics.synteny.

Usage:
    python 04_make_blocks_layout.py --bed-dir beds/ --out blocks.layout
    python 04_make_blocks_layout.py --order speciesA,speciesB,... --out blocks.layout
"""

import argparse
import glob
import os

TRACK_COLORS = ["m", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"]


def write_layout(species_list, outfile):
    n = len(species_list)
    y_top, y_bottom = 0.9, 0.1
    step = (y_top - y_bottom) / (n - 1) if n > 1 else 0
    with open(outfile, "w") as f:
        f.write("# x,   y, rotation,   ha,     va,   color, ratio,     label\n")
        for i, species in enumerate(species_list):
            f.write("{0}, {1}, {2:8}, {3:4}, {4:6}, {5:>7}, {6:5}, {7}\n".format(
                0.6, y_top - i * step, 0, "leftalign", "center",
                TRACK_COLORS[i % len(TRACK_COLORS)], 1, species.capitalize()))
        f.write("# edges\n")
        for i in range(1, n):
            f.write(f"e, {i - 1}, {i}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bed-dir", default="beds")
    ap.add_argument("--order", help="comma-separated species order (top to bottom)")
    ap.add_argument("--out", default="blocks.layout")
    args = ap.parse_args()

    if args.order:
        species_list = args.order.split(",")
    else:
        species_list = [os.path.splitext(os.path.basename(p))[0]
                        for p in sorted(glob.glob(os.path.join(args.bed_dir, "*.bed")))]
    write_layout(species_list, args.out)
    print(f"[DONE] {len(species_list)} tracks -> {args.out}")


if __name__ == "__main__":
    main()
