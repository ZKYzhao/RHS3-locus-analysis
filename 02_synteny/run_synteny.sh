#!/bin/bash
# Synteny pipeline: marker hits -> GFF -> BED -> orthology blocks -> colors -> JCVI plot.
set -euo pipefail

MATCHES=${1:-best_matches_results.csv}   # marker hit table from module 01

python 01_matches_to_gff.py        --matches "$MATCHES" --out-dir gffs
python 02_gff_to_bed_and_merge.py  --gff-dir gffs --bed-dir beds --merged merged_all.bed
python 03_build_orthology_blocks.py --bed-dir beds --out blocks
python 04_make_blocks_layout.py    --bed-dir beds --out blocks.layout
python 05_make_color_config.py     --blocks blocks --out colors.txt

python -m jcvi.graphics.synteny blocks merged_all.bed blocks.layout \
    --format svg --glyphstyle=arrow --font Arial --dpi 500 --figsize 16x14 \
    --genelabelrotation 45 --glyphcolor custom --colorfile colors.txt
