#!/bin/bash
# Region phylogeny: collinear alignment (progressiveMauve) -> rename -> IQ-TREE.
set -euo pipefail

REGION=${1:-region.fa}          # merged region FASTA from step 01
ALIGNED=aligned.xmfa
RENAMED=aligned_renamed.fasta

progressiveMauve --seed-weight=15 --weight=10000 --collinear \
    --output="$ALIGNED" "$REGION"

python 02_rename_mauve_alignment.py --input "$REGION" --aligned "$ALIGNED" --out "$RENAMED"

iqtree3 -s "$RENAMED" -m MFP -bb 2000 -nm 2000 -redo -nt AUTO -keep-ident
