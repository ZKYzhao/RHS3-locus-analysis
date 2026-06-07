#!/bin/bash
# Align CDS (high-accuracy MAFFT) and build the tree with IQ-TREE (ModelFinder + ultrafast bootstrap).
set -euo pipefail

CDS=${1:-CDSs.fa}
ALN=${CDS%.fa}.fas
THREADS=${2:-AUTO}

mafft --maxiterate 500 --genafpair --thread 24 "$CDS" > "$ALN"
iqtree3 -s "$ALN" -m MFP -bb 2000 -nm 5000 -redo -nt "$THREADS" -keep-ident
