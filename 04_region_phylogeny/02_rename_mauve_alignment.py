#!/usr/bin/env python3
"""Tidy progressiveMauve output: keep the main LCB, drop '=' separators, restore sample names from input headers.

Usage:
    python 02_rename_mauve_alignment.py --input region.fa --aligned aligned.xmfa --out aligned_renamed.fasta

Renaming relies on input order == alignment order; a count mismatch aborts to avoid mislabeled tree tips.
"""

import argparse

from Bio import SeqIO, SeqRecord, Seq


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="pre-alignment region FASTA (naming order)")
    ap.add_argument("--aligned", required=True, help="progressiveMauve output (xmfa)")
    ap.add_argument("--out", default="aligned_renamed.fasta")
    args = ap.parse_args()

    names = [rec.id for rec in SeqIO.parse(args.input, "fasta")]
    records = list(SeqIO.parse(args.aligned, "fasta-pearson"))
    if len(records) < len(names):
        raise SystemExit(f"[ERROR] alignment has {len(records)} < input {len(names)} sequences; "
                         f"check progressiveMauve output before renaming.")

    fixed = []
    for i, name in enumerate(names):
        clean = str(records[i].seq).replace("=", "")
        fixed.append(SeqRecord.SeqRecord(Seq.Seq(clean), id=name, description=""))
        print(f"  [{i}] {records[i].id} -> {name}  len={len(clean)}")  # mapping for manual check

    SeqIO.write(fixed, args.out, "fasta")
    print(f"[DONE] {len(fixed)} sequences -> {args.out}")


if __name__ == "__main__":
    main()
