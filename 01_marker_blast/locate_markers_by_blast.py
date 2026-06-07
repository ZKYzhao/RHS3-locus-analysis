#!/usr/bin/env python3
"""Locate RHS3 ORF/marker sequences in each genome by BLAST; report position, strand and coverage-weighted match.

Usage:
    python locate_markers_by_blast.py --query markers.fa --genome-dir genomes/ --out-dir blast_results/

Per genome a *_blast.csv is written; MatchingDegree = (length - mismatches - gaps) / qlen * 100.
"""

import argparse
import os
import subprocess

import pandas as pd

OUTFMT_FIELDS = [
    "qseqid", "sseqid", "sseq", "sstrand", "length", "qlen",
    "mismatch", "gaps", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]

# Relaxed gap penalties to recover divergent / indel-rich homologous segments.
BLAST_PARAMS = [
    "-sorthsps", "4", "-dust", "no",
    "-penalty", "-2", "-gapopen", "1", "-gapextend", "2", "-reward", "1",
    "-xdrop_gap", "500", "-xdrop_gap_final", "3000", "-culling_limit", "3",
]

MIN_MATCH_DEGREE = 20  # drop hits below this coverage-weighted match


def make_blast_db(genome_fasta):
    subprocess.run(["makeblastdb", "-in", genome_fasta, "-dbtype", "nucl", "-parse_seqids"],
                   check=True, capture_output=True, text=True)


def run_blastn(query, genome_fasta, out_file, threads):
    cmd = ["blastn", "-query", query, "-db", genome_fasta, "-out", out_file,
           "-outfmt", "10 " + " ".join(OUTFMT_FIELDS), "-num_threads", str(threads)] + BLAST_PARAMS
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def parse_hits(raw_csv, genome_label):
    if os.path.getsize(raw_csv) == 0:
        return pd.DataFrame()
    df = pd.read_csv(raw_csv, header=None)
    df.columns = ["marker", "chrome", "sseq", "sstrand", "length", "qlen", "mismatches",
                  "gaps", "q_start", "q_end", "s_start", "s_end", "evalue", "bitscore"]
    df["line"] = genome_label
    df[["s_start", "s_end"]] = df.apply(
        lambda r: pd.Series([min(r["s_start"], r["s_end"]), max(r["s_start"], r["s_end"])]), axis=1)
    df["MatchingDegree"] = (df["length"] - df["mismatches"] - df["gaps"]) / df["qlen"] * 100
    df = df.sort_values("MatchingDegree", ascending=False)
    df.to_csv(raw_csv, index=False)
    return df[df["MatchingDegree"] >= MIN_MATCH_DEGREE]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--query", required=True, help="marker sequences (FASTA)")
    ap.add_argument("--genome-dir", required=True, help="genome assemblies")
    ap.add_argument("--out-dir", default="blast_results")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    all_hits = []
    for fname in os.listdir(args.genome_dir):
        if not fname.endswith((".fasta", ".fa", ".fna")):
            continue
        genome = os.path.join(args.genome_dir, fname)
        label = os.path.splitext(fname)[0]
        out_csv = os.path.join(args.out_dir, f"{label}_blast.csv")
        make_blast_db(genome)
        run_blastn(args.query, genome, out_csv, args.threads)
        hits = parse_hits(out_csv, label)
        if not hits.empty:
            all_hits.append(hits)
        print(f"[OK] {label}: {len(hits)} hits >= {MIN_MATCH_DEGREE}")

    if all_hits:
        merged = pd.concat(all_hits, ignore_index=True)
        merged.to_csv(os.path.join(args.out_dir, "best_matches_results.csv"), index=False)
        print(f"[DONE] best_matches_results.csv ({len(merged)} rows)")


if __name__ == "__main__":
    main()
