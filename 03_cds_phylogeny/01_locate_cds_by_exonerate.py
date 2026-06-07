#!/usr/bin/env python3
"""Predict the reference CDS in each genome's target chromosome with exonerate (est2genome), in parallel.

Usage:
    python 01_locate_cds_by_exonerate.py --ref-cds ref_CDS.fa --genome-dir genomes/ \
        --out-dir exonerate_result/ --chr-pattern "Chromosome 3" --threads 10

For each genome: index it, extract the chromosome whose header matches --chr-pattern,
then run exonerate with the reference CDS as query, writing *.exResult.gff.
"""

import argparse
import os
import subprocess
import threading
from queue import Queue

EXONERATE_OPTS = ["--model", "est2genome", "--showtargetgff",
                  "--maxintron", "10000", "--percent", "80"]


def process_genome(genome_path, ref_cds, out_dir, chr_pattern):
    label = os.path.basename(genome_path).split(".")[0]

    index_file = f"{genome_path}.index"
    if not os.path.exists(index_file):
        subprocess.run(["fastaindex", genome_path, index_file], check=True)

    headers = subprocess.run(["grep", "^>", genome_path], capture_output=True, text=True).stdout
    target_header = next((h for h in headers.splitlines() if chr_pattern in h), "")
    if not target_header:
        print(f"[SKIP] {label}: no chromosome matching '{chr_pattern}'")
        return
    chr_name = target_header.replace(">", "").split()[0]

    chr_file = os.path.join(out_dir, f"temp_{label}.fa")
    with open(chr_file, "w") as out:
        subprocess.run(["fastafetch", genome_path, index_file, chr_name], stdout=out, check=True)

    out_gff = os.path.join(out_dir, f"{label}.exResult.gff")
    with open(out_gff, "w") as out:
        subprocess.run(["exonerate"] + EXONERATE_OPTS +
                       ["--query", ref_cds, "--target", chr_file], stdout=out, check=True)
    os.remove(chr_file)
    print(f"[OK] {label}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref-cds", required=True, help="reference CDS (FASTA)")
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--out-dir", default="exonerate_result")
    ap.add_argument("--chr-pattern", default="Chromosome 3", help="target chromosome header substring")
    ap.add_argument("--threads", type=int, default=10)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    queue = Queue()

    def worker():
        while True:
            item = queue.get()
            if item is None:
                break
            process_genome(item, args.ref_cds, args.out_dir, args.chr_pattern)
            queue.task_done()

    workers = [threading.Thread(target=worker) for _ in range(args.threads)]
    for t in workers:
        t.start()
    for fname in os.listdir(args.genome_dir):
        if fname.endswith((".fasta", ".fa", ".fna")):
            queue.put(os.path.join(args.genome_dir, fname))
    queue.join()
    for _ in workers:
        queue.put(None)
    for t in workers:
        t.join()
    print("[DONE]")


if __name__ == "__main__":
    main()
