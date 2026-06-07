# Comparative genomic and phylogenetic analysis of the RHS3 locus

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20579120.svg)](https://doi.org/10.5281/zenodo.20579120)

Code accompanying **X. He *et al.*, "A tripartite genetic conflict system controls hybrid sterility in
rice," *Science* (2026)**.

Downstream bioinformatics code for analyzing the **RHS3 locus** across multiple *Oryza* genomes.
Every script is parameterized (no hard-coded private paths or sample names) so the workflow can be
reproduced on public data; each file carries a one-line purpose header.

## Workflow

Inputs: *Oryza* genome assemblies (FASTA) and reference ORF/marker sequences of the RHS3 locus.

```
                       reference ORF/marker sequences (FASTA)
                                │
        ┌───────────────────────┴───────────────────────┐
        ▼                                                ▼
 [01] marker BLAST                                [03] CDS phylogeny
 locate each RHS3 ORF in every genome             exonerate CDS → MAFFT → IQ-TREE
        │
        ├──────────────► [02] synteny
        │                hits → GFF → BED → orthology blocks → JCVI plot
        │
        └──────────────► [04] region phylogeny
                         extract locus span by markers → progressiveMauve → IQ-TREE
```

Modules are independent; modules 02 and 04 reuse the marker hits / locus span from module 01.

## Modules

| Directory | Purpose | Main software |
|-----------|---------|---------------|
| `01_marker_blast/` | Locate RHS3 ORFs/markers in each genome; output a hit table | BLAST+ |
| `02_synteny/` | Build cross-species orthology blocks and draw the synteny plot | JCVI |
| `03_cds_phylogeny/` | Predict locus CDS per genome, align and build a tree | exonerate, MAFFT, IQ-TREE |
| `04_region_phylogeny/` | Extract the locus span, collinear-align and build a tree | progressiveMauve, IQ-TREE |

Numbered scripts run in order; `run_*.sh` chains a module end to end.

**Transposable element / repeat / miRNA annotation.** Transposable elements were identified by BLAST
against the Rice Genome Annotation Project database (https://rice.uga.edu/analyses_search_blast.shtml),
with the locus span as query; miRNAs and repeats were annotated with Geneious. These steps use external
/ interactive tools and are reported here for completeness.

## Environment

Tools and Python libraries are listed in [`environment.yml`](environment.yml) (compatible reference
versions are pinned; set `mauve` and `iqtree` to the versions you actually used).

```bash
conda env create -f environment.yml
conda activate rhs3
```

Species names follow valid *Oryza* nomenclature (e.g. *Oryza glumaepatula*); keep code, figures and
text consistent.

## Data availability

This repository contains code only. Genome assemblies are from public databases (accessions in the
paper and supplement); bulky intermediates are not version-controlled.

## Citation

> X. He *et al.*, "A tripartite genetic conflict system controls hybrid sterility in rice," *Science* (2026).

Archived code: Zenodo (2026); https://doi.org/10.5281/zenodo.20579120.

## License

MIT, see [`LICENSE`](LICENSE).
