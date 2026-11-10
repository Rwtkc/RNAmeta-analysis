# RNAmeta analysis source and manuscript reproducibility

This repository contains the shared RNAmeta analysis implementations, the RNAenrich comparison core and scripts and compact annotation outputs used for the manuscript analyses.

RNAmeta: https://rnainformatics.cn/RNAmeta/

## Contents

| Directory | Contents |
| --- | --- |
| `online/` | R source for site annotation, transcript and exon features, metagene and boundary profiles, motif discovery, transcript mapping and structure analysis, extracted from the deployed modular application |
| `api/meta_plot_api.R` | RNAenrich metagene service using the shared Online implementations |
| `frontend/` | RNAenrich exact-coordinate, interval and background statistics, and React source components for RNAenrich, Meta Plot and RNAplot |
| `analyses/` | Original RBP analysis scripts, saved annotation/profile inputs, statistical outputs, feature/performance measurements and cross-environment consistency results |
| `resources/` | Reference annotation and public dataset manifests |
| `scripts/` | R syntax checking, RNAenrich core checks and a standalone statistical reproduction command |

The source snapshot was extracted on 10 November 2026. `SOURCE_MANIFEST.json` records the SHA-256 hash of every published source/input file. No reference genomes, user uploads, credentials, server configuration or third-party executables are included. These are analysis sources and selected interface components, not a complete website deployment bundle. RNAplot rendering-service source is not included in this snapshot. The runnable local analysis interfaces are distributed separately:

- R package: https://github.com/Rwtkc/RNAmeta_R
- Desktop: https://github.com/Rwtkc/RNAmeta_desktop

## Reproduce the key RBP statistics

Use an environment containing the Python dependencies in `requirements.txt`. From the repository root:

```text
python scripts/reproduce_statistics.py
```

The command recomputes regional Jensen-Shannon distances, equal-weight four-view Hellinger distances, 10000 protein-label permutations, all 114 functional association tests and FAM120A gene-paired bootstrap intervals. It verifies these against the archived results and writes fresh outputs to `outputs/statistics/`. It uses the included per-dataset annotations, transcript-profile matrices and per-gene region counts, without repeating peak calling or reference annotation.

Original scripts under `analyses/` preserve their original relative input paths. Some original scripts additionally analyse RNA-seq responses or regenerate complete figures and therefore require inputs beyond the compact data included here. The standalone command above reproduces the key manuscript statistics from the bundled inputs.

## Analysis dependencies and reference inputs

The Online source uses R 4.3.0 and Shiny 1.7.4.1, GenomicRanges 1.52.0, GenomicFeatures 1.52.1, rtracklayer 1.60.0, genomation 1.32.0, Biostrings 2.68.1 and Rsamtools 2.16.0. Sequence tools are MEME Suite 5.5.5, HOMER 4.11 and RNAfold 2.4.7. Other imported packages are declared in the source files. Public precomputation used R 4.4.3.

Reference resources: https://rnainformatics.cn/RNAmeta/downloads

Public BED files and precomputed outputs: https://rnainformatics.org.cn/RNAmod2/download.php

Reference identifiers must agree between BED coordinates, annotation objects and FASTA files. The Online code expects reference subdirectories such as `data/txdb`, `data/gff`, `data/txlens` and `data/fa`. Reference files are provided separately because of their size. `RNAMETA_ONLINE_ROOT` selects the source directory for the metagene API, and `METAPLOT_API_PORT` selects its listening port. For example, in a shell with the R dependencies and references prepared:

```sh
RNAMETA_ONLINE_ROOT="$PWD/online" Rscript api/meta_plot_api.R
```

This starts the analysis service. Production routing, authentication, storage and worker configuration belong to the deployment environment and are not provisioned by this repository.

## Source checks

```text
Rscript scripts/check_r_sources.R .
node scripts/check_rnaenrich.cjs
```

The Node check requires a locally available `typescript` module. Alternatively, set `TYPESCRIPT_MODULE` to its absolute module path. It tests chromosome aliases, duplicate removal, invalid intervals, exact matching, touching versus overlapping intervals and background-aware statistics.

## Manuscript analysis definitions

Regional composition distances use Jensen-Shannon distance with natural logarithms. The four-view distance is the root mean square of Hellinger distances over regional composition, mRNA and 5′ and 3′ splice-boundary profiles. Each view has equal weight. The other-protein reference averages the row and column medians after excluding the matched protein.

Functional comparisons use protein-level averages across available cell lines. Nineteen functions with at least five annotated proteins and five other proteins are compared across six regions, with two-sided Mann-Whitney tests and one Benjamini-Hochberg correction across all 114 tests. Cliff’s delta is `2U/(n1*n2)-1`.

Peak-count sensitivity uses thresholds of 0, 100, 250, 500, 1000, 2000 and 5000 peaks. Equal-depth sampling retains 43 peaks per dataset, the minimum among the 146 datasets for 73 shared proteins, without replacement for 200 repetitions. These are sensitivity analyses, not an inclusion requirement for the primary 223-dataset analysis.

FAM120A expression matching uses baseline ENCODE RNA-seq files ENCFF286GLL and ENCFF986DBN for K562 and ENCFF533XPJ and ENCFF321JIT for HepG2. Genes require mean TPM ≥1 in both cell lines and a higher-to-lower expression ratio ≤2. Shared targets additionally have at least one same-strand peak in each cell line. Per-gene regional fractions give genes equal weight. The 95% intervals use the 2.5th and 97.5th percentiles from 10000 gene-level bootstrap resamples, with seed 20261109.

Benchmark measurements describe the respective deployed platforms and their native task scopes. They compare practical runtime against RNAmod, with cache states, thread settings and individual runs documented in the included tables.

## Licensing and external tools

RNAmeta code in this repository is provided under GPL-3.0-or-later. External R packages, MEME Suite, HOMER, ViennaRNA, IGV and JavaScript dependencies retain their own licenses and are not relicensed here. Public data remain subject to their original source terms and should be cited through the source manifests.
