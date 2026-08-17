# Doridida spicule and phylogeny

This repository contains lightweight data files and analysis code associated with the manuscript:

**Micro-CT and molecular phylogenetics suggest evolutionary patterns in spicule and shell architecture within Doridida (Gastropoda: Heterobranchia)**

Raw Micro-CT volumes and specimen images are stored outside this repository because they are large research data files. During review, they are available through temporary MorphoSource reviewer links:

- Micro-CT data: <https://www.morphosource.org/projects/000859624/temporary_link/4AAgNEtfpN6SUk9WqeJ39kYn?locale=en>
- Specimen images: <https://www.morphosource.org/projects/000874416/temporary_link/MQQEzPNZbJ8a4zBtjZPdHoZK?locale=en>

Permanent MorphoSource identifiers should replace these temporary links after publication or repository finalisation.

## Repository structure

- `R/` - R scripts for morphometric statistics and final Fig. 3 plotting.
- `data/` - lightweight input data used by the R analyses.
- `metadata/` - Micro-CT scan metadata prepared for MorphoSource.
- `scripts/` - auxiliary audit scripts used when the R environment is unavailable.
- `results/` - small reproducible result tables and summaries retained for audit.
- `trees/` - local notes for tree-file placement; original tree and alignment files are tracked in the top-level `../Molecular_Phylogeny/` folder.

## Morphometric analysis

The main morphometric script is:

```bash
Rscript R/analysis_microct_spicule_arrangement.R
```

The script reads `data/analysis.csv`, standardises the Micro-CT morphometric variables, and analyses the same Euclidean distance matrix using PERMANOVA and PERMDISP. It also writes PCA summaries, pairwise PERMANOVA comparisons, Kruskal-Wallis summaries, and optional Dunn post-hoc tests when the `dunn.test` package is available.

Required R package:

- `vegan`

Optional R package:

- `dunn.test`

The script writes generated output files to `results/`.


### Final Fig. 3 plotting workflow

The final ggplot2 script used to recreate the corrected Fig. 3 layout is:

```bash
Rscript R/plot_fig3_ggplot2_final_20260807.R
```

The script reads `data/analysis.csv` and writes PNG/PDF outputs to:

```text
results/fig3_20260807/
```

Required R packages:

- `ggplot2`
- `grid` (base R)

### Python audit workflow

When R is unavailable, the morphometric analyses can be audited with:

```bash
python scripts/morphometric_analysis_python_audit.py
```

The Python audit reads `data/analysis.csv` and writes PCA summaries, PERMANOVA, pairwise PERMANOVA, PERMDISP, Kruskal-Wallis, and Dunn-Holm tables to:

```text
results/morphometric_analysis_20260806_R002_R003_D/
```

These outputs are intended for final-stage consistency checks. The R/vegan workflow should remain the primary analysis source for final manuscript p-values where available.

## Data notes

`data/analysis.csv` contains four spicule-bearing arrangement categories coded in the manuscript data table. The current coding treats R002 and R003 as `D` (`Trabecular`), consistent with their Cadlinidae assignment. The analysis does not include the `None` category because shell-state analyses and non-spicule-bearing states were not part of the Micro-CT morphometric PERMANOVA/PERMDISP tests.

`metadata/` contains MorphoSource-oriented Micro-CT scan metadata extracted from scanner logs and associated notes.

## Tree files

The manuscript presents the Bayesian tree as the main phylogenetic result and provides the ML topology for comparison. To avoid duplicating large or parallel tree folders, the original alignment and tree files are tracked in the repository-level `../Molecular_Phylogeny/` folder, including the IQ-TREE output folder and the MrBayes consensus/MCC tree files. The local `trees/` folder is retained as a pointer for readers browsing only this lightweight package.
