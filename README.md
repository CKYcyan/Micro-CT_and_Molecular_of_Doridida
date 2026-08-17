# Micro-CT and Molecular Data for Doridida

This repository contains lightweight analysis files associated with the manuscript:

**Micro-CT and molecular phylogenetics suggest evolutionary patterns in spicule and shell architecture within Doridida (Gastropoda: Heterobranchia)**

Raw Micro-CT volumes and specimen images are archived separately in MorphoSource.

## Repository Layout

- `Doridida spicule and phylogeny/` - public lightweight analysis package, including morphometric input data, scripts, metadata, and reproducible result tables.
- `Molecular_Phylogeny/` - source files and outputs for ML and BI phylogenetic analyses.
- `Trace_History/` - ancestral-state reconstruction and trace-history source outputs.
- `Morphological_Analysis/` - historical morphometric-analysis files retained for audit and comparison.

## Current Final-Stage Note

As of 2026-08-17, the active morphometric input is:

- `Doridida spicule and phylogeny/data/analysis.csv`

The current Python audit outputs are:

- `Doridida spicule and phylogeny/results/morphometric_analysis_20260806_R002_R003_D/`

The public R workflows are:

- `Doridida spicule and phylogeny/R/analysis_microct_spicule_arrangement.R` for PCA, PERMANOVA, PERMDISP, pairwise PERMANOVA, Kruskal-Wallis and Dunn-style post-hoc tables.
- `Doridida spicule and phylogeny/R/plot_fig3_ggplot2_final_20260807.R` for the final Fig. 3 ggplot2 layout based on the corrected R002/R003 coding.

Alignment files and BI/ML tree files are tracked in the top-level `Molecular_Phylogeny/` folder. Ancestral-state reconstruction files are tracked in `Trace_History/`.
