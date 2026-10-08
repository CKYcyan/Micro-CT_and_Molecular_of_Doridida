# Doridida spicule and phylogeny

Analysis package for **Micro-CT and molecular phylogenetics suggest evolutionary patterns in spicule and shell architecture within Doridida (Gastropoda: Heterobranchia)**.

## R-only workflow

From this package directory, with Rscript available:

~~~sh
Rscript --vanilla R/install_dependencies.R
Rscript --vanilla tests/test_workflow.R
Rscript --vanilla R/analysis_microct_spicule_arrangement.R
Rscript --vanilla R/analysis_paired_regions.R
Rscript --vanilla R/plot_fig3_ggplot2_final_20260807.R
Rscript --vanilla R/validate_phylogeny_assets.R
~~~

The installer installs missing packages; it does not downgrade packages or freeze a library. Exact tested versions are recorded in metadata/R_environment_20261008.csv and each run's sessionInfo.txt. This release used R 4.5.3, vegan 2.7-6, permute 0.9-10, dunn.test 1.4.2, ggplot2 4.0.3, multcompView 0.1-12 and ape 5.8-1. Figure export requires R's Cairo graphics capability.

Scripts locate the package from their own path. Analysis/figure options are --input PATH, --output PATH, --seed INTEGER and --permutations INTEGER. Defaults: seed 20260707; 9999 permutations. Plotting is deterministic and records zero permutations.

Windows paths containing Chinese characters require a supported UTF-8 locale. If Windows R is configured with the unsupported C.UTF-8 locale, set LC_ALL, LC_CTYPE and LANG to English_United States.utf8 before starting Rscript. BOM and BOM-free CSV reading were also tested under the C locale using ASCII paths.

## Data and methods

The canonical input, data/analysis.csv, contains 51 unique specimens and 12 finite morphometric variables: SpA1=21, SpA2=11, SpA3=15, SpA4=4. R002/R003 are SpA3 (Trabecular). The loader accepts B-E aliases, mapping them to SpA1-SpA4. Invalid IDs, unknown labels, invalid numeric values and constant responses produce errors.

- analysis_microct_spicule_arrangement.R: standardized 12-variable PCA; Euclidean PERMANOVA; median PERMDISP; six pairwise PERMANOVA contrasts with Holm adjustment; 12 Kruskal-Wallis tests with Holm adjustment; two-sided Dunn tests with Holm adjustment within each variable's six contrasts.
- analysis_paired_regions.R: new R implementation of the recorded paired design. Six metrics are standardized jointly across 102 region measurements. Specimen averages test SpA; mantle-minus-foot differences test SpA x Region; whole-specimen sign flips test Region. The independent unit remains the specimen. These separate-space R2 values are not additive.
- plot_fig3_ggplot2_final_20260807.R: retains the entry-point name, uses SpA labels and recalculates significance letters. Negative pattern-factor values are shown. Cairo PDF export embeds fonts for micrometre symbols.
- validate_phylogeny_assets.R: alignment/tree structure audit. Optional positional arguments: repository root and output directory. REVIEW rows are unresolved items; zero exit means the audit completed, not that those items passed.

The paired script is newly implemented, not a recovered original. See [provenance](metadata/ANALYSIS_PROVENANCE.md).

## Outputs and historical records

Current outputs: results/r_onefactor_20261008/, results/r_paired_20261008/, results/fig3_R_20261008/ and results/r_phylogeny_validation_20261008/. Runs record input MD5, seed, permutations, package versions and R session information.

See [validation report](metadata/VALIDATION_20261008.md) and metadata/manuscript_reconciliation.csv before updating the manuscript.

Older result folders and the auxiliary Python file in scripts/ are historical records, excluded from this R-only workflow and validation evidence. archive/R/ preserves three earlier local R figure scripts byte-for-byte; their older logic and machine-specific paths are not supported entry points.

Tree/alignment inputs are in ../Molecular_Phylogeny/; trace files are in ../Trace_History/. A complete rerun of IQ-TREE, MrBayes, MCC generation or ancestral-state reconstruction is outside this R morphometric validation.

## Raw data access

Existing temporary review links are retained:

- [Micro-CT volumes](https://www.morphosource.org/projects/000859624/temporary_link/4AAgNEtfpN6SUk9WqeJ39kYn?locale=en)
- [Specimen images](https://www.morphosource.org/projects/000874416/temporary_link/MQQEzPNZbJ8a4zBtjZPdHoZK?locale=en)

These links were not tested in this code audit. Replace temporary links with permanent MorphoSource identifiers when available. Scan metadata remain in metadata/.
