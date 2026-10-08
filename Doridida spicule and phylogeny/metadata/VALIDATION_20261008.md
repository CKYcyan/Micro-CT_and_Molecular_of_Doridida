# Validation report — 2026-10-08

## Scope and outcome

The canonical R morphometric workflow was corrected and executed using R 4.5.3. The release contains both the executable R code and dated R-generated results. No Python analysis or validation was used for this release. Earlier audit/Python artifacts are historical, not evidence for these results.

## Checks completed

- 51 unique specimens; group sizes 21/11/15/4; 12 finite response columns. R002/R003 remain SpA3. Numerical measurements and specimen IDs match the pre-release R input; only B-E labels became SpA1-SpA4.
- 77 R regression checks: BOM/no-BOM and B-E/SpA equivalence; invalid input/CLI rejection; independent rank/tie-corrected Dunn Z, two-sided p-values, Holm adjustment and all 72 compact-letter decisions; pooled scaling and pair reconstruction; independent sums-of-squares formulas; paired t relationship; zero-difference and exhaustive small sign-flip examples; seed repeatability; R syntax.
- Six comparisons against the earlier R execution: PCA variance, Kruskal-Wallis, pairwise PERMANOVA, PCA scores/loadings and distances to group medians. Shared numerical columns agree within 1e-10 tolerance.
- Main and paired analyses: 9999 permutations with seed 20260707. Tests, analyses, figure export and tree audit also ran directly from the repository path with spaces/Chinese characters under a supported UTF-8 locale. All 28 CSV outputs reproduced the staged run within 1e-12 tolerance; both repository input CSVs are byte-identical.
- New Fig. 3 PNG and PDF visually inspected. PDF export changed to Cairo to preserve micrometre symbols; negative pattern-factor observations are displayed.
- Three earlier local R figure scripts were archived byte-for-byte with SHA-256 hashes and BOM-aware syntax checks. They are historical sources, not supported runtime entry points.

## R results

| Analysis | F | R2 | p |
|---|---:|---:|---:|
| One-factor SpA | 12.043584 | 0.434626 | 0.0001 |
| One-factor PERMDISP | 0.115535 | — | 0.9566 |
| Paired SpA | 13.483215 | 0.462548 | 0.0001 |
| Paired Region | 0.590060 | 0.011664 | 0.6092 |
| Paired SpA x Region | 3.401235 | 0.178375 | 0.0036 |
| Average-space PERMDISP | 0.157812 | — | 0.9288 |
| Difference-space PERMDISP | 0.322911 | — | 0.8121 |

Paired F/R2 reproduce the manuscript snapshot to its displayed precision. The effects use separate spaces. R2 values are not additive.

## Items still requiring scientific reconciliation

1. Ten of 19 displayed p-value entries differ from Supplement S11-S14 in the saved manuscript snapshot. See manuscript_reconciliation.csv for exact new values and R's three-decimal display. The Methods R version is also now 4.5.3 instead of the saved 4.4.2. The Word manuscript has not been rewritten. Do not mix the old manuscript p-values with these new output tables.
2. The paired script is newly implemented from the recorded method; an original standalone R script was not found.
3. Tree auditing checked a 218-taxon x 6133-site alignment and eight primary trees with matching tip sets and finite nonnegative branch lengths. Of 13 tree assets, seven have STRUCTURE_PASS and six REVIEW status. Two MCC files have labels up to 30001; three Mesquite files need verified taxon mapping; one Mesquite file could not be parsed by the selected ape reader. REVIEW is not a pass.
4. No full IQ-TREE/MrBayes rerun, MCMC convergence certification, MCC support correction or ancestral-state reconstruction validation was performed. Source histories and large sampled trees remain outside Git under existing ignore rules.
5. Package versions are recorded; the installer is not a dependency lock. Exact numerical reproduction should use the recorded environment.
6. External MorphoSource accessibility/cloud sync and manuscript figure placement were not checked.

## Evidence files

- regression_checks_20261008.txt
- r_baseline_comparison.csv
- manuscript_reconciliation.csv
- R_environment_20261008.csv
- source_inventory_20261008.csv
- ../results/r_onefactor_20261008/
- ../results/r_paired_20261008/
- ../results/fig3_R_20261008/
- ../results/r_phylogeny_validation_20261008/

The repository status, final commit and remote verification are recorded in the local audit/Obsidian follow-up so the report does not need to embed its own commit ID.
