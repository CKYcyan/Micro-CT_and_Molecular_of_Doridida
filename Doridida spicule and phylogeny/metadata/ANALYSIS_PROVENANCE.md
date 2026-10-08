# R analysis provenance — 2026-10-08

The manuscript workflow is R-only. Historical non-R scripts/results remain unchanged for traceability and are excluded from the release's numerical validation.

## Inputs and lineage

The previously edited Morphological_Analysis/analysis.csv replaces B-E with SpA1-SpA4. This release synchronizes the canonical CSV to those bytes. Sample IDs and 12 measurement columns are unchanged. There are 51 independent specimens; the 102 region rows are paired measurements.

The one-factor and Fig. 3 scripts revise tracked R code. Earlier scripts found in the local manuscript's 0.Documant/Publich/ folder are archived byte-for-byte, with SHA-256 hashes in source_inventory_20261008.csv. Local audit copies and R session histories are not formal analysis sources.

No original standalone paired-analysis R script was located. The new script implements Methods paragraph 27 and Supplement S13-S14 in the saved 2026-09-23 manuscript snapshot. The 2026-09-04 handoff described earlier results as a direct-permutation audit awaiting R/vegan verification. This release identifies the paired R workflow as newly implemented.

## Paired calculation

Center and scale each corresponding metric over the pooled 102 mantle/foot measurements. Define A_i=(M_i+F_i)/2 and D_i=M_i-F_i in standardized space.

- SpA: Euclidean PERMANOVA of A_i, permuting group labels across specimens.
- Region: SS_region=n*sum(mean(D)^2); SS_residual=sum((D-mean(D))^2); F=SS_region/(SS_residual/(n-1)); R2=SS_region/(SS_region+SS_residual). Flip one sign per complete six-variable specimen vector. Under the null, this assumes within-specimen region swaps are exchangeable.
- SpA x Region: Euclidean PERMANOVA of D_i, permuting specimen group labels.
- Dispersion: separate spatial-median PERMDISP tests in A and D spaces, bias.adjust=FALSE, with vegan residual permutations.

The three R2 values use separate reference spaces and must not be summed as a partition of one model. Unbalanced group sizes and dispersion results should be considered when interpreting group effects.

RNG: Mersenne-Twister/Inversion/Rejection; seed 20260707. One-factor order: global PERMANOVA, six pairs in combn order, PERMDISP. Paired order: SpA, Region sign flips, interaction, average dispersion, difference dispersion. Each uses 9999 random permutations with the plus-one p-value convention. Seed and order were recorded for reproducibility, not tuned to match older p-values.

## Two-sided inference

Dunn's default half-tail representation uses a corresponding half-alpha threshold. The revised script explicitly uses altp=TRUE and base R Holm adjustment within each variable's six contrasts. This removes ambiguity without implying that the earlier half-tail results were necessarily interpreted incorrectly. All 72 contrasts were checked independently using ranks, tie-corrected variance, normal tails and Holm adjustment. Figure letters use the calculated results.

References: [dunn.test manual](https://cran.r-project.org/web/packages/dunn.test/dunn.test.pdf), [vegan adonis2](https://vegandevs.github.io/vegan/reference/adonis.html), [vegan betadisper](https://vegandevs.github.io/vegan/reference/betadisper.html).

## Evidence boundaries

Paired F and R2 agree with the manuscript's displayed precision. Some p-values and the R version differ; see manuscript_reconciliation.csv. The Word manuscript was not modified.

Eight primary trees have 218 unique alignment-matched tips and finite nonnegative branch lengths. This does not certify inference, convergence, support values or ancestral states. MCC labels reach 30001 and cannot be read directly as posterior probabilities. Some Mesquite files use numeric or renamed taxa; one is not readable with the selected ape reader. Mapping and provenance remain unresolved. No speculative relabeling or support normalization was performed.

Large sampled MrBayes trees, raw volumes, backups and .Rhistory remain excluded. An exploratory session history does not substitute for a validated MCC or ancestral-state workflow.
