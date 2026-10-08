#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
root <- if (length(script_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), ".."), winslash = "/", mustWork = TRUE)
} else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "workflow_helpers.R"))

require_packages(c("vegan", "dunn.test"))
opt <- parse_options(root, "r_onefactor_20261008")
dat <- read_morphometrics(opt$input)
responses <- dat[response_names]
scaled_responses <- scale(responses)
distance_matrix <- dist(scaled_responses, method = "euclidean")
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
start_rng(opt$seed)

pca <- prcomp(scaled_responses, center = FALSE, scale. = FALSE)
variance <- pca$sdev^2 / sum(pca$sdev^2)
write.csv(data.frame(Sample = dat$Sample, Network = dat$Network, pca$x, check.names = FALSE),
          file.path(opt$output, "pca_scores.csv"), row.names = FALSE)
write.csv(data.frame(Variable = rownames(pca$rotation), pca$rotation, check.names = FALSE),
          file.path(opt$output, "pca_loadings.csv"), row.names = FALSE)
write.csv(data.frame(PC = paste0("PC", seq_along(variance)), StandardDeviation = pca$sdev,
                     ProportionVariance = variance, CumulativeVariance = cumsum(variance)),
          file.path(opt$output, "pca_variance.csv"), row.names = FALSE)

group <- dat$Network
permanova <- vegan::adonis2(distance_matrix ~ group, permutations = opt$permutations)
write.csv(adonis_row(permanova, "SpA", "12 standardized region-specific variables"),
          file.path(opt$output, "permanova_global.csv"), row.names = FALSE)

pairwise <- do.call(rbind, lapply(combn(levels(group), 2, simplify = FALSE), function(pair) {
  keep <- group %in% pair
  pair_distance <- as.dist(as.matrix(distance_matrix)[keep, keep])
  pair_group <- droplevels(group[keep])
  fit <- vegan::adonis2(pair_distance ~ pair_group, permutations = opt$permutations)
  data.frame(Group1 = pair[1], Group2 = pair[2], Df = fit$Df[1], Df_residual = fit$Df[2],
             SumOfSqs = fit$SumOfSqs[1], R2 = fit$R2[1], F = fit$F[1], P = fit[["Pr(>F)"]][1])
}))
pairwise$P_adjusted_holm <- p.adjust(pairwise$P, method = "holm")
write.csv(pairwise, file.path(opt$output, "pairwise_permanova.csv"), row.names = FALSE)

dispersion <- dispersion_fit(distance_matrix, group, opt$permutations)
write.csv(dispersion$table, file.path(opt$output, "permdisp_global.csv"), row.names = FALSE)
write.csv(data.frame(Sample = dat$Sample, Network = group,
                     DistanceToMedian = dispersion$model$distances),
          file.path(opt$output, "permdisp_distances_to_median.csv"), row.names = FALSE)

kruskal <- do.call(rbind, lapply(response_names, function(variable) {
  fit <- kruskal.test(dat[[variable]] ~ group)
  data.frame(Variable = variable, Statistic = unname(fit$statistic),
             Df = unname(fit$parameter), P = fit$p.value)
}))
kruskal$P_adjusted_holm <- p.adjust(kruskal$P, method = "holm")
write.csv(kruskal, file.path(opt$output, "kruskal_wallis_by_variable.csv"), row.names = FALSE)
write.csv(dunn_tables(dat), file.path(opt$output, "dunn_tests_by_variable.csv"), row.names = FALSE)

write_run_metadata(opt, dat, c("vegan", "permute", "dunn.test"), "One-factor morphometrics")
capture.output({
  cat("R-only one-factor analysis; specimen is the independent unit.\n")
  cat("Input:", basename(opt$input), "\nSamples:", nrow(dat), "\nPermutations:", opt$permutations, "\n")
  cat("Dunn P: two-sided; Holm adjustment within six contrasts for each variable.\n")
  cat("Kruskal-Wallis Holm family: 12 variables; pairwise PERMANOVA Holm family: six contrasts.\n\n")
  print(table(group)); print(permanova); print(pairwise); print(dispersion$test); print(kruskal)
}, file = file.path(opt$output, "summary.txt"))
message("R one-factor outputs: ", opt$output)
