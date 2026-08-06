#!/usr/bin/env Rscript

# Morphometric analyses for spicule-bearing arrangement types.
# The workflow uses standardised Micro-CT morphometric variables and a
# Euclidean distance matrix for PERMANOVA and PERMDISP.

required_packages <- c("vegan")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    ". Install them before running this script.",
    call. = FALSE
  )
}

get_repo_root <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", command_args, value = TRUE)

  if (length(file_arg) > 0) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

repo_root <- get_repo_root()
data_path <- file.path(repo_root, "data", "analysis.csv")
results_dir <- file.path(repo_root, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(data_path, check.names = FALSE, strip.white = TRUE, stringsAsFactors = FALSE)

required_columns <- c("Sample", "Network")
if (!all(required_columns %in% names(dat))) {
  stop("Input data must contain columns named 'Sample' and 'Network'.", call. = FALSE)
}

response_columns <- setdiff(names(dat), required_columns)
if (length(response_columns) == 0) {
  stop("No morphometric response columns were found.", call. = FALSE)
}

responses <- as.data.frame(lapply(dat[response_columns], function(x) as.numeric(trimws(x))))
rownames(responses) <- dat$Sample
network <- factor(dat$Network)

complete_rows <- complete.cases(responses) & !is.na(network)
responses <- responses[complete_rows, , drop = FALSE]
network <- droplevels(network[complete_rows])

scaled_responses <- scale(responses)
distance_matrix <- dist(scaled_responses, method = "euclidean")

set.seed(20260707)
n_permutations <- 9999

pca <- prcomp(scaled_responses, center = FALSE, scale. = FALSE)
pca_variance <- (pca$sdev^2) / sum(pca$sdev^2)
pca_summary <- data.frame(
  PC = paste0("PC", seq_along(pca_variance)),
  ProportionVariance = pca_variance,
  CumulativeVariance = cumsum(pca_variance),
  row.names = NULL
)

write.csv(
  data.frame(Sample = rownames(scaled_responses), Network = network, pca$x, row.names = NULL),
  file.path(results_dir, "pca_scores.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(Variable = rownames(pca$rotation), pca$rotation, row.names = NULL),
  file.path(results_dir, "pca_loadings.csv"),
  row.names = FALSE
)
write.csv(pca_summary, file.path(results_dir, "pca_variance.csv"), row.names = FALSE)

permanova <- vegan::adonis2(distance_matrix ~ network, permutations = n_permutations)

pairwise_permanova <- function(distance_matrix, grouping, permutations = 9999) {
  groups <- levels(grouping)
  pairs <- combn(groups, 2, simplify = FALSE)

  pair_results <- lapply(pairs, function(pair) {
    keep <- grouping %in% pair
    pair_distance <- as.dist(as.matrix(distance_matrix)[keep, keep])
    pair_group <- droplevels(grouping[keep])
    fit <- vegan::adonis2(pair_distance ~ pair_group, permutations = permutations)

    data.frame(
      Group1 = pair[1],
      Group2 = pair[2],
      Df = fit$Df[1],
      SumOfSqs = fit$SumOfSqs[1],
      R2 = fit$R2[1],
      F = fit$F[1],
      P = fit$`Pr(>F)`[1],
      row.names = NULL
    )
  })

  pair_results <- do.call(rbind, pair_results)
  pair_results$P_adjusted_holm <- p.adjust(pair_results$P, method = "holm")
  pair_results
}

pairwise_results <- pairwise_permanova(distance_matrix, network, n_permutations)
write.csv(pairwise_results, file.path(results_dir, "pairwise_permanova.csv"), row.names = FALSE)

# PERMDISP using spatial medians follows the manuscript notes and uses the same
# standardised Euclidean distance matrix as PERMANOVA.
dispersion <- vegan::betadisper(distance_matrix, network, type = "median")
permdisp <- vegan::permutest(dispersion, permutations = n_permutations)
write.csv(
  data.frame(
    Sample = rownames(scaled_responses),
    Network = network,
    DistanceToMedian = dispersion$distances,
    row.names = NULL
  ),
  file.path(results_dir, "permdisp_distances_to_median.csv"),
  row.names = FALSE
)

kruskal_results <- do.call(
  rbind,
  lapply(names(responses), function(variable) {
    fit <- kruskal.test(responses[[variable]] ~ network)
    data.frame(
      Variable = variable,
      Statistic = unname(fit$statistic),
      Df = unname(fit$parameter),
      P = fit$p.value,
      row.names = NULL
    )
  })
)
kruskal_results$P_adjusted_holm <- p.adjust(kruskal_results$P, method = "holm")
write.csv(kruskal_results, file.path(results_dir, "kruskal_wallis_by_variable.csv"), row.names = FALSE)

if (requireNamespace("dunn.test", quietly = TRUE)) {
  dunn_results <- lapply(names(responses), function(variable) {
    fit <- dunn.test::dunn.test(
      x = responses[[variable]],
      g = network,
      method = "holm",
      kw = FALSE,
      list = TRUE
    )

    data.frame(
      Variable = variable,
      Comparison = fit$comparisons,
      Z = fit$Z,
      P = fit$P,
      P_adjusted_holm = fit$P.adjusted,
      row.names = NULL
    )
  })
  write.csv(do.call(rbind, dunn_results), file.path(results_dir, "dunn_tests_by_variable.csv"), row.names = FALSE)
}

sink(file.path(results_dir, "microct_spicule_arrangement_analysis_summary.txt"))
cat("Input file:", data_path, "\n")
cat("Samples included:", nrow(responses), "\n")
cat("Variables included:", length(response_columns), "\n")
cat("Network group sizes:\n")
print(table(network))

cat("\nPCA variance summary:\n")
print(pca_summary)

cat("\nPERMANOVA on standardised Euclidean distance matrix:\n")
print(permanova)

cat("\nPairwise PERMANOVA comparisons:\n")
print(pairwise_results)

cat("\nPERMDISP using spatial medians on the same distance matrix:\n")
print(permdisp)

cat("\nKruskal-Wallis tests by variable:\n")
print(kruskal_results)

cat("\nSession information:\n")
print(sessionInfo())
sink()

message("Analysis complete. Outputs written to: ", results_dir)
