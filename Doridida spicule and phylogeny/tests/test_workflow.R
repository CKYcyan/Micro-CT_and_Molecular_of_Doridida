args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
root <- if (length(script_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), ".."), winslash = "/", mustWork = TRUE)
} else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "workflow_helpers.R"))


require_packages(c("vegan", "dunn.test", "multcompView", "ggplot2", "ape"))
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAIL: ", label, call. = FALSE)
  checks <<- checks + 1L
  cat("PASS: ", label, "\n", sep = "")
}
equal <- function(a, b, tolerance = 1e-10) isTRUE(all.equal(unname(a), unname(b), tolerance = tolerance, check.attributes = FALSE))
dat <- read_morphometrics(file.path(root, "data", "analysis.csv"))
check(nrow(dat) == 51L && identical(as.integer(table(dat$Network)), c(21L, 11L, 15L, 4L)), "specimen count and group sizes")
check(all(dat$Network[dat$Sample %in% c("R002", "R003")] == "SpA3"), "R002/R003 coding")
scratch <- tempfile("doridida-r-tests-")
dir.create(scratch)
# Fixtures remain in the R temporary directory; no project files are deleted.
write_fixture <- function(x, name, bom = FALSE) {
  path <- file.path(scratch, name)
  text <- capture.output(write.csv(x, row.names = FALSE, quote = TRUE))
  bytes <- charToRaw(paste0(paste(text, collapse = "\n"), "\n"))
  if (bom) bytes <- c(as.raw(c(239, 187, 191)), bytes)
  writeBin(bytes, path)
  path
}
canonical <- dat
canonical$Network <- as.character(canonical$Network)
legacy <- canonical
legacy$Network <- setNames(c("B", "C", "D", "E"), arrangement_levels)[legacy$Network]
paths <- c(write_fixture(canonical, "spa.csv"), write_fixture(canonical, "spa_bom.csv", TRUE),
           write_fixture(legacy, "legacy.csv"), write_fixture(legacy, "legacy_bom.csv", TRUE))
for (path in paths) check(equal(read_morphometrics(path), dat), paste("BOM/label parity:", basename(path)))
reject <- function(x, name) {
  path <- write_fixture(x, paste0(name, ".csv"))
  check(inherits(tryCatch(read_morphometrics(path), error = identity), "error"), paste("reject", name))
}
bad <- canonical; bad$Sample[2] <- bad$Sample[1]; reject(bad, "duplicate_sample")
bad <- canonical; bad$Sample[1] <- ""; reject(bad, "empty_sample")
bad <- canonical; bad$Network[1] <- "SpA0"; reject(bad, "unknown_group")
bad <- canonical; bad[[response_names[1]]][1] <- Inf; reject(bad, "infinite_value")
bad <- canonical; bad[[response_names[1]]][1] <- NA; reject(bad, "missing_value")
bad <- canonical; bad[[response_names[1]]] <- 1; reject(bad, "constant_response")
bad <- canonical; bad[[response_names[1]]] <- NULL; reject(bad, "missing_column")
bad <- canonical; bad[[response_names[1]]][1] <- "invalid"; reject(bad, "nonnumeric_value")
for (args in list(c("--seed", "-1"), c("--permutations", "0"), c("--other", "1"), "--input")) {
  check(inherits(tryCatch(parse_options(root, "unused", args), error = identity), "error"),
        paste("reject CLI:", paste(args, collapse = " ")))
}

dunn <- dunn_tables(dat)
letters <- figure_letters(dunn)
check(nrow(dunn) == 72L && nrow(letters) == 48L, "Dunn and compact-letter dimensions")
for (variable in response_names) {
  x <- dat[[variable]]
  n <- length(x)
  ranks <- rank(x)
  ties <- as.numeric(table(x))
  variance <- n * (n + 1) / 12 - sum(ties^3 - ties) / (12 * (n - 1))
  means <- tapply(ranks, dat$Network, mean)
  sizes <- table(dat$Network)
  rows <- dunn[dunn$Variable == variable, ]
  z <- (means[rows$Group1] - means[rows$Group2]) /
       sqrt(variance * (1 / sizes[rows$Group1] + 1 / sizes[rows$Group2]))
  independent_p <- 2 * pnorm(-abs(as.numeric(z)))
  check(equal(abs(rows$Z), as.numeric(abs(z))), paste("rank/tie-corrected Dunn Z:", variable))
  check(equal(rows$P, independent_p) && equal(rows$P_adjusted_holm, p.adjust(independent_p, "holm")),
        paste("two-sided P and six-comparison Holm family:", variable))
  cld <- subset(letters, Variable == variable)
  bygroup <- setNames(cld$Letter, cld$Network)
  matches <- vapply(seq_len(nrow(rows)), function(i) {
    shared <- length(intersect(strsplit(bygroup[rows$Group1[i]], "")[[1]],
                               strsplit(bygroup[rows$Group2[i]], "")[[1]])) > 0
    shared == (rows$P_adjusted_holm[i] >= 0.05)
  }, logical(1))
  check(all(matches), paste("letters match six pairwise decisions:", variable))
}

spaces <- paired_spaces(dat)
check(max(abs(colMeans(spaces$scaled))) < 1e-12 &&
      max(abs(apply(spaces$scaled, 2, sd) - 1)) < 1e-12, "joint standardization of 102 measurements")
n <- nrow(dat)
check(equal(spaces$averages + spaces$differences / 2, spaces$scaled[seq_len(n), ]) &&
      equal(spaces$averages - spaces$differences / 2, spaces$scaled[n + seq_len(n), ]),
      "average/difference reconstruct original paired observations")
group_statistic <- function(x, g) {
  grand <- colMeans(x)
  ss_within <- sum(vapply(split(seq_len(nrow(x)), g), function(i) {
    sum(sweep(x[i, , drop = FALSE], 2, colMeans(x[i, , drop = FALSE]), "-")^2)
  }, numeric(1)))
  total <- sum(sweep(x, 2, grand, "-")^2)
  between <- total - ss_within
  c(F = (between / (nlevels(g) - 1)) / (ss_within / (length(g) - nlevels(g))),
    R2 = between / total)
}
for (space in c("averages", "differences")) {
  g <- dat$Network
  fit <- vegan::adonis2(dist(spaces[[space]]) ~ g, permutations = 0)
  direct <- group_statistic(spaces[[space]], g)
  check(equal(c(fit$F[1], fit$R2[1]), direct), paste("independent sums of squares:", space))
}
check(equal(region_statistic(matrix(spaces$differences[, 1], ncol = 1))[["F"]],
            unname(t.test(spaces$differences[, 1])$statistic)^2),
      "univariate Region F equals paired t squared")
zero <- matrix(0, nrow = 4, ncol = 2)
check(region_statistic(zero)[["F"]] == 0 && sign_flip_test(zero, 19)[["P"]] == 1,
      "zero-difference edge case")
start_rng(15); a <- sign_flip_test(spaces$differences, 199)
start_rng(15); b <- sign_flip_test(spaces$differences, 199)
check(identical(a, b), "sign-flip seed repeatability")
# Exhaust all 2^3 sign configurations independently; compare the statistic against
# a formula based on raw sums and check whole-vector sign changes.
toy <- matrix(c(1, 2, 4, 3, -1, 2), nrow = 3)
signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), 3)))
for (i in seq_len(nrow(signs))) {
  flipped <- toy * signs[i, ]
  numerator <- sum(colSums(flipped)^2) / 3
  denominator <- (sum(flipped^2) - numerator) / 2
  check(equal(region_statistic(flipped)[["F"]], numerator / denominator),
        paste("exhaustive sign-flip formula", i))
}
check(identical(dim(spaces$differences), c(51L, 6L)), "51 independent specimen vectors")
for (path in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) {
  parse(file = path)
  check(TRUE, paste("R parse:", basename(path)))
}
cat("SUCCESS: ", checks, " checks passed.\n", sep = "")
