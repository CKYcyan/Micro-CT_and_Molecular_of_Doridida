# Shared R-only input, inference and provenance functions.
arrangement_levels <- paste0("SpA", 1:4)
arrangement_names <- setNames(c("Separate", "Cobweb", "Trabecular", "Sheets"), arrangement_levels)
metric_names <- c("BV/TV", "BS/BV", "Tb.Pf", "SMI", "Tb.Th", "Tb.N")
response_names <- c(paste0("M.", metric_names), paste0("F.", metric_names))

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "),
                            ". Run R/install_dependencies.R.", call. = FALSE)
}

parse_options <- function(root, output_name, args = commandArgs(trailingOnly = TRUE)) {
  opt <- list(input = file.path(root, "data", "analysis.csv"),
              output = file.path(root, "results", output_name),
              seed = 20260707L, permutations = 9999L)
  if (length(args) %% 2L) stop("Use --input PATH --output PATH --seed INTEGER --permutations INTEGER.")
  if (length(args)) for (i in seq.int(1L, length(args), by = 2L)) {
    key <- sub("^--", "", args[i])
    if (!startsWith(args[i], "--") || !key %in% names(opt)) stop("Unknown option: ", args[i])
    opt[[key]] <- args[i + 1L]
  }
  for (key in c("seed", "permutations")) {
    value <- suppressWarnings(as.numeric(opt[[key]]))
    minimum <- if (key == "seed") 0 else 1
    if (length(value) != 1L || !is.finite(value) || value < minimum ||
        value != floor(value) || value > .Machine$integer.max) stop("Invalid --", key)
    opt[[key]] <- as.integer(value)
  }
  opt$input <- normalizePath(opt$input, winslash = "/", mustWork = TRUE)
  opt
}

read_morphometrics <- function(path) {
  dat <- read.csv(path, fileEncoding = "UTF-8-BOM", check.names = FALSE,
                  strip.white = TRUE, stringsAsFactors = FALSE)
  expected <- c("Sample", "Network", response_names)
  if (anyDuplicated(names(dat)) || !setequal(names(dat), expected)) {
    stop("Input must contain Sample, Network and the 12 documented morphometric columns.", call. = FALSE)
  }
  dat <- dat[expected]
  dat$Sample <- trimws(dat$Sample)
  if (anyNA(dat$Sample) || any(dat$Sample == "") || anyDuplicated(dat$Sample)) {
    stop("Sample identifiers must be nonempty and unique.", call. = FALSE)
  }
  codes <- trimws(dat$Network)
  legacy <- setNames(arrangement_levels, c("B", "C", "D", "E"))
  old <- codes %in% names(legacy)
  codes[old] <- unname(legacy[codes[old]])
  if (anyNA(codes) || any(!codes %in% arrangement_levels)) {
    stop("Network must use SpA1-SpA4 or the legacy aliases B-E.", call. = FALSE)
  }
  dat$Network <- factor(codes, levels = arrangement_levels)
  if (nrow(dat) < 8L || any(table(dat$Network) < 2L)) {
    stop("Each of SpA1-SpA4 needs at least two complete specimens.", call. = FALSE)
  }
  for (variable in response_names) {
    values <- suppressWarnings(as.numeric(trimws(as.character(dat[[variable]]))))
    if (any(!is.finite(values))) stop("Missing, nonnumeric or nonfinite values in ", variable, call. = FALSE)
    if (sd(values) == 0) stop("Zero-variance response: ", variable, call. = FALSE)
    dat[[variable]] <- values
  }
  rownames(dat) <- dat$Sample
  dat
}

dunn_tables <- function(dat) {
  require_packages("dunn.test")
  do.call(rbind, lapply(response_names, function(variable) {
    invisible(capture.output(fit <- dunn.test::dunn.test(
      dat[[variable]], dat$Network, method = "none", kw = FALSE,
      altp = TRUE, table = FALSE, list = FALSE)))
    groups <- strsplit(fit$comparisons, " - ", fixed = TRUE)
    data.frame(Variable = variable,
               Group1 = vapply(groups, "[", character(1), 1L),
               Group2 = vapply(groups, "[", character(1), 2L),
               Comparison = fit$comparisons, Z = fit$Z, P = fit$altP,
               P_adjusted_holm = p.adjust(fit$altP, method = "holm"),
               row.names = NULL, check.names = FALSE)
  }))
}

figure_letters <- function(dunn) {
  require_packages("multcompView")
  out <- lapply(unique(dunn$Variable), function(variable) {
    rows <- dunn[dunn$Variable == variable, ]
    p <- setNames(rows$P_adjusted_holm, paste(rows$Group1, rows$Group2, sep = "-"))
    compact <- multcompView::multcompLetters(p, threshold = 0.05)$Letters
    data.frame(Variable = variable, Network = arrangement_levels,
               Letter = unname(compact[arrangement_levels]), row.names = NULL)
  })
  out <- do.call(rbind, out)
  if (anyNA(out$Letter)) stop("Cannot construct significance letters for all groups.")
  out
}

paired_spaces <- function(dat) {
  mantle <- as.matrix(dat[paste0("M.", metric_names)])
  foot <- as.matrix(dat[paste0("F.", metric_names)])
  colnames(mantle) <- colnames(foot) <- metric_names
  pooled <- rbind(mantle, foot)
  scaled <- scale(pooled)
  n <- nrow(dat)
  m <- scaled[seq_len(n), , drop = FALSE]
  f <- scaled[n + seq_len(n), , drop = FALSE]
  averages <- (m + f) / 2
  differences <- m - f
  rownames(averages) <- rownames(differences) <- dat$Sample
  list(averages = averages, differences = differences,
       center = attr(scaled, "scaled:center"), scale = attr(scaled, "scaled:scale"),
       scaled = scaled)
}

region_statistic <- function(differences) {
  n <- nrow(differences)
  ss_region <- n * sum(colMeans(differences)^2)
  ss_residual <- sum(sweep(differences, 2, colMeans(differences), "-")^2)
  total <- ss_region + ss_residual
  f <- if (total == 0) 0 else if (ss_residual == 0) Inf else ss_region / (ss_residual / (n - 1))
  c(F = f, R2 = if (total == 0) 0 else ss_region / total,
    SumOfSqs = ss_region, SumOfSqs_residual = ss_residual)
}

sign_flip_test <- function(differences, permutations) {
  observed <- region_statistic(differences)
  # One sign for the whole six-variable vector per specimen preserves pairing.
  exceed <- 0L
  for (i in seq_len(permutations)) {
    signs <- sample(c(-1, 1), nrow(differences), replace = TRUE)
    f <- region_statistic(differences * signs)[["F"]]
    if (f >= observed[["F"]] - sqrt(.Machine$double.eps)) exceed <- exceed + 1L
  }
  c(observed, P = (exceed + 1) / (permutations + 1))
}

adonis_row <- function(fit, effect, space) {
  data.frame(Effect = effect, Space = space, Df = fit$Df[1],
             Df_residual = fit$Df[2], SumOfSqs = fit$SumOfSqs[1],
             SumOfSqs_residual = fit$SumOfSqs[2], R2 = fit$R2[1],
             F = fit$F[1], P = fit[["Pr(>F)"]][1], row.names = NULL)
}

dispersion_fit <- function(distance, group, permutations) {
  dispersion <- vegan::betadisper(distance, group, type = "median", bias.adjust = FALSE)
  test <- vegan::permutest(dispersion, permutations = permutations)
  tab <- test$tab
  list(model = dispersion, test = test,
       table = data.frame(Df = tab$Df[1], Df_residual = tab$Df[2],
                          SumOfSqs = tab[["Sum Sq"]][1], SumOfSqs_residual = tab[["Sum Sq"]][2],
                          F = tab$F[1], P = tab[["Pr(>F)"]][1]))
}

start_rng <- function(seed) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
}

write_run_metadata <- function(opt, dat, packages, analysis) {
  capture.output(sessionInfo(), file = file.path(opt$output, "sessionInfo.txt"))
  versions <- data.frame(Package = packages,
                         Version = vapply(packages, function(p) as.character(packageVersion(p)), character(1)))
  write.csv(versions, file.path(opt$output, "package_versions.csv"), row.names = FALSE)
  details <- data.frame(Analysis = analysis, Input = basename(opt$input),
                        Input_MD5 = unname(tools::md5sum(opt$input)),
                        Specimens = nrow(dat), Seed = opt$seed, Permutations = opt$permutations,
                        R = as.character(getRversion()), RNG = paste(RNGkind(), collapse = "/"))
  write.csv(details, file.path(opt$output, "run_metadata.csv"), row.names = FALSE)
  write.csv(data.frame(Network = names(table(dat$Network)), N = as.integer(table(dat$Network))),
            file.path(opt$output, "group_sizes.csv"), row.names = FALSE)
}
