#!/usr/bin/env Rscript

# Added R implementation of the documented paired design; see metadata/ANALYSIS_PROVENANCE.md.
args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
root <- if (length(script_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), ".."), winslash = "/", mustWork = TRUE)
} else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "workflow_helpers.R"))

require_packages("vegan")
opt <- parse_options(root, "r_paired_20261008")
dat <- read_morphometrics(opt$input)
spaces <- paired_spaces(dat)
group <- dat$Network
dir.create(opt$output, showWarnings = FALSE, recursive = TRUE)
start_rng(opt$seed)

average_distance <- dist(spaces$averages)
difference_distance <- dist(spaces$differences)
spa <- vegan::adonis2(average_distance ~ group, permutations = opt$permutations)
region <- sign_flip_test(spaces$differences, opt$permutations)
interaction <- vegan::adonis2(difference_distance ~ group, permutations = opt$permutations)
tests <- rbind(
  adonis_row(spa, "SpA", "Specimen averages"),
  data.frame(Effect = "Region", Space = "Paired mantle-minus-foot differences",
             Df = 1, Df_residual = nrow(dat) - 1,
             SumOfSqs = unname(region["SumOfSqs"]), SumOfSqs_residual = unname(region["SumOfSqs_residual"]),
             R2 = unname(region["R2"]), F = unname(region["F"]), P = unname(region["P"])),
  adonis_row(interaction, "SpA x Region", "Paired mantle-minus-foot differences"))
tests$Permutation <- c("Specimen-level SpA labels", "Independent sign flip of each specimen vector",
                       "Specimen-level SpA labels")
tests$N_specimens <- nrow(dat)
tests$N_permutations <- opt$permutations
write.csv(tests, file.path(opt$output, "paired_permutation_tests.csv"), row.names = FALSE)

disp_average <- dispersion_fit(average_distance, group, opt$permutations)
disp_difference <- dispersion_fit(difference_distance, group, opt$permutations)
disp <- rbind(data.frame(Space = "Specimen averages", disp_average$table),
              data.frame(Space = "Paired differences", disp_difference$table))
write.csv(disp, file.path(opt$output, "paired_permdisp.csv"), row.names = FALSE)
write.csv(data.frame(Sample = dat$Sample, Network = group, spaces$averages, check.names = FALSE),
          file.path(opt$output, "specimen_averages.csv"), row.names = FALSE)
write.csv(data.frame(Sample = dat$Sample, Network = group, spaces$differences, check.names = FALSE),
          file.path(opt$output, "paired_differences.csv"), row.names = FALSE)
write.csv(data.frame(Variable = metric_names, Center = spaces$center, SD = spaces$scale),
          file.path(opt$output, "pooled_standardization.csv"), row.names = FALSE)
write.csv(data.frame(Sample = dat$Sample, Network = group,
                     AverageDistanceToMedian = disp_average$model$distances,
                     DifferenceDistanceToMedian = disp_difference$model$distances),
          file.path(opt$output, "paired_dispersion_distances.csv"), row.names = FALSE)
write_run_metadata(opt, dat, c("vegan", "permute"), "Paired specimen-average/difference analysis")
capture.output({
  cat("New R implementation of the recorded paired method; no independent treatment of 102 rows.\n")
  cat("Six corresponding variables scaled jointly across 102 region measurements before pairing.\n")
  cat("Region SS = n * sum(mean difference vector squared); residual SS = centered difference SS.\n")
  cat("Region pseudo-F = SS_region / (SS_residual/(n-1)); sign flips preserve each entire vector.\n")
  cat("The three effects use separate reference spaces; their R2 values are not additive.\n")
  cat("PERMDISP uses vegan residual permutations, spatial medians and bias.adjust=FALSE.\n\n")
  print(tests); print(disp)
}, file = file.path(opt$output, "summary.txt"))
message("R paired outputs: ", opt$output)
