#!/usr/bin/env Rscript
# Installs missing packages into the active R library; consult recorded versions
# in results/*/package_versions.csv when reproducing the archived run.
packages <- c("vegan", "dunn.test", "ggplot2", "multcompView", "ape")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
for (package in packages) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Package installation failed: ", package)
  cat(package, as.character(packageVersion(package)), "\n")
}
