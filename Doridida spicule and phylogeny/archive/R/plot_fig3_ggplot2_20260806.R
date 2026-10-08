#!/usr/bin/env Rscript

# Recreate a Fig. 3-style morphometric boxplot after the R002/R003 recoding.
# This script uses only base R plus ggplot2.

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

default_input <- file.path(
  "G:", "???脩垢蝖祉?", "0.RPG_System_Center", "Raw", "0.Publish",
  "1.Doridida spicule and phylogeny", "MicroCT_data", "analysis.csv"
)
default_outdir <- file.path(
  "G:", "???脩垢蝖祉?", "0.RPG_System_Center", "Raw", "0.Publish",
  "1.Doridida spicule and phylogeny", "0.Documant", "Publich"
)

input_csv <- if (length(args) >= 1) args[[1]] else default_input
outdir <- if (length(args) >= 2) args[[2]] else default_outdir

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

raw <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
raw <- raw[raw$Network %in% c("B", "C", "D", "E"), ]
raw$Network <- factor(
  raw$Network,
  levels = c("B", "C", "D", "E"),
  labels = c("Separate", "Cobweb", "Trabecular", "Sheets")
)

spec <- data.frame(
  Column = c(
    "M.BV/TV", "F.BV/TV",
    "M.Tb.Th", "F.Tb.Th",
    "M.Tb.N", "F.Tb.N",
    "M.Tb.Pf", "F.Tb.Pf",
    "M.SMI", "F.SMI",
    "M.BS/BV", "F.BS/BV"
  ),
  Variable = rep(c("SV/TV", "S.Th", "S.N", "S.Pf", "SMI", "SS/SV"), each = 2),
  Region = rep(c("Mantle", "Foot"), times = 6),
  stringsAsFactors = FALSE
)

long <- do.call(rbind, lapply(seq_len(nrow(spec)), function(i) {
  data.frame(
    Sample = raw$Sample,
    Network = raw$Network,
    Region = spec$Region[[i]],
    Variable = spec$Variable[[i]],
    Value = raw[[spec$Column[[i]]]],
    stringsAsFactors = FALSE
  )
}))

long$Region <- factor(long$Region, levels = c("Mantle", "Foot"))
long$Variable <- factor(
  long$Variable,
  levels = c("SV/TV", "S.Th", "S.N", "S.Pf", "SMI", "SS/SV")
)

letters <- data.frame(
  Region = rep(c(rep("Mantle", 6), rep("Foot", 6)), each = 4),
  Variable = rep(rep(c("SV/TV", "S.Th", "S.N", "S.Pf", "SMI", "SS/SV"), each = 4), 2),
  Network = rep(c("Separate", "Cobweb", "Trabecular", "Sheets"), times = 12),
  Letter = c(
    "a", "ab", "b", "b",
    "a", "ab", "b", "ab",
    "a", "ab", "b", "b",
    "a", "ab", "c", "bc",
    "a", "a", "b", "b",
    "a", "ab", "b", "ab",
    "a", "b", "b", "b",
    "a", "b", "b", "ab",
    "a", "b", "b", "b",
    "a", "b", "b", "b",
    "a", "b", "b", "b",
    "a", "b", "b", "ab"
  ),
  stringsAsFactors = FALSE
)
letters$Region <- factor(letters$Region, levels = c("Mantle", "Foot"))
letters$Variable <- factor(letters$Variable, levels = levels(long$Variable))
letters$Network <- factor(letters$Network, levels = levels(long$Network))

max_by_group <- aggregate(
  Value ~ Variable + Network + Region,
  data = long,
  FUN = max,
  na.rm = TRUE
)
range_by_variable <- aggregate(
  Value ~ Variable,
  data = long,
  FUN = function(x) diff(range(x, na.rm = TRUE))
)
names(range_by_variable)[2] <- "Range"
letters <- merge(letters, max_by_group, by = c("Variable", "Network", "Region"), all.x = TRUE)
letters <- merge(letters, range_by_variable, by = "Variable", all.x = TRUE)
letters$Range[letters$Range == 0 | is.na(letters$Range)] <- 1
letters$Y <- letters$Value + letters$Range * 0.08

p <- ggplot(long, aes(x = Network, y = Value, fill = Region)) +
  geom_boxplot(
    width = 0.62,
    position = position_dodge(width = 0.75),
    outlier.shape = 21,
    outlier.size = 1.4,
    linewidth = 0.45
  ) +
  geom_text(
    data = letters,
    aes(x = Network, y = Y, label = Letter, group = Region),
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    size = 3.2,
    vjust = 0
  ) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Mantle" = "black", "Foot" = "white")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "top",
    legend.key.size = grid::unit(0.55, "cm"),
    panel.spacing = grid::unit(0.7, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  )

png_file <- file.path(outdir, "Fig3_ggplot2_R002_R003_D_preview.png")
pdf_file <- file.path(outdir, "Fig3_ggplot2_R002_R003_D_preview.pdf")

ggsave(png_file, p, width = 10, height = 7, units = "in", dpi = 600, bg = "white")
ggsave(pdf_file, p, width = 10, height = 7, units = "in", bg = "white")

message("Wrote: ", png_file)
message("Wrote: ", pdf_file)

