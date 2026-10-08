#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

input_csv <- file.path(
  "1.Doridida spicule and phylogeny",
  "MicroCT_data",
  "analysis.csv"
)
outdir <- file.path(
  "1.Doridida spicule and phylogeny",
  "0.Documant",
  "Publich"
)

if (!file.exists(input_csv)) {
  stop("Input CSV not found. Run this script from the 0.Publish directory.")
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
dat <- dat[dat$Network %in% c("B", "C", "D", "E"), ]
dat$NetworkCode <- factor(dat$Network, levels = c("B", "C", "D", "E"))
dat$Arrangement <- factor(
  dat$Network,
  levels = c("B", "C", "D", "E"),
  labels = c("Separate", "Cobweb", "Trabecular", "Sheets")
)

arrangement_cols <- c(
  "Separate" = "orange",
  "Cobweb" = "green",
  "Trabecular" = "blue",
  "Sheets" = "pink"
)

variable_specs <- list(
  list(panel = "A", label = "SV/TV", m = "M.BV/TV", f = "F.BV/TV"),
  list(panel = "B", label = "S.Th", m = "M.Tb.Th", f = "F.Tb.Th"),
  list(panel = "C", label = "S.N", m = "M.Tb.N", f = "F.Tb.N"),
  list(panel = "D", label = "S.Pf", m = "M.Tb.Pf", f = "F.Tb.Pf", ymin = 0),
  list(panel = "E", label = "SMI", m = "M.SMI", f = "F.SMI", smi_axis = TRUE),
  list(panel = "F", label = "SS/SV", m = "M.BS/BV", f = "F.BS/BV")
)

dunn_letters <- data.frame(
  Region = rep(c(rep("Mantle", 6), rep("Foot", 6)), each = 4),
  Variable = rep(rep(c("SV/TV", "S.Th", "S.N", "S.Pf", "SMI", "SS/SV"), each = 4), 2),
  Arrangement = rep(c("Separate", "Cobweb", "Trabecular", "Sheets"), times = 12),
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
dunn_letters$Arrangement <- factor(
  dunn_letters$Arrangement,
  levels = c("Separate", "Cobweb", "Trabecular", "Sheets")
)
dunn_letters$Region <- factor(dunn_letters$Region, levels = c("Mantle", "Foot"))

long_for_spec <- function(spec) {
  rbind(
    data.frame(
      Sample = dat$Sample,
      NetworkCode = dat$NetworkCode,
      Arrangement = dat$Arrangement,
      Region = "Mantle",
      Variable = spec$label,
      Value = dat[[spec$m]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      Sample = dat$Sample,
      NetworkCode = dat$NetworkCode,
      Arrangement = dat$Arrangement,
      Region = "Foot",
      Variable = spec$label,
      Value = dat[[spec$f]],
      stringsAsFactors = FALSE
    )
  )
}

make_panel <- function(spec, show_x = FALSE) {
  df <- long_for_spec(spec)
  df$Arrangement <- factor(df$Arrangement, levels = names(arrangement_cols))
  df$Region <- factor(df$Region, levels = c("Mantle", "Foot"))
  df$Variable <- spec$label

  letters_df <- dunn_letters[dunn_letters$Variable == spec$label, ]
  max_df <- aggregate(
    Value ~ Arrangement + Region,
    data = df,
    FUN = max,
    na.rm = TRUE
  )
  letters_df <- merge(letters_df, max_df, by = c("Arrangement", "Region"), all.x = TRUE)
  panel_range <- diff(range(df$Value, na.rm = TRUE))
  if (panel_range == 0) {
    panel_range <- 1
  }
  letters_df$Y <- letters_df$Value + panel_range * 0.08

  ymin <- min(df$Value, na.rm = TRUE)
  ymax <- max(c(df$Value, letters_df$Y), na.rm = TRUE)
  if (!is.null(spec$ymin)) {
    ymin <- spec$ymin
  } else {
    ymin <- ymin - panel_range * 0.10
  }
  ymax <- ymax + panel_range * 0.12

  y_scale <- scale_y_continuous()
  if (!is.null(spec$smi_axis) && isTRUE(spec$smi_axis)) {
    y_scale <- scale_y_continuous(labels = function(x) sprintf("%.1f", x))
  }

  ggplot(df, aes(x = Arrangement, y = Value)) +
    geom_boxplot(
      aes(fill = Arrangement, alpha = Region),
      width = 0.84,
      position = position_dodge2(width = 0.78, preserve = "single", padding = 0.04),
      outlier.shape = 21,
      outlier.size = 1.2,
      outlier.stroke = 0.30,
      linewidth = 0.45,
      color = "black"
    ) +
    geom_text(
      data = letters_df,
      aes(x = Arrangement, y = Y, label = Letter, group = Region),
      position = position_dodge2(width = 0.78, preserve = "single", padding = 0.04),
      inherit.aes = FALSE,
      size = 3.25,
      family = "sans",
      vjust = 0
    ) +
    scale_fill_manual(values = arrangement_cols, name = "Spicule arrangement") +
    scale_alpha_manual(
      values = c("Mantle" = 0.95, "Foot" = 0.40),
      name = "Region"
    ) +
    y_scale +
    coord_cartesian(ylim = c(ymin, ymax), clip = "off") +
    labs(
      x = NULL,
      y = spec$label,
      title = paste0("(", spec$panel, ") ", spec$label)
    ) +
    theme_classic(base_family = "sans", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      axis.text.x = if (show_x) {
        element_text(size = 9, angle = 25, hjust = 1)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) element_line(linewidth = 0.35) else element_blank(),
      legend.position = "none",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      legend.key.width = unit(0.8, "lines"),
      legend.key.height = unit(0.8, "lines"),
      plot.margin = margin(6, 7, 6, 7)
    )
}

plots <- list(
  make_panel(variable_specs[[1]], show_x = FALSE),
  make_panel(variable_specs[[2]], show_x = FALSE),
  make_panel(variable_specs[[3]], show_x = FALSE),
  make_panel(variable_specs[[4]], show_x = TRUE),
  make_panel(variable_specs[[5]], show_x = TRUE),
  make_panel(variable_specs[[6]], show_x = TRUE)
)

png_file <- file.path(outdir, "Fig3_ggplot2_R002_R003_D_final.png")
pdf_file <- file.path(outdir, "Fig3_ggplot2_R002_R003_D_final.pdf")

draw_manual_legend <- function() {
  pushViewport(viewport(xscale = c(0, 1), yscale = c(0, 1)))
  grid.text(
    "Spicule arrangement",
    x = 0.06, y = 0.62,
    just = c("left", "center"),
    gp = gpar(fontfamily = "sans", fontsize = 10, fontface = "bold")
  )
  x0 <- 0.25
  for (i in seq_along(arrangement_cols)) {
    xx <- x0 + (i - 1) * 0.135
    grid.rect(
      x = xx, y = 0.62,
      width = 0.022, height = 0.22,
      gp = gpar(fill = arrangement_cols[[i]], col = "black", alpha = 0.95, lwd = 0.5)
    )
    grid.text(
      names(arrangement_cols)[[i]],
      x = xx + 0.018, y = 0.62,
      just = c("left", "center"),
      gp = gpar(fontfamily = "sans", fontsize = 9)
    )
  }
  grid.rect(
    x = 0.83, y = 0.62,
    width = 0.022, height = 0.22,
    gp = gpar(fill = "gray35", col = "black", alpha = 0.95, lwd = 0.5)
  )
  grid.text(
    "Mantle",
    x = 0.85, y = 0.62,
    just = c("left", "center"),
    gp = gpar(fontfamily = "sans", fontsize = 9)
  )
  grid.rect(
    x = 0.925, y = 0.62,
    width = 0.022, height = 0.22,
    gp = gpar(fill = "gray35", col = "black", alpha = 0.40, lwd = 0.5)
  )
  grid.text(
    "Foot",
    x = 0.945, y = 0.62,
    just = c("left", "center"),
    gp = gpar(fontfamily = "sans", fontsize = 9)
  )
  popViewport()
}

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    3, 3,
    heights = unit(c(0.10, 0.45, 0.45), "npc")
  )))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  draw_manual_legend()
  popViewport()
  for (i in seq_along(plots)) {
    row <- ceiling(i / 3) + 1
    col <- ((i - 1) %% 3) + 1
    print(plots[[i]], vp = viewport(layout.pos.row = row, layout.pos.col = col))
  }
  popViewport()
}

png(filename = png_file, width = 3300, height = 2300, res = 300, type = "cairo")
draw_figure()
dev.off()

pdf(file = pdf_file, width = 11, height = 7.67, useDingbats = FALSE)
draw_figure()
dev.off()

message("Wrote: ", png_file)
message("Wrote: ", pdf_file)


