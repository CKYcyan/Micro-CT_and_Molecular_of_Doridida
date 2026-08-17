#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

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
input_csv <- file.path(repo_root, "data", "analysis.csv")
outdir <- file.path(repo_root, "results", "fig3_20260807")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
dat <- dat[dat$Network %in% c("B", "C", "D", "E"), ]
dat$NetworkCode <- factor(dat$Network, levels = c("B", "C", "D", "E"))
dat$Arrangement <- factor(
  dat$Network,
  levels = c("B", "C", "D", "E"),
  labels = c("Separate", "Cobweb", "Trabecular", "Sheets")
)

# Use the previous Foot-like pale tones as the base colors for both regions.
arrangement_cols <- c(
  "Separate" = "#FFDB99",
  "Cobweb" = "#99FF99",
  "Trabecular" = "#9999FF",
  "Sheets" = "#FFD9E0"
)

x_axis_labels <- c("(B) Separate", "(C) Cobweb", "(D) Trabecular", "(E) Sheets")

variable_specs <- list(
  list(
    panel = "A", key = "SV/TV",
    axis_label = "Percent spicule volume (%)",
    m = "M.BV/TV", f = "F.BV/TV"
  ),
  list(
    panel = "B", key = "S.Th",
    axis_label = expression("Spicule arrangement thickness (" * mu * "m)"),
    m = "M.Tb.Th", f = "F.Tb.Th"
  ),
  list(
    panel = "C", key = "S.N",
    axis_label = expression("Spicule number (" * mu * "m"^-1 * ")"),
    m = "M.Tb.N", f = "F.Tb.N"
  ),
  list(
    panel = "D", key = "S.Pf",
    axis_label = expression("Spicule pattern factor (" * mu * "m"^-1 * ")"),
    m = "M.Tb.Pf", f = "F.Tb.Pf", ymin = 0
  ),
  list(
    panel = "E", key = "SMI",
    axis_label = "Structure model index (dimensionless)",
    m = "M.SMI", f = "F.SMI", smi_axis = TRUE, tick_by = 1.0
  ),
  list(
    panel = "F", key = "SS/SV",
    axis_label = expression("Spicule surface/volume ratio (" * mu * "m"^-1 * ")"),
    m = "M.BS/BV", f = "F.BS/BV"
  )
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
      Variable = spec$key,
      Value = dat[[spec$m]],
      stringsAsFactors = FALSE
    ),
    data.frame(
      Sample = dat$Sample,
      NetworkCode = dat$NetworkCode,
      Arrangement = dat$Arrangement,
      Region = "Foot",
      Variable = spec$key,
      Value = dat[[spec$f]],
      stringsAsFactors = FALSE
    )
  )
}

box_stats_for_panel <- function(df) {
  out <- lapply(split(df, list(df$Arrangement, df$Region), drop = TRUE), function(d) {
    x <- d$Value[is.finite(d$Value)]
    qs <- as.numeric(quantile(x, probs = c(0.25, 0.50, 0.75), na.rm = TRUE))
    iqr <- qs[3] - qs[1]
    lower_fence <- qs[1] - 1.5 * iqr
    upper_fence <- qs[3] + 1.5 * iqr
    inlier <- x[x >= lower_fence & x <= upper_fence]
    data.frame(
      Arrangement = d$Arrangement[1],
      Region = d$Region[1],
      ymin = min(inlier, na.rm = TRUE),
      lower = qs[1],
      middle = qs[2],
      upper = qs[3],
      ymax = max(inlier, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  box_df <- do.call(rbind, out)
  box_df$Arrangement <- factor(box_df$Arrangement, levels = names(arrangement_cols))
  box_df$Region <- factor(box_df$Region, levels = c("Mantle", "Foot"))
  box_df$x_center <- as.numeric(box_df$Arrangement) +
    ifelse(box_df$Region == "Mantle", -0.18, 0.18)
  box_df$xmin <- box_df$x_center - 0.15
  box_df$xmax <- box_df$x_center + 0.15
  box_df
}

hatch_segments_for_boxes <- function(box_df) {
  hatched <- box_df[box_df$Region == "Mantle", ]
  rows <- list()
  n <- 1
  for (i in seq_len(nrow(hatched))) {
    row <- hatched[i, ]
    box_w <- row$xmax - row$xmin
    box_h <- row$upper - row$lower
    if (!is.finite(box_h) || box_h <= 0) {
      next
    }
    starts <- seq(row$xmin - box_w, row$xmax, length.out = 7)
    for (x0 in starts) {
      x1 <- x0 + box_w
      x_start <- x0
      y_start <- row$lower
      x_end <- x1
      y_end <- row$upper
      if (x_start < row$xmin) {
        y_start <- row$lower + ((row$xmin - x_start) / box_w) * box_h
        x_start <- row$xmin
      }
      if (x_end > row$xmax) {
        y_end <- row$upper - ((x_end - row$xmax) / box_w) * box_h
        x_end <- row$xmax
      }
      if (x_start <= row$xmax && x_end >= row$xmin && y_start <= row$upper && y_end >= row$lower) {
        rows[[n]] <- data.frame(
          x = x_start, xend = x_end,
          y = y_start, yend = y_end,
          stringsAsFactors = FALSE
        )
        n <- n + 1
      }
    }
  }
  if (length(rows) == 0) {
    return(data.frame(x = numeric(), xend = numeric(), y = numeric(), yend = numeric()))
  }
  do.call(rbind, rows)
}

outliers_for_panel <- function(df, box_df) {
  rows <- list()
  n <- 1
  for (i in seq_len(nrow(box_df))) {
    row <- box_df[i, ]
    d <- df[df$Arrangement == row$Arrangement & df$Region == row$Region, ]
    x <- d$Value[is.finite(d$Value)]
    out <- x[x < row$ymin | x > row$ymax]
    if (length(out) > 0) {
      rows[[n]] <- data.frame(
        x = row$x_center,
        y = out,
        Arrangement = row$Arrangement,
        stringsAsFactors = FALSE
      )
      n <- n + 1
    }
  }
  if (length(rows) == 0) {
    return(data.frame(x = numeric(), y = numeric(), Arrangement = character()))
  }
  do.call(rbind, rows)
}
axis_breaks_for_panel <- function(ymin, ymax, tick_by = NULL) {
  if (!is.null(tick_by) && is.finite(tick_by) && tick_by > 0) {
    start <- floor(ymin / tick_by) * tick_by
    end <- ceiling(ymax / tick_by) * tick_by
    return(seq(start, end, by = tick_by))
  }
  pretty(c(ymin, ymax), n = 7)
}

make_panel <- function(spec, show_x = FALSE) {
  df <- long_for_spec(spec)
  df$Arrangement <- factor(df$Arrangement, levels = names(arrangement_cols))
  df$Region <- factor(df$Region, levels = c("Mantle", "Foot"))
  df$Variable <- spec$key

  box_df <- box_stats_for_panel(df)
  hatch_df <- hatch_segments_for_boxes(box_df)
  outlier_df <- outliers_for_panel(df, box_df)
  letters_df <- dunn_letters[dunn_letters$Variable == spec$key, ]
  letters_df <- merge(
    letters_df,
    box_df[, c("Arrangement", "Region", "ymax", "x_center")],
    by = c("Arrangement", "Region"),
    all.x = TRUE
  )
  panel_range <- diff(range(df$Value, na.rm = TRUE))
  if (panel_range == 0) {
    panel_range <- 1
  }
  letters_df$Y <- letters_df$ymax + panel_range * 0.08

  ymin <- min(df$Value, na.rm = TRUE)
  ymax <- max(c(df$Value, letters_df$Y), na.rm = TRUE)
  if (!is.null(spec$ymin)) {
    ymin <- spec$ymin
  } else {
    ymin <- ymin - panel_range * 0.10
  }
  ymax <- ymax + panel_range * 0.12

  y_breaks <- axis_breaks_for_panel(ymin, ymax, spec$tick_by)
  y_scale <- scale_y_continuous(breaks = y_breaks)
  if (!is.null(spec$smi_axis) && isTRUE(spec$smi_axis)) {
    y_scale <- scale_y_continuous(
      breaks = y_breaks,
      labels = function(x) sprintf("%.1f", x)
    )
  }

  ggplot() +
    geom_segment(
      data = box_df,
      aes(x = x_center, xend = x_center, y = ymin, yend = lower),
      linewidth = 0.40,
      color = "black"
    ) +
    geom_segment(
      data = box_df,
      aes(x = x_center, xend = x_center, y = upper, yend = ymax),
      linewidth = 0.40,
      color = "black"
    ) +
    geom_segment(
      data = box_df,
      aes(x = x_center - 0.09, xend = x_center + 0.09, y = ymin, yend = ymin),
      linewidth = 0.40,
      color = "black"
    ) +
    geom_segment(
      data = box_df,
      aes(x = x_center - 0.09, xend = x_center + 0.09, y = ymax, yend = ymax),
      linewidth = 0.40,
      color = "black"
    ) +
    geom_rect(
      data = box_df,
      aes(xmin = xmin, xmax = xmax, ymin = lower, ymax = upper, fill = Arrangement),
      color = "black",
      linewidth = 0.45
    ) +
    geom_segment(
      data = hatch_df,
      aes(x = x, xend = xend, y = y, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.33,
      color = "gray25"
    ) +
    geom_segment(
      data = box_df,
      aes(x = xmin, xend = xmax, y = middle, yend = middle),
      linewidth = 0.55,
      color = "black"
    ) +
    geom_point(
      data = outlier_df,
      aes(x = x, y = y, fill = Arrangement),
      shape = 21,
      size = 1.2,
      stroke = 0.30,
      color = "black"
    ) +
    geom_text(
      data = letters_df,
      aes(x = x_center, y = Y, label = Letter),
      inherit.aes = FALSE,
      size = 3.25,
      family = "sans",
      vjust = 0
    ) +
    scale_x_continuous(
      breaks = seq_along(arrangement_cols),
      labels = x_axis_labels,
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    scale_fill_manual(values = arrangement_cols, guide = "none") +
    y_scale +
    coord_cartesian(ylim = c(ymin, ymax), clip = "off") +
    labs(
      x = NULL,
      y = spec$axis_label,
      title = spec$panel
    ) +
    theme_classic(base_family = "sans", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      axis.title.y = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.text.x = if (show_x) {
        element_text(size = 9, angle = 25, hjust = 1)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) element_line(linewidth = 0.35) else element_blank(),
      legend.position = "none",
      plot.margin = margin(8, 8, 8, 8)
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
align_plot_grobs <- function(plot_list) {
  grobs <- lapply(plot_list, ggplotGrob)
  max_widths <- do.call(unit.pmax, lapply(grobs, function(g) g$widths))
  max_heights <- do.call(unit.pmax, lapply(grobs, function(g) g$heights))
  lapply(grobs, function(g) {
    g$widths <- max_widths
    g$heights <- max_heights
    g
  })
}

png_file <- file.path(outdir, "Fig3_ggplot2_R002_R003_D_20260807.png")
pdf_file <- file.path(outdir, "Fig3_ggplot2_R002_R003_D_20260807.pdf")


draw_hatched_key <- function(x, y, fill, label) {
  w <- 0.025
  h <- 0.18
  grid.rect(
    x = x, y = y,
    width = w, height = h,
    gp = gpar(fill = fill, col = "black", lwd = 0.5)
  )
  pushViewport(viewport(x = x, y = y, width = w, height = h, clip = "on"))
  for (offset in seq(-0.55, 0.55, length.out = 5)) {
    grid.segments(
      x0 = unit(offset, "npc"),
      y0 = unit(0, "npc"),
      x1 = unit(offset + 1, "npc"),
      y1 = unit(1, "npc"),
      gp = gpar(col = "gray25", lwd = 0.6)
    )
  }
  popViewport()
  grid.rect(
    x = x, y = y,
    width = w, height = h,
    gp = gpar(fill = NA, col = "black", lwd = 0.5)
  )
  grid.text(
    label,
    x = x + 0.022, y = y,
    just = c("left", "center"),
    gp = gpar(fontfamily = "sans", fontsize = 8.5)
  )
}

draw_plain_key <- function(x, y, fill, label) {
  w <- 0.025
  h <- 0.18
  grid.rect(
    x = x, y = y,
    width = w, height = h,
    gp = gpar(fill = fill, col = "black", lwd = 0.5)
  )
  grid.text(
    label,
    x = x + 0.022, y = y,
    just = c("left", "center"),
    gp = gpar(fontfamily = "sans", fontsize = 8.5)
  )
}

draw_manual_legend <- function() {
  pushViewport(viewport(xscale = c(0, 1), yscale = c(0, 1)))
  grid.text(
    "Spicule arrangement",
    x = 0.035, y = 0.56,
    just = c("left", "center"),
    gp = gpar(fontfamily = "sans", fontsize = 9.2, fontface = "bold")
  )
  x0 <- 0.18
  for (i in seq_along(arrangement_cols)) {
    xx <- x0 + (i - 1) * 0.13
    draw_plain_key(xx, 0.56, arrangement_cols[[i]], x_axis_labels[[i]])
  }
  grid.text(
    "Region",
    x = 0.73, y = 0.56,
    just = c("left", "center"),
    gp = gpar(fontfamily = "sans", fontsize = 9.2, fontface = "bold")
  )
  draw_hatched_key(0.805, 0.56, "gray92", "Mantle")
  draw_plain_key(0.900, 0.56, "gray92", "Foot")
  popViewport()
}
draw_figure <- function() {
  grid.newpage()
  aligned_grobs <- align_plot_grobs(plots)
  pushViewport(viewport(layout = grid.layout(
    3, 3,
    heights = unit(c(0.10, 0.45, 0.45), "npc")
  )))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  draw_manual_legend()
  popViewport()
  for (i in seq_along(aligned_grobs)) {
    row <- ceiling(i / 3) + 1
    col <- ((i - 1) %% 3) + 1
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    grid.draw(aligned_grobs[[i]])
    popViewport()
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
