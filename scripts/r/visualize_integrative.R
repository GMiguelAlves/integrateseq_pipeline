#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The R package ggplot2 is required for visualization.")
  }
})

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = "") {
  idx <- match(flag, args)
  if (is.na(idx) || idx >= length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

project_dir <- normalizePath(get_arg("--project-dir", "."), mustWork = FALSE)
outdir <- normalizePath(get_arg("--outdir", file.path(project_dir, "090-visualizations")), mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

manifest <- data.frame(figure = character(), status = character(), stringsAsFactors = FALSE)

add_manifest <- function(figure, status = "created") {
  manifest <<- rbind(manifest, data.frame(figure = figure, status = status, stringsAsFactors = FALSE))
}

read_tsv <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(stringsAsFactors = FALSE))
  }
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  read.delim(
    con,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

as_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_col <- function(df, name, default = "") {
  if (name %in% names(df)) {
    return(df[[name]])
  }
  rep(default, nrow(df))
}

trim_label <- function(x, n = 42) {
  x <- as.character(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 3), "..."), x)
}

save_plot <- function(plot, stem, width = 8, height = 5) {
  outputs <- c("png", "pdf", "svg")
  for (ext in outputs) {
    path <- file.path(outdir, paste0(stem, ".", ext))
    ok <- tryCatch({
      if (ext == "png") {
        png(path, width = width, height = height, units = "in", res = 160)
      } else if (ext == "pdf") {
        pdf(path, width = width, height = height)
      } else {
        svg(path, width = width, height = height)
      }
      print(plot)
      dev.off()
      TRUE
    }, error = function(e) {
      if (dev.cur() > 1) {
        dev.off()
      }
      warning(sprintf("Could not write %s: %s", basename(path), e$message))
      FALSE
    })
    add_manifest(basename(path), if (ok) "created" else "error")
  }
}

empty_plot <- function(title, message) {
  ggplot() +
    annotate("text", x = 0, y = 0.1, label = title, fontface = "bold", size = 5) +
    annotate("text", x = 0, y = -0.08, label = message, size = 3.5) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void()
}

theme_integrative <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", color = "#111827"),
      plot.subtitle = element_text(color = "#4b5563"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(color = "#374151"),
      axis.text = element_text(color = "#374151"),
      strip.text = element_text(face = "bold", color = "#111827"),
      legend.position = "bottom"
    )
}

palette_marks <- c(
  "#2563eb", "#16a34a", "#f59e0b", "#dc2626", "#7c3aed",
  "#0f766e", "#db2777", "#64748b", "#0891b2", "#ea580c"
)

class_counts <- read_tsv(file.path(project_dir, "070-integrated-tables", "integrative_class_counts.tsv"))
candidate_scores <- read_tsv(file.path(project_dir, "080-candidate-scoring", "candidate_gene_scores.tsv"))
catalog <- read_tsv(file.path(project_dir, "030-id-harmonization", "epigenetic_machinery_catalog.tsv"))
chip_mark_stage <- read_tsv(file.path(project_dir, "060-chipseq-summary", "chip_mark_stage_metadata.tsv"))
gene_mark_summary <- read_tsv(file.path(project_dir, "070-integrated-tables", "gene_mark_stage_summary.tsv"))
gene_mark_links <- read_tsv(file.path(project_dir, "070-integrated-tables", "gene_mark_stage_links.tsv"))
stage_mark_comparison <- read_tsv(file.path(project_dir, "080-candidate-scoring", "stage_mark_comparison.tsv"))

if (nrow(class_counts) > 0 && all(c("integrative_class", "n_genes") %in% names(class_counts))) {
  class_counts$n_genes <- as_num(class_counts$n_genes)
  p <- ggplot(class_counts, aes(x = reorder(integrative_class, n_genes), y = n_genes, fill = integrative_class)) +
    geom_col(width = 0.72, show.legend = FALSE) +
    coord_flip() +
    scale_fill_manual(values = rep(palette_marks, length.out = nrow(class_counts))) +
    labs(title = "Integrative classes", x = NULL, y = "Genes") +
    theme_integrative()
} else {
  p <- empty_plot("Integrative classes", "No integrated class table was found.")
}
save_plot(p, "barplot_integrative_classes", 8, 4.8)

if (nrow(candidate_scores) > 0 && "candidate_score" %in% names(candidate_scores)) {
  candidate_scores$candidate_score <- as_num(candidate_scores$candidate_score)
  candidate_scores <- candidate_scores[order(-candidate_scores$candidate_score), , drop = FALSE]
  top <- head(candidate_scores, 30)
  top$gene_label <- trim_label(ifelse(safe_col(top, "gene_name") != "", safe_col(top, "gene_name"), safe_col(top, "gene_id")))
  p <- ggplot(top, aes(x = reorder(gene_label, candidate_score), y = candidate_score)) +
    geom_col(fill = "#16a34a", width = 0.72) +
    coord_flip() +
    labs(title = "Top candidate genes", x = NULL, y = "Candidate score") +
    theme_integrative()
} else {
  p <- empty_plot("Top candidate genes", "No candidate scoring table was found.")
}
save_plot(p, "top_candidate_scores", 8.5, 6)

if (nrow(catalog) > 0) {
  group <- safe_col(catalog, "machinery_group", "unknown")
  groups <- as.data.frame(sort(table(group), decreasing = TRUE), stringsAsFactors = FALSE)
  names(groups) <- c("machinery_group", "n_genes")
  p <- ggplot(groups, aes(x = reorder(machinery_group, n_genes), y = n_genes)) +
    geom_col(fill = "#0f766e", width = 0.72) +
    coord_flip() +
    labs(title = "Epigenetic machinery catalog", x = NULL, y = "Genes") +
    theme_integrative()
} else {
  p <- empty_plot("Epigenetic machinery catalog", "No epigenetic machinery catalog was found.")
}
save_plot(p, "epigenetic_catalog_groups", 8, 5)

plot_heatmap <- function(df, row_col, col_col, value_col, title, fill_label, empty_message) {
  if (nrow(df) == 0 || !(row_col %in% names(df)) || !(col_col %in% names(df))) {
    return(empty_plot(title, empty_message))
  }
  if (value_col %in% names(df)) {
    value <- as_num(df[[value_col]])
  } else {
    value <- rep(1, nrow(df))
  }
  plot_df <- data.frame(row = df[[row_col]], col = df[[col_col]], value = value, stringsAsFactors = FALSE)
  plot_df$row[plot_df$row == ""] <- "unknown"
  plot_df$col[plot_df$col == ""] <- "unknown"
  plot_df <- aggregate(value ~ row + col, plot_df, sum, na.rm = TRUE)
  ggplot(plot_df, aes(x = col, y = row, fill = value)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = value), size = 3, color = "#111827") +
    scale_fill_gradient(low = "#e0f2fe", high = "#075985") +
    labs(title = title, x = "Stage or condition", y = "Mark or factor", fill = fill_label) +
    theme_integrative() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}

p <- plot_heatmap(
  chip_mark_stage,
  "mark_or_factor",
  "stage_or_condition",
  "n_samples",
  "ChIP-seq marks by life-cycle stage",
  "Samples",
  "No ChIP mark-stage metadata was found."
)
save_plot(p, "chip_mark_stage_matrix", 8.5, 6)

p <- plot_heatmap(
  gene_mark_summary,
  "mark_or_factor",
  "stage_or_condition",
  "n_peaks",
  "Gene-mark-stage links",
  "Peaks",
  "No gene-mark-stage links were found."
)
save_plot(p, "gene_mark_stage_matrix", 8.5, 6)

if (nrow(stage_mark_comparison) > 0 && all(c("mark_or_factor", "stage_or_condition") %in% names(stage_mark_comparison))) {
  stage_mark_comparison$n_deg_linked_genes <- as_num(safe_col(stage_mark_comparison, "n_deg_linked_genes", 0))
  stage_mark_comparison$n_epigenetic_machinery_genes <- as_num(safe_col(stage_mark_comparison, "n_epigenetic_machinery_genes", 0))
  stage_mark_comparison$n_wgcna_hits <- as_num(safe_col(stage_mark_comparison, "n_wgcna_hits", 0))
  stage_mark_comparison$n_mfuzz_hits <- as_num(safe_col(stage_mark_comparison, "n_mfuzz_hits", 0))
  stage_mark_comparison$n_dtu_hits <- as_num(safe_col(stage_mark_comparison, "n_dtu_hits", 0))
  stage_mark_comparison$n_splicing_hits <- as_num(safe_col(stage_mark_comparison, "n_splicing_hits", 0))
  stage_mark_comparison$integrated_evidence <- with(
    stage_mark_comparison,
    n_deg_linked_genes + n_epigenetic_machinery_genes + n_wgcna_hits + n_mfuzz_hits + n_dtu_hits + n_splicing_hits
  )
  p <- ggplot(stage_mark_comparison, aes(x = stage_or_condition, y = mark_or_factor, fill = integrated_evidence)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = integrated_evidence), size = 3, color = "#111827") +
    scale_fill_gradient(low = "#ecfeff", high = "#0f766e") +
    labs(
      title = "Integrated evidence by stage and mark",
      subtitle = "Counts combine DE-linked genes, epigenetic machinery, WGCNA, Mfuzz, DTU, and splicing evidence.",
      x = "Stage or condition",
      y = "Mark or factor",
      fill = "Evidence"
    ) +
    theme_integrative() +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
} else {
  p <- empty_plot("Integrated evidence by stage and mark", "Run scoring to populate stage-mark comparisons.")
}
save_plot(p, "stage_mark_integrated_evidence", 8.5, 6)

if (nrow(gene_mark_links) > 0) {
  links <- gene_mark_links
  links$mark_or_factor <- safe_col(links, "mark_or_factor", "unknown")
  links$stage_or_condition <- safe_col(links, "stage_or_condition", "unknown")
  links$promoter_flag <- safe_col(links, "promoter_flag", "false")
  links$peak_chrom <- safe_col(links, "peak_chrom", "")
  links$peak_start <- as_num(safe_col(links, "peak_start", NA))
  links$peak_end <- as_num(safe_col(links, "peak_end", NA))
  links$peak_midpoint <- as_num(safe_col(links, "peak_midpoint", NA))
  links$distance_to_tss <- as_num(safe_col(links, "distance_to_tss", NA))
  links$gene_label <- ifelse(
    safe_col(links, "gene_name") != "" & safe_col(links, "gene_name") != safe_col(links, "gene_id"),
    paste0(safe_col(links, "gene_name"), " (", safe_col(links, "gene_id"), ")"),
    safe_col(links, "gene_id")
  )
  links$gene_label <- trim_label(links$gene_label, 38)
  gene_counts <- sort(table(links$gene_label), decreasing = TRUE)
  keep_genes <- names(head(gene_counts, 25))
  links <- links[links$gene_label %in% keep_genes, , drop = FALSE]
  links$gene_label <- factor(links$gene_label, levels = rev(keep_genes))
  links$peak_width <- pmax(1, links$peak_end - links$peak_start)

  has_distance <- any(!is.na(links$distance_to_tss))
  if (has_distance) {
    plot_df <- links[!is.na(links$distance_to_tss), , drop = FALSE]
    p <- ggplot(plot_df, aes(x = distance_to_tss, y = gene_label, color = mark_or_factor)) +
      geom_vline(xintercept = 0, color = "#374151", linewidth = 0.4) +
      geom_point(aes(shape = promoter_flag, size = peak_width), alpha = 0.82) +
      facet_wrap(~stage_or_condition, scales = "free_x") +
      scale_color_manual(values = rep(palette_marks, length.out = length(unique(plot_df$mark_or_factor)))) +
      scale_size_continuous(range = c(1.8, 5.5), guide = "none") +
      labs(
        title = "Gene-position-mark associations",
        subtitle = "Each point is a ChIP peak linked to a gene; x = peak distance to TSS.",
        x = "Peak distance to TSS (bp)",
        y = "Gene",
        color = "Mark",
        shape = "Promoter"
      ) +
      theme_integrative()
  } else if (any(!is.na(links$peak_start)) && any(!is.na(links$peak_end))) {
    plot_df <- links[!is.na(links$peak_start) & !is.na(links$peak_end), , drop = FALSE]
    p <- ggplot(plot_df, aes(y = gene_label, color = mark_or_factor)) +
      geom_segment(aes(x = peak_start, xend = peak_end, yend = gene_label), linewidth = 2.4, alpha = 0.8) +
      geom_point(aes(x = peak_midpoint, shape = promoter_flag), size = 2.3) +
      facet_grid(peak_chrom ~ stage_or_condition, scales = "free_x", space = "free_x") +
      scale_color_manual(values = rep(palette_marks, length.out = length(unique(plot_df$mark_or_factor)))) +
      labs(
        title = "Gene-position-mark associations",
        subtitle = "Segments show peak genomic intervals linked to each gene.",
        x = "Peak genomic coordinate (bp)",
        y = "Gene",
        color = "Mark",
        shape = "Promoter"
      ) +
      theme_integrative()
  } else {
    p <- empty_plot("Gene-position-mark associations", "The link table exists, but peak positions are missing.")
  }
} else {
  p <- empty_plot("Gene-position-mark associations", "Run integration after peak-gene mapping to populate this figure.")
}
save_plot(p, "gene_position_mark_map", 10.5, 7.2)

workflow <- data.frame(
  xmin = c(0.05, 0.05, 0.38, 0.72),
  xmax = c(0.27, 0.27, 0.60, 0.94),
  ymin = c(0.60, 0.18, 0.39, 0.39),
  ymax = c(0.82, 0.40, 0.66, 0.66),
  title = c("RNA-seq", "ChIP-seq", "Integration", "Outputs"),
  body = c("expression + DEG", "marks + peaks", "gene x position x mark", "figures + HTML report"),
  fill = c("#dbeafe", "#dcfce7", "#fef3c7", "#ede9fe"),
  stringsAsFactors = FALSE
)

p <- ggplot() +
  geom_rect(
    data = workflow,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    color = "#374151",
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  annotate("text", x = (workflow$xmin + workflow$xmax) / 2, y = workflow$ymax - 0.07, label = workflow$title, fontface = "bold", size = 4) +
  annotate("text", x = (workflow$xmin + workflow$xmax) / 2, y = workflow$ymin + 0.07, label = workflow$body, size = 3.2) +
  annotate("segment", x = 0.27, xend = 0.38, y = 0.71, yend = 0.54, arrow = grid::arrow(length = grid::unit(0.018, "npc")), color = "#374151") +
  annotate("segment", x = 0.27, xend = 0.38, y = 0.29, yend = 0.51, arrow = grid::arrow(length = grid::unit(0.018, "npc")), color = "#374151") +
  annotate("segment", x = 0.60, xend = 0.72, y = 0.525, yend = 0.525, arrow = grid::arrow(length = grid::unit(0.018, "npc")), color = "#374151") +
  labs(title = "Integrative analysis workflow") +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, color = "#111827"))
save_plot(p, "integrative_workflow_overview", 9, 3.6)

write.table(
  manifest,
  file = file.path(outdir, "visualization_manifest.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
