# ============================================================
# utils2.R — HD Gene Expression Explorer
# Pure helper functions — no Shiny code here
# ============================================================

# ── Color palette (pink/purple theme) ───────────────────────
COL_PRIMARY   <- "#D4537E"   # deep rose — HD / Up
COL_SECONDARY <- "#7F77DD"   # soft purple — Control / Down
COL_ACCENT    <- "#ED93B1"   # light pink — highlights
COL_NS        <- "#C4C2BA"   # muted grey — not significant / filtered
COL_UP        <- "#D4537E"   # rose — upregulated
COL_DN        <- "#7F77DD"   # purple — downregulated
COL_HD        <- "#D4537E"   # rose — HD samples
COL_CTRL      <- "#7F77DD"   # purple — control samples
COL_PASSING   <- "#D4537E"   # rose — passing filter
COL_POS       <- "#D4537E"   # positive NES
COL_NEG       <- "#7F77DD"   # negative NES

# ── Null coalescing operator ─────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && a != "") a else b

# ── Shared clean theme ───────────────────────────────────────
# Applied to all ggplot outputs: white background, no grid, bold title
theme_hd <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      # White / clean background
      panel.background  = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background   = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major  = ggplot2::element_blank(),
      panel.grid.minor  = ggplot2::element_blank(),
      panel.border      = ggplot2::element_rect(colour = "#DDDDDD", fill = NA, linewidth = 0.6),
      # Bold, larger plot title
      plot.title        = ggplot2::element_text(face = "bold", size = base_size + 2,
                                                colour = "#3A1524", hjust = 0),
      plot.subtitle     = ggplot2::element_text(size = base_size - 1, colour = "#7A5060"),
      # Axis
      axis.title        = ggplot2::element_text(size = base_size, colour = "#4A2535"),
      axis.text         = ggplot2::element_text(size = base_size - 1, colour = "#555555"),
      axis.ticks        = ggplot2::element_line(colour = "#CCCCCC"),
      # Legend
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key        = ggplot2::element_rect(fill = "white", colour = NA),
      legend.title      = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text       = ggplot2::element_text(size = base_size - 1),
      legend.position   = "bottom",
      # Strip (facets, if any)
      strip.background  = ggplot2::element_rect(fill = "#FFF0F5", colour = "#F4C0D1"),
      strip.text        = ggplot2::element_text(face = "bold", colour = "#72243E")
    )
}

# ============================================================
# TAB 1 — Sample metadata
# ============================================================

#' Summary table: one row per column, type + mean/sd or distinct values
summarize_metadata <- function(df) {
  data.frame(
    Column  = names(df),
    Type    = sapply(df, function(x) class(x)[1]),
    Summary = sapply(df, function(col) {
      if (is.numeric(col)) {
        sprintf("%.2f (\u00b1 %.2f)", mean(col, na.rm = TRUE), sd(col, na.rm = TRUE))
      } else {
        vals <- unique(na.omit(as.character(col)))
        if (length(vals) <= 8) paste(vals, collapse = ", ")
        else paste0(paste(vals[1:5], collapse = ", "), " ... (", length(vals), " levels)")
      }
    }),
    stringsAsFactors = FALSE
  )
}

#' Histogram (no group) or violin + jitter (grouped)
plot_metadata_distribution <- function(df, col, grp = "None") {
  
  if (is.null(col) || col == "" || !col %in% names(df)) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = "Select a column to plot") +
             theme_hd())
  }
  if (!is.numeric(df[[col]])) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = paste(col, "is not numeric")) +
             theme_hd())
  }
  
  use_group <- !is.null(grp) && grp != "None" &&
    grp != "" && grp %in% names(df)
  
  if (use_group) {
    df[[grp]] <- as.factor(df[[grp]])
    ggplot2::ggplot(df, ggplot2::aes(x     = .data[[grp]],
                                     y     = .data[[col]],
                                     fill  = .data[[grp]],
                                     color = .data[[grp]])) +
      ggplot2::geom_violin(alpha = 0.45, trim = FALSE) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.6, size = 1.8) +
      ggplot2::scale_fill_manual(values  = c(COL_PRIMARY, COL_SECONDARY,
                                             COL_ACCENT, "#9FE1CB", "#FAC775", "#85B7EB")) +
      ggplot2::scale_color_manual(values = c(COL_PRIMARY, COL_SECONDARY,
                                             COL_ACCENT, "#9FE1CB", "#FAC775", "#85B7EB")) +
      ggplot2::labs(title = paste(col, "by", grp), x = grp, y = col) +
      theme_hd() +
      ggplot2::theme(legend.position = "none")
  } else {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[col]])) +
      ggplot2::geom_histogram(bins = 25, fill = COL_PRIMARY,
                              color = "#FBEAF0", alpha = 0.85) +
      ggplot2::labs(title = paste("Distribution of", col), x = col, y = "Count") +
      theme_hd()
  }
}

# ============================================================
# TAB 2 — Counts explorer
# ============================================================

#' Per-gene stats: median, variance, n_zero
compute_gene_stats <- function(counts_mat) {
  data.frame(
    gene     = rownames(counts_mat),
    median   = apply(counts_mat, 1, median),
    variance = apply(counts_mat, 1, var),
    n_zero   = apply(counts_mat, 1, function(x) sum(x == 0)),
    stringsAsFactors = FALSE
  )
}

#' Filter genes by variance percentile and min non-zero samples
filter_counts <- function(counts_mat, var_pct, nonzero_min) {
  stats      <- compute_gene_stats(counts_mat)
  var_thresh <- quantile(stats$variance, var_pct / 100)
  n_samples  <- ncol(counts_mat)
  passing    <- stats$gene[
    stats$variance >= var_thresh &
      (n_samples - stats$n_zero) >= nonzero_min
  ]
  counts_mat[passing, , drop = FALSE]
}

#' Scatter: median vs variance, passing genes highlighted
plot_scatter_variance <- function(counts_mat, filtered_mat) {
  stats        <- compute_gene_stats(counts_mat)
  stats$status <- ifelse(stats$gene %in% rownames(filtered_mat), "Passing", "Filtered")
  
  ggplot2::ggplot(stats, ggplot2::aes(x = log10(median + 1),
                                      y = log10(variance + 1),
                                      color = status)) +
    ggplot2::geom_point(alpha = 0.35, size = 0.7) +
    ggplot2::scale_color_manual(values = c("Passing" = COL_PASSING, "Filtered" = COL_NS)) +
    ggplot2::labs(title = "Median count vs variance",
                  x = "log10(median + 1)", y = "log10(variance + 1)", color = NULL) +
    theme_hd()
}

#' Scatter: median vs zero count, passing genes highlighted
plot_scatter_zeros <- function(counts_mat, filtered_mat) {
  stats        <- compute_gene_stats(counts_mat)
  stats$status <- ifelse(stats$gene %in% rownames(filtered_mat), "Passing", "Filtered")
  
  ggplot2::ggplot(stats, ggplot2::aes(x = log10(median + 1),
                                      y = n_zero,
                                      color = status)) +
    ggplot2::geom_point(alpha = 0.35, size = 0.7) +
    ggplot2::scale_color_manual(values = c("Passing" = COL_PASSING, "Filtered" = COL_NS)) +
    ggplot2::labs(title = "Median count vs zero count",
                  x = "log10(median + 1)", y = "Samples with zero counts", color = NULL) +
    theme_hd()
}

#' Clustered heatmap of top N genes by variance
plot_heatmap <- function(filtered_mat, top_n = 500, log_transform = TRUE) {
  n         <- min(top_n, nrow(filtered_mat))
  top_genes <- names(sort(apply(filtered_mat, 1, var), decreasing = TRUE))[1:n]
  mat       <- as.matrix(filtered_mat[top_genes, ])
  if (log_transform) mat <- log2(mat + 1)
  mat <- t(scale(t(mat)))
  mat[is.nan(mat)] <- 0
  
  pheatmap::pheatmap(
    mat,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    cluster_rows      = TRUE,
    cluster_cols      = TRUE,
    clustering_method = "ward.D2",
    color             = colorRampPalette(c("#7F77DD", "#FFF8FA", "#D4537E"))(100),
    main              = paste0("Top ", n, " genes by variance (z-scored)"),
    fontsize          = 10,
    border_color      = NA
  )
}

#' Run PCA on filtered counts
run_pca <- function(filtered_mat) {
  prcomp(t(as.matrix(filtered_mat)), scale. = TRUE)
}

#' PCA scatter plot, optionally colored by diagnosis
#' meta: data frame with columns sample_id and diagnosis (from the same
#'   metadata.csv the user uploads in Tab 1 — no separate file needed)
plot_pca <- function(pca_result, pc_x = "PC1", pc_y = "PC2", meta = NULL) {
  pct     <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 1)
  ix      <- as.integer(sub("PC", "", pc_x))
  iy      <- as.integer(sub("PC", "", pc_y))
  plot_df <- as.data.frame(pca_result$x) %>%
    tibble::rownames_to_column("sample_id")
  
  if (!is.null(meta) &&
      "sample_id" %in% names(meta) &&
      "diagnosis" %in% names(meta)) {
    plot_df <- dplyr::left_join(
      plot_df, meta[, c("sample_id", "diagnosis")], by = "sample_id"
    )
  }
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[pc_x]], y = .data[[pc_y]]))
  
  if ("diagnosis" %in% names(plot_df)) {
    # Normalise case: accept both capitalisation variants in the CSV
    plot_df$diagnosis <- dplyr::case_when(
      tolower(plot_df$diagnosis) %in% c("neurologically normal",
                                        "neurologically_normal",
                                        "control", "ctrl") ~ "Neurologically normal",
      tolower(plot_df$diagnosis) %in% c("huntington's disease",
                                        "huntington disease",
                                        "hd") ~ "Huntington's Disease",
      TRUE ~ plot_df$diagnosis
    )
    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = .data[[pc_x]], y = .data[[pc_y]])) +
      ggplot2::geom_point(ggplot2::aes(color = diagnosis), size = 2.8, alpha = 0.85) +
      ggplot2::scale_color_manual(
        values = c("Huntington's Disease"  = COL_HD,
                   "Neurologically normal" = COL_CTRL),
        na.value = COL_NS
      )
  } else {
    p <- p + ggplot2::geom_point(color = COL_PRIMARY, size = 2.8, alpha = 0.85)
  }
  
  p + ggplot2::labs(
    title = paste("PCA \u2014", pc_x, "vs", pc_y),
    x     = sprintf("%s (%.1f%%)", pc_x, pct[ix]),
    y     = sprintf("%s (%.1f%%)", pc_y, pct[iy]),
    color = "Diagnosis"
  ) + theme_hd()
}

# ============================================================
# TAB 3 — Differential expression
# ============================================================

#' Label DE genes as Up / Down / NS
label_de_direction <- function(de_df, padj_thresh, lfc_thresh) {
  de_df %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    dplyr::mutate(direction = dplyr::case_when(
      padj < padj_thresh & log2FoldChange >  lfc_thresh ~ "Up",
      padj < padj_thresh & log2FoldChange < -lfc_thresh ~ "Down",
      TRUE ~ "NS"
    ))
}

#' Volcano plot with top N labelled genes
plot_volcano <- function(de_df, padj_thresh, lfc_thresh, n_labels = 10) {
  top_genes <- de_df %>%
    dplyr::filter(direction != "NS") %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = n_labels)
  
  ggplot2::ggplot(de_df, ggplot2::aes(x = log2FoldChange,
                                      y = -log10(padj),
                                      color = direction)) +
    ggplot2::geom_point(alpha = 0.45, size = 0.9) +
    ggplot2::geom_point(data = top_genes, size = 2.2) +
    ggrepel::geom_text_repel(
      data          = top_genes,
      ggplot2::aes(label = gene),
      size          = 3,
      max.overlaps  = 15,
      color         = "#4B1528",
      segment.color = "#ED93B1"
    ) +
    ggplot2::scale_color_manual(
      values = c("Up" = COL_UP, "Down" = COL_DN, "NS" = COL_NS)
    ) +
    ggplot2::geom_vline(xintercept = c(-lfc_thresh, lfc_thresh),
                        linetype = "dashed", color = "#ED93B1", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(padj_thresh),
                        linetype = "dashed", color = "#ED93B1", linewidth = 0.5) +
    ggplot2::labs(title = "Volcano plot \u2014 HD vs Control",
                  x = "log2 fold change", y = "-log10(padj)", color = NULL) +
    theme_hd()
}

# ============================================================
# TAB 4 — GSEA
# ============================================================

#' Strip HALLMARK_ prefix and underscores for readable labels
shorten_pathway_names <- function(pathway_vec) {
  pathway_vec %>%
    stringr::str_replace("HALLMARK_", "") %>%
    stringr::str_replace_all("_", " ") %>%
    tolower()
}

#' Top N pathways sorted by padj, prepped for barplot
prep_gsea_barplot_df <- function(gsea_df, top_n) {
  gsea_df %>%
    dplyr::filter(!is.na(padj), !is.na(NES)) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::arrange(NES) %>%
    dplyr::mutate(
      pathway_short = shorten_pathway_names(pathway),
      pathway_short = factor(pathway_short, levels = pathway_short),
      direction     = ifelse(NES > 0, "Positive", "Negative")
    )
}

#' NES barplot
plot_gsea_barplot <- function(gsea_top_df, top_n) {
  ggplot2::ggplot(gsea_top_df,
                  ggplot2::aes(x = NES, y = pathway_short, fill = direction)) +
    ggplot2::geom_col(alpha = 0.88) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "#993556") +
    ggplot2::scale_fill_manual(values = c("Positive" = COL_POS, "Negative" = COL_NEG)) +
    ggplot2::labs(
      title = paste("Top", top_n, "pathways by adjusted p-value"),
      x     = "Normalized enrichment score (NES)",
      y     = NULL,
      fill  = "NES direction"
    ) +
    theme_hd() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9, color = "#72243E")
    )
}

#' Filter GSEA table by padj and NES direction
filter_gsea_table <- function(gsea_df, padj_max, nes_filter = "All") {
  df <- gsea_df %>%
    dplyr::filter(!is.na(padj), padj <= padj_max)
  if (nes_filter == "Positive") df <- dplyr::filter(df, NES > 0)
  if (nes_filter == "Negative") df <- dplyr::filter(df, NES < 0)
  df %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))
}

#' NES vs -log10(padj) scatter
plot_gsea_scatter <- function(gsea_df, padj_thresh) {
  df <- gsea_df %>%
    dplyr::filter(!is.na(padj), !is.na(NES)) %>%
    dplyr::mutate(
      log10p    = -log10(padj),
      highlight = padj <= padj_thresh,
      label     = ifelse(highlight, shorten_pathway_names(pathway), NA)
    )
  
  ggplot2::ggplot(df, ggplot2::aes(x = NES, y = log10p, color = highlight)) +
    ggplot2::geom_point(alpha = 0.75, size = 2) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label    = label),
      size         = 2.8,
      na.rm        = TRUE,
      max.overlaps = 20,
      color        = "#4B1528",
      segment.color = "#ED93B1"
    ) +
    ggplot2::geom_hline(yintercept = -log10(padj_thresh),
                        linetype = "dashed", color = "#ED93B1", linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = 0,
                        linetype = "solid", color = "#F4C0D1", linewidth = 0.4) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = COL_NS, "TRUE" = COL_PRIMARY),
      labels = c("Not significant", paste0("padj \u2264 ", padj_thresh))
    ) +
    ggplot2::labs(
      title = "GSEA \u2014 NES vs significance",
      x     = "Normalized enrichment score (NES)",
      y     = "-log10(adjusted p-value)",
      color = NULL
    ) +
    theme_hd()
}