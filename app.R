# ============================================================
# app.R — HD Gene Expression Explorer
# BF591 Final Project | Dataset: GSE64810
# ============================================================

library(shiny)
library(shinythemes)
library(DT)
library(ggplot2)
library(dplyr)
library(tibble)
library(ggrepel)
library(stringr)

source("utils2.R")

# Raise the default 5 MB per-upload limit to 200 MB so users can
# upload large normalized counts matrices without an error.
options(shiny.maxRequestSize = 200 * 1024^2)

# ============================================================
# Shared CSS
# ============================================================
app_css <- "
/* ── Global ───────────────────────────────────────── */
body { background-color: #FAFAFA; }

/* ── Navbar ─────────────────────────────────────────*/
/* Override flatly defaults with the project rose palette */
.navbar-default { background-color: #72243E !important; border-color: #5A1B30 !important; }
.navbar-default .navbar-brand,
.navbar-default .navbar-nav > li > a { color: #FFE4EF !important; font-weight: 600; }
.navbar-default .navbar-nav > .active > a,
.navbar-default .navbar-nav > li > a:hover {
  background-color: #5A1B30 !important; color: #fff !important;
}

/* ── Sidebar ─────────────────────────────────────── */
/* Shiny wraps sidebarPanel() content in a <div class='well'> */
.well {
  background-color: #FFF5F8;
  border: 1px solid #F4C0D1;
  border-radius: 8px;
  padding: 14px;
  /* overflow-wrap + word-break prevent long labels (e.g. 'Downregulated')
     from spilling outside the sidebar panel */
  overflow-wrap: break-word;
  word-break: break-word;
}

/* ── Sidebar headings ──────────────────────────────*/
.well h5 { color: #72243E; font-weight: 700; margin-bottom: 6px; }
.well hr  { border-color: #F4C0D1; margin: 10px 0; }

/* ── Stat boxes ─────────────────────────────────── */
/* Used by renderUI() blocks to display summary counts
   (samples, total genes, passing genes, etc.) */
.stat-box {
  background: white;
  border: 1px solid #F4C0D1;
  border-radius: 8px;
  padding: 10px 6px;
  text-align: center;
  margin-bottom: 10px;
}
.stat-val { font-size: 22px; font-weight: 700; color: #D4537E; line-height: 1.2; }
.stat-lbl { font-size: 11px; color: #8A5060; margin-top: 2px; }

/* ── Upregulated / Downregulated badges ───────────*/
/* Full-width badges replace the previous two-column stat-box
   layout, which was too narrow and clipped the text */
.dir-badge {
  display: block;
  border-radius: 6px;
  padding: 8px 4px;
  text-align: center;
  margin-bottom: 6px;
  font-size: 12px;
  font-weight: 600;
  white-space: normal;   /* allow text to wrap inside badge */
  word-break: break-word;
}
.badge-up   { background: #FBEAF0; color: #D4537E; border: 1px solid #F4C0D1; }
.badge-down { background: #EEEDF9; color: #7F77DD; border: 1px solid #C5C2EE; }

/* ── Download button ────────────────────────────── */
/* width:100% + white-space:normal prevent the button label
   from overflowing the sidebar panel at narrow widths */
.btn-download {
  width: 100%;
  background-color: #D4537E;
  color: white;
  border: none;
  border-radius: 6px;
  padding: 7px 10px;
  font-size: 13px;
  margin-top: 6px;
  white-space: normal;
  word-break: break-word;
}
.btn-download:hover { background-color: #B83C65; color: white; }

/* ── Tab panels ─────────────────────────────────── */
.nav-tabs > li.active > a,
.nav-tabs > li.active > a:hover { color: #72243E; border-bottom-color: white; font-weight: 600; }
.nav-tabs > li > a:hover { color: #D4537E; }

/* ── Instruction / note text ────────────────────── */
/* .instruct : rose-coloured guidance line at top of each tab */
/* .note-text: italic helper hint inside the sidebar          */
.instruct { color: #993556; font-size: 13px; margin-bottom: 6px; }
.note-text { color: #7A5060; font-size: 12px; font-style: italic; }
"

# ============================================================
# UI
# ============================================================
ui <- navbarPage(
  title = "HD Gene Expression Explorer",
  theme = shinytheme("flatly"),                      # Bootstrap "flatly" base
  header = tags$head(tags$style(HTML(app_css))),     # inject custom CSS
  
  # ── TAB 1: Sample info ──────────────────────────────────────
  # Users upload a metadata CSV and explore it through three
  # sub-tabs: a column summary, the full data table, and plots.
  tabPanel("Sample info",
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h5("Upload data"),
               p("Upload a metadata CSV file with one row per sample.",
                 class = "note-text"),

               fileInput("sample_file", "Sample metadata (.csv)",
                         accept = ".csv", placeholder = "metadata.csv"),
               hr(),
               h5("Plot controls"),
               p("Select a numeric column to visualise its distribution.
           Optionally group by a categorical column to compare groups.",
                 class = "note-text"),

               uiOutput("hist_col_ui"),    # numeric column to plot
               uiOutput("group_col_ui")    # optional categorical group-by
             ),
             mainPanel(
               width = 9,
               tabsetPanel(
                 # Summary: row/column counts + per-column type summary table
                 tabPanel("Summary",
                          br(),
                          p("High-level counts and a per-column summary of the uploaded metadata.",
                            class = "instruct"),
                          uiOutput("sample_meta_counts"),    # stat boxes: N samples, N columns
                          br(),
                          tableOutput("sample_summary_tbl") # column | type | mean±sd or distinct values
                 ),
                 # Table: full interactive DT table with search and sort
                 tabPanel("Table",
                          br(),
                          p("Browse, search, and sort the full metadata table.",
                            class = "instruct"),
                          DT::dataTableOutput("sample_dt")
                 ),
                 # Plots: histogram (ungrouped) or violin+jitter (grouped)
                 tabPanel("Plots",
                          br(),
                          p("Distribution of the selected metadata column.
               Choose a grouping variable in the sidebar to switch to a violin plot.",
                            class = "instruct"),
                          plotOutput("sample_plot", height = "420px")
                 )
               )
             )
           )
  ),
  
  # ── TAB 2: Counts explorer ──────────────────────────────────
  # Users upload a normalized counts matrix, adjust two gene-
  # filtering sliders, and explore the filtered data through
  # four sub-tabs: filter summary, scatter plots, heatmap, PCA.
  tabPanel("Counts explorer",
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h5("Upload data"),
               # \u00d7 renders as the multiplication sign ×
               p("Upload a normalized counts matrix (genes \u00d7 samples).",
                 class = "note-text"),
               fileInput("counts_file", "Normalized counts (.csv)",
                         accept = ".csv", placeholder = "normalized_counts.csv"),
               hr(),
               h5("Gene filters"),
               p("Genes below the variance percentile or with too many zeros will be excluded.",
                 class = "note-text"),
               # Variance filter: retain genes at or above this percentile
               sliderInput("var_pct", "Min variance percentile",
                           min = 0, max = 100, value = 50, step = 1, post = "th"),
               # Zero filter: retain genes with at least this many non-zero samples
               sliderInput("nonzero_n", "Min non-zero samples",
                           min = 0, max = 69, value = 10, step = 1),
               hr(),
               h5("PCA controls"),
               p("Choose which principal components to display on each axis.",
                 class = "note-text"),
               # PC selectors: rendered dynamically after a file is loaded
               # because the number of available PCs depends on sample count
               uiOutput("pc_x_ui"),    # X-axis PC
               uiOutput("pc_y_ui"),    # Y-axis PC
               # When checked, a second file input appears so users can supply
               # metadata for coloring PCA points by diagnosis
               checkboxInput("pca_color_meta",
                             "Color samples by diagnosis", value = FALSE),
               # conditionalPanel evaluates a JavaScript expression client-side;
               # the panel is only visible when the checkbox is TRUE
               conditionalPanel(
                 "input.pca_color_meta",
                 p("Re-upload the same metadata CSV used in the Sample Info tab to
             color samples by diagnosis.", class = "note-text"),
                 fileInput("pca_meta_file", "Metadata for coloring (.csv)",
                           accept = ".csv")
               )
             ),
             mainPanel(
               width = 9,
               tabsetPanel(
                 # Filter summary: stat boxes + table showing pass/fail gene counts
                 tabPanel("Filter summary",
                          br(),
                          p("Effect of the current variance and zero-count filters on gene count.",
                            class = "instruct"),
                          uiOutput("filter_stats"),         # four stat boxes
                          br(),
                          tableOutput("filter_summary_tbl") # tabular version of the same info
                 ),
                 # Scatter plots: two side-by-side diagnostic plots
                 tabPanel("Scatter plots",
                          br(),
                          p("Each point is a gene. Highlighted genes pass the current filters.
               Left: median vs variance; right: median vs zero-count.",
                            class = "instruct"),
                          # fluidRow + column(6) places both plots side by side
                          fluidRow(
                            column(6, plotOutput("scatter_var",   height = "320px")),
                            column(6, plotOutput("scatter_zeros", height = "320px"))
                          )
                 ),
                 # Heatmap: clustered heatmap of top N variable genes after filtering
                 tabPanel("Heatmap",
                          br(),
                          p("Clustered heatmap of the top N most variable genes after filtering.
               Counts are z-scored per gene.",
                            class = "instruct"),
                          # Controls placed above the plot so users adjust before it renders
                          checkboxInput("hm_log", "Log\u2082-transform counts", value = TRUE),
                          numericInput("hm_topn", "Top N genes by variance",
                                       value = 500, min = 50, max = 5000, step = 50),
                          plotOutput("heatmap_plot", height = "500px")
                 ),
                 # PCA: scatter of two user-selected principal components
                 tabPanel("PCA",
                          br(),
                          # \u2018 / \u2019 are typographic left/right single quotes
                          p("Principal component analysis on filtered, scaled counts.
               Toggle \u2018Color samples by diagnosis\u2019 in the sidebar and
               upload the metadata CSV to colour by group.",
                            class = "instruct"),
                          plotOutput("pca_plot", height = "420px")
                 )
               )
             )
           )
  ),
  
  # ── TAB 3: Differential expression ─────────────────────────
  # Users upload DESeq2 results and interactively explore them.
  # Two threshold sliders (padj, |log2FC|) reactively update the
  # volcano plot and the up/down-regulated gene count badges.
  tabPanel("Differential expression",
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h5("Upload data"),
               p("Upload a DESeq2 results CSV with columns: gene, log2FoldChange, padj.",
                 class = "note-text"),
               fileInput("de_file", "DESeq2 results (.csv)",
                         accept = ".csv", placeholder = "deseq2_results.csv"),
               hr(),
               h5("Volcano controls"),
               p("Adjust thresholds to change which genes are called significant.",
                 class = "note-text"),
               # padj threshold: genes below this adjusted p-value are considered significant
               sliderInput("padj_thresh", "padj threshold",
                           min = 0.001, max = 0.2, value = 0.05, step = 0.001),
               # |log2FC| threshold: genes must also exceed this fold-change magnitude
               sliderInput("lfc_thresh", "|log\u2082FC| threshold",
                           min = 0, max = 3, value = 1, step = 0.1),
               hr(),
               # Reactive up/down count badges — rebuilt whenever thresholds change
               uiOutput("de_sig_counts")
             ),
             mainPanel(
               width = 9,
               tabsetPanel(
                 # Volcano plot: colored scatter of all genes with threshold lines
                 tabPanel("Volcano plot",
                          br(),
                          p("Each point is a gene. Dashed lines mark the current padj and
               log\u2082FC thresholds. Top significant genes are labelled.",
                            class = "instruct"),
                          plotOutput("volcano_plot", height = "450px")
                 ),
                 # Results table: full DESeq2 output with fold-change cell shading
                 tabPanel("Results table",
                          br(),
                          p("Full DESeq2 results table, sortable and searchable. Cells are
               shaded by fold-change direction.",
                            class = "instruct"),
                          DT::dataTableOutput("de_dt")
                 )
               )
             )
           )
  ),
  
  # ── TAB 4: GSEA ─────────────────────────────────────────────
  # Displays fgsea pathway enrichment results across three sub-tabs:
  #   Barplot — NES bar chart; click a bar to see pathway details
  #   Table   — filterable/downloadable results table
  #   Scatter — NES vs -log10(padj) with significance highlighting
  #
  # The sidebar uses conditionalPanel() so only the controls
  # relevant to the currently active sub-tab are shown.
  tabPanel("GSEA",
           sidebarLayout(
             sidebarPanel(
               width = 3,
               h5("Upload data"),
               p("Upload an fgsea results CSV with columns: pathway, NES, pval, padj, size.",
                 class = "note-text"),
               fileInput("gsea_file", "fgsea results (.csv)",
                         accept = ".csv", placeholder = "fgsea_results.csv"),
               hr(),
               
               # ── Barplot controls (visible only when Barplot tab is active) ──
               conditionalPanel("input.gsea_tabs == 'Barplot'",
                                h5("Barplot controls"),
                                p("Show the N most significant pathways ranked by adjusted p-value.",
                                  class = "note-text"),
                                sliderInput("gsea_top_n", "Top N pathways",
                                            min = 5, max = 50, value = 10, step = 1)
               ),
               
               # ── Table controls (visible only when Table tab is active) ──────
               conditionalPanel("input.gsea_tabs == 'Table'",
                                h5("Table controls"),
                                p("Filter pathways by adjusted p-value and NES direction.",
                                  class = "note-text"),
                                sliderInput("gsea_padj_tbl", "Max padj",
                                            min = 0.001, max = 1, value = 0.25, step = 0.001),
                                # Named vector maps user-friendly labels to server-side values:
                                # "Upregulated" -> "Positive", "Downregulated" -> "Negative"
                                tags$div(
                                  style = "margin-bottom: 6px;",
                                  tags$strong("NES direction", style = "font-size:13px; color:#4A2535;"),
                                  br(),
                                  radioButtons("gsea_nes_filter", label = NULL,
                                               choices  = c("All", "Upregulated" = "Positive",
                                                            "Downregulated" = "Negative"),
                                               selected = "All")
                                ),
                                br(),
                                # btn-download CSS class gives this button a full-width rose style
                                # and prevents the label from overflowing the sidebar
                                downloadButton("gsea_download", "Download filtered table",
                                               class = "btn-download")
               ),
               
               # ── Scatter controls (visible only when Scatter tab is active) ──
               conditionalPanel("input.gsea_tabs == 'Scatter'",
                                h5("Scatter controls"),
                                p("Pathways below this padj threshold are highlighted and labelled.",
                                  class = "note-text"),
                                sliderInput("gsea_padj_scatter", "Highlight padj threshold",
                                            min = 0.001, max = 1, value = 0.25, step = 0.001)
               )
             ),
             
             mainPanel(
               width = 9,
               tabsetPanel(id = "gsea_tabs",
                           
                           # Barplot sub-tab
                           tabPanel("Barplot",
                                    br(),
                                    p("Bar length = Normalized Enrichment Score (NES). Positive NES
               (rose) = pathway upregulated in HD; negative NES (purple) =
               downregulated. Click a bar to see details.",
                                      class = "instruct"),
                                    # click = "gsea_bar_click" sends click coordinates to the server
                                    # as input$gsea_bar_click; used to identify which bar was clicked
                                    plotOutput("gsea_barplot", height = "420px",
                                               click = "gsea_bar_click"),
                                    br(),
                                    uiOutput("gsea_clicked_info")  # detail card, appears after a click
                           ),
                           
                           # Table sub-tab
                           tabPanel("Table",
                                    br(),
                                    p("Filtered fgsea results table. Use the sidebar to restrict by
               padj and NES direction, then download.",
                                      class = "instruct"),
                                    DT::dataTableOutput("gsea_dt")
                           ),
                           
                           # Scatter sub-tab
                           tabPanel("Scatter",
                                    br(),
                                    p("NES vs -log10(padj). Highlighted gene sets are below the padj
               threshold set in the sidebar.",
                                      class = "instruct"),
                                    plotOutput("gsea_scatter", height = "420px")
                           )
               )
             )
           )
  )
)

# ============================================================
# SERVER
# ============================================================
server <- function(input, output, session) {
  
  # ── TAB 1: Sample metadata ────────────────────────────────────
  
  # Read the uploaded metadata CSV and cache it. 
  sample_data <- reactive({
    req(input$sample_file)   # wait until the user has selected a file
    df <- read.csv(input$sample_file$datapath, stringsAsFactors = FALSE)
    validate(need(nrow(df) > 0, "Uploaded file is empty"))
    df
  })
  
  # Dynamically build the "Column to plot" dropdown after upload.
  output$hist_col_ui <- renderUI({
    req(input$sample_file)
    df       <- req(sample_data())
    num_cols <- names(df)[sapply(df, is.numeric)]
    validate(need(length(num_cols) > 0, "No numeric columns found"))
    selectInput("hist_col", "Column to plot",
                choices  = num_cols,
                selected = num_cols[1])
  })
  
  # Dynamically build the "Group by" dropdown after upload.
  output$group_col_ui <- renderUI({
    req(input$sample_file)
    df       <- req(sample_data())
    cat_cols <- names(df)[sapply(df, function(x) is.character(x) | is.factor(x))]
    cat_cols <- cat_cols[!cat_cols %in% c("sample_id", "Sample_ID", "SampleID")]
    selectInput("group_col", "Group by (optional)",
                choices  = c("None", cat_cols),
                selected = "None")
  })
  
  # Render two stat boxes showing the total number of samples (rows)
  # and metadata columns.
  output$sample_meta_counts <- renderUI({
    req(input$sample_file)
    df <- req(sample_data())
    fluidRow(
      column(3, div(class = "stat-box",
                    div(class = "stat-val", nrow(df)),
                    div(class = "stat-lbl", "samples"))),
      column(3, div(class = "stat-box",
                    div(class = "stat-val", ncol(df)),
                    div(class = "stat-lbl", "columns")))
    )
  })
  
  # Render the per-column summary table.
  output$sample_summary_tbl <- renderTable({
    req(input$sample_file)
    summarize_metadata(req(sample_data()))
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  # Render the full metadata as an interactive DT table.
  output$sample_dt <- DT::renderDataTable({
    req(input$sample_file)
    DT::datatable(req(sample_data()),
                  options  = list(pageLength = 15, scrollX = TRUE),
                  rownames = FALSE, filter = "top")
  })
  
  # Render the distribution plot by calling plot_metadata_distribution()
  output$sample_plot <- renderPlot({
    req(input$sample_file)
    df  <- req(sample_data())
    col <- req(input$hist_col)
    # Normalise NULL / empty string to "None" so the helper function
    # receives a consistent sentinel value
    grp <- if (is.null(input$group_col) || input$group_col == "") "None"
    else input$group_col
    plot_metadata_distribution(df, col, grp)
  })
  
  
  # ── TAB 2: Counts explorer ────────────────────────────────────
  
  # Read the normalized counts matrix
  counts_data <- reactive({
    req(input$counts_file)
    df <- read.csv(input$counts_file$datapath,
                   stringsAsFactors = FALSE, row.names = 1)
    validate(need(is.numeric(df[, 1]), "Count columns must be numeric"))
    df
  })
  
  # Apply variance-percentile and min-non-zero filters from the sidebar sliders
  filtered_counts_data <- reactive({
    filter_counts(req(counts_data()), input$var_pct, input$nonzero_n)
  })
  
  # Build the X-axis PC selector dynamically once a counts file is loaded.
  output$pc_x_ui <- renderUI({
    req(input$counts_file)
    n_pcs <- min(10, ncol(req(counts_data())) - 1)
    selectInput("pc_x", "X axis (PC)",
                choices = paste0("PC", 1:n_pcs), selected = "PC1")
  })
  
  # Same logic for the Y-axis PC selector; defaults to PC2.
  output$pc_y_ui <- renderUI({
    req(input$counts_file)
    n_pcs <- min(10, ncol(req(counts_data())) - 1)
    selectInput("pc_y", "Y axis (PC)",
                choices = paste0("PC", 1:n_pcs), selected = "PC2")
  })
  
  # Render four stat boxes summarising the filter effect:
  #   samples | total genes | passing genes (%) | filtered genes (%)
  output$filter_stats <- renderUI({
    req(input$counts_file)
    total   <- nrow(req(counts_data()))
    passing <- nrow(filtered_counts_data())
    fluidRow(
      column(3, div(class = "stat-box",
                    div(class = "stat-val", ncol(counts_data())),
                    div(class = "stat-lbl", "samples"))),
      column(3, div(class = "stat-box",
                    div(class = "stat-val", format(total, big.mark = ",")),
                    div(class = "stat-lbl", "total genes"))),
      column(3, div(class = "stat-box",
                    div(class = "stat-val", format(passing, big.mark = ",")),
                    div(class = "stat-lbl",
                        paste0("passing (",
                               round(100 * passing / total, 1), "%)")))),
      column(3, div(class = "stat-box",
                    div(class = "stat-val", format(total - passing, big.mark = ",")),
                    div(class = "stat-lbl",
                        paste0("filtered (",
                               round(100 * (total - passing) / total, 1), "%)"))))
    )
  })
  
  # Tabular version of the same filter summary (easier to copy/paste).
  output$filter_summary_tbl <- renderTable({
    req(input$counts_file)
    total   <- nrow(req(counts_data()))
    passing <- nrow(filtered_counts_data())
    data.frame(
      Metric  = c("Samples", "Total genes", "Genes passing", "Genes filtered"),
      Count   = c(ncol(counts_data()), total, passing, total - passing),
      Percent = c("--", "100%",
                  paste0(round(100 * passing / total, 1), "%"),
                  paste0(round(100 * (total - passing) / total, 1), "%"))
    )
  }, striped = TRUE, bordered = TRUE)
  
  # Scatter plot 1: log10(median count + 1) vs log10(variance + 1).
  # Rose = passing genes; grey = filtered genes.
  output$scatter_var <- renderPlot({
    req(input$counts_file)
    plot_scatter_variance(req(counts_data()), filtered_counts_data())
  })
  
  # Scatter plot 2: log10(median count + 1) vs number of zero-count samples.
  output$scatter_zeros <- renderPlot({
    req(input$counts_file)
    plot_scatter_zeros(req(counts_data()), filtered_counts_data())
  })
  
  # Clustered heatmap of the top N most variable genes after filtering.
  output$heatmap_plot <- renderPlot({
    req(input$counts_file)
    plot_heatmap(req(filtered_counts_data()),
                 top_n = input$hm_topn, log_transform = input$hm_log)
  })
  
  # Run PCA once on the filtered, scaled matrix and cache the result.
  pca_result <- reactive({ run_pca(req(filtered_counts_data())) })
  
  # Render the PCA scatter plot
  output$pca_plot <- renderPlot({
    req(input$counts_file)
    meta <- NULL
    if (isTRUE(input$pca_color_meta) && !is.null(input$pca_meta_file)) {
      meta <- tryCatch(
        read.csv(input$pca_meta_file$datapath, stringsAsFactors = FALSE),
        error = function(e) NULL   # return NULL on any read error
      )
    }
    plot_pca(pca_result(),
             pc_x = req(input$pc_x),
             pc_y = req(input$pc_y),
             meta = meta)
  })
  
  
  # ── TAB 3: Differential expression ───────────────────────────
  
  # Read the DESeq2 results CSV. validate(need()) checks for the three
  # columns required by downstream functions and displays a clear inline
  # message if any are missing.
  de_data <- reactive({
    req(input$de_file)
    df <- read.csv(input$de_file$datapath, stringsAsFactors = FALSE)
    validate(need("gene"           %in% names(df), "Missing column: gene"))
    validate(need("log2FoldChange" %in% names(df), "Missing column: log2FoldChange"))
    validate(need("padj"           %in% names(df), "Missing column: padj"))
    df
  })
  
  # Add a "direction" column (Up / Down / NS) to each gene based on the
  # current padj and |log2FC| slider values.
  de_labeled <- reactive({
    label_de_direction(req(de_data()), input$padj_thresh, input$lfc_thresh)
  })
  
  # Render the up/down-regulated count badges in the sidebar.
  output$de_sig_counts <- renderUI({
    req(input$de_file)
    df   <- req(de_labeled())
    n_up <- sum(df$direction == "Up")
    n_dn <- sum(df$direction == "Down")
    tagList(
      tags$h5("Significant genes",
              style = "color:#72243E; font-weight:700; margin-bottom:8px;"),
      div(class = "dir-badge badge-up",
          tags$span("\u2191 Upregulated", style = "display:block; font-size:11px;"),
          tags$span(n_up, style = "font-size:20px; font-weight:700;")),
      div(class = "dir-badge badge-down",
          tags$span("\u2193 Downregulated", style = "display:block; font-size:11px;"),
          tags$span(n_dn, style = "font-size:20px; font-weight:700;"))
    )
  })
  
  # Render the volcano plot
  output$volcano_plot <- renderPlot({
    req(input$de_file)
    plot_volcano(req(de_labeled()), input$padj_thresh, input$lfc_thresh)
  })
  
  # Render the full DE results as an interactive DT table.
  output$de_dt <- DT::renderDataTable({
    req(input$de_file)
    df <- req(de_data()) %>%
      mutate(across(where(is.numeric), ~ round(.x, 4)))
    DT::datatable(df, filter = "top", rownames = FALSE,
                  options = list(pageLength = 20, scrollX = TRUE)) %>%
      DT::formatStyle("log2FoldChange",
                      backgroundColor = DT::styleInterval(
                        c(-1, 1), c("#F4C0D1", "white", "#FBEAF0")))
  })
  
  
  # ── TAB 4: GSEA ──────────────────────────────────────────────
  
  # Read the fgsea results CSV
  gsea_data <- reactive({
    req(input$gsea_file)
    df <- read.csv(input$gsea_file$datapath, stringsAsFactors = FALSE)
    validate(need("pathway" %in% names(df), "Missing column: pathway"))
    validate(need("NES"     %in% names(df), "Missing column: NES"))
    validate(need("padj"    %in% names(df), "Missing column: padj"))
    df
  })
  
  # Prepare the top N pathways for the barplot
  gsea_top_df <- reactive({
    prep_gsea_barplot_df(req(gsea_data()), input$gsea_top_n)
  })
  
  # Render the NES horizontal bar chart.
  output$gsea_barplot <- renderPlot({
    req(input$gsea_file)
    plot_gsea_barplot(req(gsea_top_df()), input$gsea_top_n)
  })
  
  # Handle barplot click events to display a detail card for the clicked pathway.
  output$gsea_clicked_info <- renderUI({
    click <- input$gsea_bar_click
    if (is.null(click)) return(NULL)   # nothing rendered before first click
    df    <- req(gsea_top_df())
    y_idx <- max(1, min(nrow(df), round(click$y)))
    row   <- df[y_idx, ]
    # Render a styled info card with NES, padj, pval, and gene set size
    div(style = "background:#FFF0F5; border:1px solid #F4C0D1; border-radius:10px;
                 padding:14px; margin-top:8px;",
        tags$h5(style = "color:#72243E; margin-bottom:10px; font-weight:700;",
                toupper(as.character(row$pathway_short))),
        fluidRow(
          column(3, div(class = "stat-box",
                        div(class = "stat-val", round(row$NES, 3)),
                        div(class = "stat-lbl", "NES"))),
          column(3, div(class = "stat-box",
                        div(class = "stat-val",
                            format(row$padj, scientific = TRUE, digits = 2)),
                        div(class = "stat-lbl", "padj"))),
          column(3, div(class = "stat-box",
                        div(class = "stat-val",
                            format(row$pval, scientific = TRUE, digits = 2)),
                        div(class = "stat-lbl", "pval"))),
          column(3, div(class = "stat-box",
                        div(class = "stat-val", row$size),
                        div(class = "stat-lbl", "size")))
        )
    )
  })
  
  # Filter the full GSEA results by the sidebar padj slider and NES
  gsea_filtered <- reactive({
    filter_gsea_table(req(gsea_data()), input$gsea_padj_tbl, input$gsea_nes_filter)
  })
  
  # Render the filtered GSEA results as an interactive DT table.
  output$gsea_dt <- DT::renderDataTable({
    req(input$gsea_file)
    DT::datatable(gsea_filtered(), filter = "top", rownames = FALSE,
                  options = list(pageLength = 15, scrollX = TRUE)) %>%
      DT::formatStyle("NES",
                      backgroundColor = DT::styleInterval(
                        0, c("#F4C0D1", "#FBEAF0")))
  })
  
  # Provide a CSV download of whatever gsea_filtered() currently contains.
  output$gsea_download <- downloadHandler(
    filename = function() paste0("gsea_filtered_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(gsea_filtered(), file, row.names = FALSE)
  )
  
  # Render the NES vs -log10(padj) scatter plot.
  output$gsea_scatter <- renderPlot({
    req(input$gsea_file)
    plot_gsea_scatter(req(gsea_data()), input$gsea_padj_scatter)
  })
  
}

shinyApp(ui = ui, server = server)