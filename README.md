# Huntington-s-Disease-Transcriptomic-Explorer
Exploration of scRNA-seq data from post-mortem Huntington's Disease patients

**This is my final project for BI 591 - Biological Data Science in R.**

I used **GEO dataset GSE64810** which profiled gene expression from 69 human subjects: 
20 who died from Huntington's Disease, and 49 neurologically normal controls using prefrontal cortex tissue. 

My submission is within the submission folder where there are 3 files relevant to this project:

  **1. preprocessing.Rmd --- takes the raw data and creates 4 csv**
  
       metadata.csv              | sample metadata (69 samples) 
       normalized_counts.csv     | DESeq2 normalized counts matrix
       fsgea_results.csv         | fgsea Hallwark pathway enrichment results
       deseq2_results.csv        | DESeq2 DE results (HD vs control)
       
  **2. app.R --- R shiny that handles interactive components of the app**
  
  **4. utils2.R --- stores helper functions that processes csv, filter, and plot**

**File structure**

**metadata.csv:**

-> sample_id = GEO sample accession

-> tissue = BA9 prefrontal cortex

-> diagnosis = Huntington's Disease or Neurologically normal

-> pmi = post-mortem interval (hours)

-> age_of_death = age at death (years)

-> rin = RNA integrity number

**normalized_counts.csv**

-> Genes x samples normalized counts using DESeq2 median-of-ratios normalization

**deseq2_results.csv**

-> gene = Ensembl ID

-> baseMean = mean normalized count across all samples

-> log2FoldChange = log2 fold change (HD/control)

-> lfcSE = standard error of log2 fold change

-> stat = Wald statistic

-> pvalue = raw p-value

-> padj = BH-adjusted p-value

**App structure**

The app is organized as four navbar tabs, each focused on a different stage of the analysis.

**Tab 1 — Sample Information**

Upload and explore the sample metadata CSV. Displays a per-column type and summary statistics table, a full searchable/sortable data table, and a distribution plot (histogram or grouped violin) for any numeric column.

**Tab 2 — Counts Explorer**

Upload a normalized counts matrix and interactively filter genes by variance percentile and minimum non-zero sample count. Explore the filtered data through:

  Filter summary — stat boxes and table showing passing vs filtered gene counts

  Scatter plots — median vs variance and median vs zero count, with passing genes highlighted

  Heatmap — hierarchically clustered heatmap of the top N most variable genes (z-scored)

  PCA — principal component scatter plot, optionally colored by diagnosis

**Tab 3 — Differential Expression**

Upload DESeq2 results and explore them with interactive padj and log2 fold change threshold sliders. Includes a volcano plot with automatic gene labeling and a full sortable results table with fold-change cell shading.

**Tab 4 — GSEA**

Upload fgsea pathway enrichment results and explore them through:

  Barplot — NES bar chart of top N pathways; click any bar for a detail card

  Table — filterable table by padj and NES direction, with CSV download

  Scatter — NES vs −log10(padj) with significance threshold highlighting
