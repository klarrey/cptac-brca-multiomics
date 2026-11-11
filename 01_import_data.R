# CPTAC breast cancer: transcriptomic and proteomic analysis
# Run from the project root, after installing dependencies (see README.md).
# Input data instructions: docs/data_access.md
# Positive subtype coefficients mean Basal minus LumA.

required_packages <- c("limma", "ggplot2", "ggrepel")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}

#inspect the first few lines of the three starting files
# Paths to the three starting datasets
clinical_file <- "data/HS_CPTAC_BRCA_2018_CLI.tsi"

rna_file <- "data/HS_CPTAC_BRCA_2018_RNA_GENE.cct"

protein_file <- "data/HS_CPTAC_BRCA_2018_Proteome_Ratio_Norm_gene_Median.cct"

# Stop early if the downloaded inputs are unavailable
input_files <- c(clinical_file, rna_file, protein_file)
if (!all(file.exists(input_files))) {
  stop("Missing input files: ", paste(input_files[!file.exists(input_files)], collapse = ", "),
       ". See docs/data_access.md.")
}

# Preview the first three lines of each file
readLines(clinical_file, n = 3)
readLines(rna_file, n = 3)
readLines(protein_file, n = 3)


# Import clinical data, including its variable-type annotation row
clinical_import <- read.delim(
  clinical_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)

# Preserve the variable-type information separately
clinical_variable_types <- clinical_import[1, , drop = FALSE]

# Keep only patient records
clinical <- clinical_import[-1, , drop = FALSE]
rownames(clinical) <- NULL

# Import RNA expression data
rna <- read.delim(
  rna_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)

# Import protein abundance data
protein <- read.delim(
  protein_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)


# Dimensions: number of rows, followed by number of columns
dim(clinical)
dim(rna)
dim(protein)

# Preview a small section of each table
clinical[1:4, 1:6]
rna[1:4, 1:5]
protein[1:4, 1:5]

# Count samples in each molecular dataset,
# excluding the first column containing gene identifiers
ncol(rna) - 1
ncol(protein) - 1

# Examine available breast cancer subtype labels
table(clinical$PAM50, useNA = "ifany")


# Extract sample identifiers
clinical_ids <- clinical$Sample.ID
rna_ids <- colnames(rna)[-1]
protein_ids <- colnames(protein)[-1]

# Check for duplicated sample identifiers
c(
  clinical = anyDuplicated(clinical_ids),
  rna = anyDuplicated(rna_ids),
  protein = anyDuplicated(protein_ids)
)

# Check whether the datasets contain the same samples
c(
  clinical_vs_rna = setequal(clinical_ids, rna_ids),
  clinical_vs_protein = setequal(clinical_ids, protein_ids)
)

# Stop if identifiers are missing, blank, duplicated, or unmatched
stopifnot(
  !anyNA(clinical_ids),
  !anyNA(rna_ids),
  !anyNA(protein_ids),
  all(nzchar(clinical_ids)),
  all(nzchar(rna_ids)),
  all(nzchar(protein_ids)),
  anyDuplicated(clinical_ids) == 0L,
  anyDuplicated(rna_ids) == 0L,
  anyDuplicated(protein_ids) == 0L,
  setequal(clinical_ids, rna_ids),
  setequal(clinical_ids, protein_ids)
)

# Use the RNA sample order as our reference
clinical_aligned <- clinical[
  match(rna_ids, clinical$Sample.ID),
  ,
  drop = FALSE
]

# Keep the gene identifier column, then reorder protein sample columns
protein_aligned <- protein[
  ,
  c(colnames(protein)[1], rna_ids),
  drop = FALSE
]

rna_aligned <- rna

# Confirm identical sample order
c(
  clinical_order = identical(
    clinical_aligned$Sample.ID,
    colnames(rna_aligned)[-1]
  ),
  protein_order = identical(
    colnames(protein_aligned)[-1],
    colnames(rna_aligned)[-1]
  )
)

# Check for duplicated gene identifiers before assigning row names
c(
  rna_duplicate_genes = sum(duplicated(rna_aligned$IDX)),
  protein_duplicate_genes = sum(duplicated(protein_aligned$IDX))
)


# Confirm that gene identifiers are present and nonblank
stopifnot(
  anyDuplicated(rna_aligned$IDX) == 0L,
  anyDuplicated(protein_aligned$IDX) == 0L,
  !anyNA(rna_aligned$IDX),
  !anyNA(protein_aligned$IDX),
  all(nzchar(rna_aligned$IDX)),
  all(nzchar(protein_aligned$IDX))
)

# Confirm that every expression column is numeric
stopifnot(
  all(vapply(rna_aligned[-1], is.numeric, logical(1))),
  all(vapply(protein_aligned[-1], is.numeric, logical(1)))
)

# Create matrices containing expression values only
rna_mat <- as.matrix(rna_aligned[, -1, drop = FALSE])
protein_mat <- as.matrix(protein_aligned[, -1, drop = FALSE])

# Use gene identifiers as row names
rownames(rna_mat) <- rna_aligned$IDX
rownames(protein_mat) <- protein_aligned$IDX

# Confirm that sample alignment was preserved
stopifnot(
  identical(colnames(rna_mat), clinical_aligned$Sample.ID),
  identical(colnames(protein_mat), clinical_aligned$Sample.ID)
)


# Overall missingness and completely missing genes
missing_summary <- data.frame(
  dataset = c("RNA", "Protein"),
  genes = c(nrow(rna_mat), nrow(protein_mat)),
  samples = c(ncol(rna_mat), ncol(protein_mat)),
  missing_percent = c(
    100 * mean(is.na(rna_mat)),
    100 * mean(is.na(protein_mat))
  ),
  entirely_missing_genes = c(
    sum(rowSums(!is.na(rna_mat)) == 0),
    sum(rowSums(!is.na(protein_mat)) == 0)
  ),
  infinite_values = c(
    sum(is.infinite(rna_mat)),
    sum(is.infinite(protein_mat))
  )
)

missing_summary

# Missing fraction for each gene
rna_gene_missing <- rowMeans(is.na(rna_mat))
protein_gene_missing <- rowMeans(is.na(protein_mat))

# Missing fraction for each sample
rna_sample_missing <- colMeans(is.na(rna_mat))
protein_sample_missing <- colMeans(is.na(protein_mat))

# Summarize the distributions
summary(rna_gene_missing)
summary(protein_gene_missing)

summary(rna_sample_missing)
summary(protein_sample_missing)


# Require measurements in at least 80% of samples
min_observed_fraction <- 0.80

# Calculate the minimum number of observed measurements
min_observed_samples <- ceiling(
  min_observed_fraction * ncol(rna_mat)
)

min_observed_samples

# Identify genes meeting the requirement
keep_rna <- rowSums(!is.na(rna_mat)) >= min_observed_samples
keep_protein <- rowSums(!is.na(protein_mat)) >= min_observed_samples

# Create filtered matrices
rna_filtered <- rna_mat[keep_rna, , drop = FALSE]
protein_filtered <- protein_mat[keep_protein, , drop = FALSE]

# Summarize the effect of filtering
filter_summary <- data.frame(
  dataset = c("RNA", "Protein"),
  genes_before = c(nrow(rna_mat), nrow(protein_mat)),
  genes_retained = c(nrow(rna_filtered), nrow(protein_filtered)),
  genes_removed = c(sum(!keep_rna), sum(!keep_protein)),
  missing_percent_after = c(
    100 * mean(is.na(rna_filtered)),
    100 * mean(is.na(protein_filtered))
  )
)

filter_summary

# Check sample missingness after filtering
summary(colMeans(is.na(rna_filtered)))
summary(colMeans(is.na(protein_filtered)))

# Count retained genes shared between the two layers
shared_genes <- intersect(
  rownames(rna_filtered),
  rownames(protein_filtered)
)

length(shared_genes)


#check whether the remaining missingness differs by subtype and inspect sample 
#expression distributions
# Confirm alignment before attaching clinical labels
stopifnot(
  identical(colnames(rna_filtered), clinical_aligned$Sample.ID),
  identical(colnames(protein_filtered), clinical_aligned$Sample.ID)
)

# Create a sample-level quality-control table
sample_qc <- data.frame(
  sample_id = clinical_aligned$Sample.ID,
  subtype = clinical_aligned$PAM50,
  rna_missing_percent = 100 * colMeans(is.na(rna_filtered)),
  protein_missing_percent = 100 * colMeans(is.na(protein_filtered)),
  rna_median = apply(rna_filtered, 2, median, na.rm = TRUE),
  protein_median = apply(protein_filtered, 2, median, na.rm = TRUE)
)

# Average missingness within each subtype
aggregate(
  cbind(rna_missing_percent, protein_missing_percent) ~ subtype,
  data = sample_qc,
  FUN = mean
)

# Count retained genes shared by RNA and protein
shared_genes <- intersect(
  rownames(rna_filtered),
  rownames(protein_filtered)
)

length(shared_genes)


# Arrange four plots in one figure
old_par <- par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))

boxplot(
  rna_missing_percent ~ subtype,
  data = sample_qc,
  main = "RNA missingness by subtype",
  ylab = "Missing measurements (%)",
  xlab = "",
  las = 2
)

boxplot(
  protein_missing_percent ~ subtype,
  data = sample_qc,
  main = "Protein missingness by subtype",
  ylab = "Missing measurements (%)",
  xlab = "",
  las = 2
)

boxplot(
  as.data.frame(rna_filtered),
  outline = FALSE,
  names = rep("", ncol(rna_filtered)),
  main = "RNA distributions by sample",
  ylab = "Normalized log2 RNA expression",
  xlab = "Samples"
)

boxplot(
  as.data.frame(protein_filtered),
  outline = FALSE,
  names = rep("", ncol(protein_filtered)),
  main = "Protein distributions by sample",
  ylab = "Normalized log2 protein ratio",
  xlab = "Samples"
)

# Restore the previous plotting settings
par(old_par)


# Select genes with measurements in every sample
rna_complete <- rna_filtered[
  rowSums(is.na(rna_filtered)) == 0,
  ,
  drop = FALSE
]

protein_complete <- protein_filtered[
  rowSums(is.na(protein_filtered)) == 0,
  ,
  drop = FALSE
]

# Remove genes with no variation across samples
rna_variance <- apply(rna_complete, 1, var)
protein_variance <- apply(protein_complete, 1, var)

rna_pca_input <- rna_complete[
  is.finite(rna_variance) & rna_variance > 0,
  ,
  drop = FALSE
]

protein_pca_input <- protein_complete[
  is.finite(protein_variance) & protein_variance > 0,
  ,
  drop = FALSE
]

# Report how many genes will contribute
c(
  RNA_genes = nrow(rna_pca_input),
  Protein_genes = nrow(protein_pca_input)
)

# PCA requires samples in rows and genes in columns
rna_pca <- prcomp(
  t(rna_pca_input),
  center = TRUE,
  scale. = FALSE
)

protein_pca <- prcomp(
  t(protein_pca_input),
  center = TRUE,
  scale. = FALSE
)

# Percentage of total variance explained by each component
rna_pve <- 100 * rna_pca$sdev^2 / sum(rna_pca$sdev^2)
protein_pve <- 100 * protein_pca$sdev^2 / sum(protein_pca$sdev^2)

round(rna_pve[1:5], 2)
round(protein_pve[1:5], 2)


# Define colors shared by both plots
subtype_colors <- c(
  Basal = "#D55E00",
  Her2 = "#CC79A7",
  LumA = "#0072B2",
  LumB = "#009E73",
  "Normal-like" = "#777777"
)

point_colors <- unname(
  subtype_colors[as.character(clinical_aligned$PAM50)]
)

# Confirm sample order and color assignments
stopifnot(
  identical(rownames(rna_pca$x), clinical_aligned$Sample.ID),
  identical(rownames(protein_pca$x), clinical_aligned$Sample.ID),
  !anyNA(point_colors)
)

old_par <- par(mfrow = c(1, 2))

plot(
  rna_pca$x[, 1:2],
  col = point_colors,
  pch = 19,
  xlab = sprintf("PC1 (%.1f%%)", rna_pve[1]),
  ylab = sprintf("PC2 (%.1f%%)", rna_pve[2]),
  main = "RNA PCA"
)

legend(
  "topright",
  legend = names(subtype_colors),
  col = subtype_colors,
  pch = 19,
  bty = "n",
  cex = 0.75
)

plot(
  protein_pca$x[, 1:2],
  col = point_colors,
  pch = 19,
  xlab = sprintf("PC1 (%.1f%%)", protein_pve[1]),
  ylab = sprintf("PC2 (%.1f%%)", protein_pve[2]),
  main = "Protein PCA"
)

legend(
  "topright",
  legend = names(subtype_colors),
  col = subtype_colors,
  pch = 19,
  bty = "n",
  cex = 0.75
)

par(old_par)


# Select Basal and LumA tumors
selected_samples <- clinical_aligned$PAM50 %in% c("LumA", "Basal")

clinical_compare <- clinical_aligned[
  selected_samples, , drop = FALSE
]

# Set LumA as the reference group
clinical_compare$subtype <- factor(
  clinical_compare$PAM50,
  levels = c("LumA", "Basal")
)

# Use the filtered matrices, which still retain their missing values
rna_compare <- rna_filtered[
  , selected_samples, drop = FALSE
]

protein_compare <- protein_filtered[
  , selected_samples, drop = FALSE
]

# Confirm sample order
stopifnot(
  identical(colnames(rna_compare), clinical_compare$Sample.ID),
  identical(colnames(protein_compare), clinical_compare$Sample.ID)
)

table(clinical_compare$subtype)

# Inspect stage distribution and missing stage information
table(
  subtype = clinical_compare$subtype,
  stage = clinical_compare$Stage,
  useNA = "ifany"
)

# List available clinical variables
names(clinical_compare)


# Identify samples belonging to each group
is_luma <- clinical_compare$subtype == "LumA"
is_basal <- clinical_compare$subtype == "Basal"

# Require at least 80% observed measurements within EACH group
min_luma <- ceiling(0.80 * sum(is_luma))
min_basal <- ceiling(0.80 * sum(is_basal))

c(min_LumA = min_luma, min_Basal = min_basal)

# Identify genes meeting both group-specific requirements
rna_group_keep <-
  rowSums(!is.na(rna_compare[, is_luma, drop = FALSE])) >= min_luma &
  rowSums(!is.na(rna_compare[, is_basal, drop = FALSE])) >= min_basal

protein_group_keep <-
  rowSums(!is.na(protein_compare[, is_luma, drop = FALSE])) >= min_luma &
  rowSums(!is.na(protein_compare[, is_basal, drop = FALSE])) >= min_basal

# Report coverage before applying this additional filter
data.frame(
  dataset = c("RNA", "Protein"),
  genes_before = c(nrow(rna_compare), nrow(protein_compare)),
  passing_both_groups = c(
    sum(rna_group_keep),
    sum(protein_group_keep)
  ),
  failing_group_coverage = c(
    sum(!rna_group_keep),
    sum(!protein_group_keep)
  )
)


# Apply the group-specific coverage filters
rna_de <- rna_compare[rna_group_keep, , drop = FALSE]
protein_de <- protein_compare[protein_group_keep, , drop = FALSE]

# Convert clinical measurements from text to numbers
# Unexpected nonnumeric values will generate a warning
clinical_compare$age_years <-
  as.numeric(clinical_compare$Age.in.Month) / 12

clinical_compare$tumor_cellularity <-
  as.numeric(clinical_compare$Tumor.Cellularity)

# Summarize availability and values within each subtype
by(
  clinical_compare[
    , c("age_years", "tumor_cellularity"),
    drop = FALSE
  ],
  clinical_compare$subtype,
  summary
)

# Inspect the study's existing simplified stage annotation
table(
  subtype = clinical_compare$subtype,
  stage_reduced = clinical_compare$Stage.reduced,
  useNA = "ifany"
)

# Check how many samples have both age and cellularity recorded
table(
  subtype = clinical_compare$subtype,
  complete_age_cellularity = complete.cases(
    clinical_compare[, c("age_years", "tumor_cellularity")]
  )
)

# Confirm the dimensions of the differential-analysis matrices
dim(rna_de)
dim(protein_de)


#Build the design matrix
library("limma")
# Confirm the reference group and sample order
stopifnot(
  identical(levels(clinical_compare$subtype), c("LumA", "Basal")),
  identical(colnames(rna_de), clinical_compare$Sample.ID),
  identical(colnames(protein_de), clinical_compare$Sample.ID)
)

# Specify the model: expression depends on subtype
design <- model.matrix(
  ~ subtype,
  data = clinical_compare
)

rownames(design) <- clinical_compare$Sample.ID

# Inspect how the groups are encoded
head(design)
colnames(design)

# Confirm that the model coefficients can be estimated
stopifnot(qr(design)$rank == ncol(design))


# Fit a separate linear model for every gene
rna_fit <- limma::lmFit(rna_de, design)
protein_fit <- limma::lmFit(protein_de, design)

# Moderate the gene-wise variance estimates
rna_fit <- limma::eBayes(rna_fit, trend = FALSE)
protein_fit <- limma::eBayes(protein_fit, trend = FALSE)

# Extract all results, ranked by statistical evidence
rna_results <- limma::topTable(
  rna_fit,
  coef = "subtypeBasal",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

protein_results <- limma::topTable(
  protein_fit,
  coef = "subtypeBasal",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

# Preserve gene identifiers as an explicit column
rna_results <- data.frame(
  gene = rownames(rna_results),
  rna_results,
  row.names = NULL
)

protein_results <- data.frame(
  gene = rownames(protein_results),
  protein_results,
  row.names = NULL
)

head(rna_results)
head(protein_results)


fdr_threshold <- 0.05

data.frame(
  dataset = c("RNA", "Protein"),
  genes_reported = c(nrow(rna_results), nrow(protein_results)),
  FDR_below_0.05 = c(
    sum(rna_results$adj.P.Val < fdr_threshold, na.rm = TRUE),
    sum(protein_results$adj.P.Val < fdr_threshold, na.rm = TRUE)
  ),
  missing_adjusted_pvalues = c(
    sum(is.na(rna_results$adj.P.Val)),
    sum(is.na(protein_results$adj.P.Val))
  )
)

# Inspect residual variability against average normalized expression
old_par <- par(mfrow = c(1, 2))

limma::plotSA(rna_fit, main = "RNA variance diagnostic")
limma::plotSA(protein_fit, main = "Protein variance diagnostic")

par(old_par)


# Define thresholds outside the function
fdr_threshold <- 0.05
rna_effect_threshold <- 1
protein_effect_threshold <- 1

# Reusable function to classify differential results
classify_results <- function(results, fdr_cutoff, effect_cutoff) {
  
  # Start with a copy of the results table
  results$status <- "Below criteria"
  
  # Identify rows with usable statistics
  valid <- is.finite(results$logFC) &
    !is.na(results$adj.P.Val)
  
  # Higher in Basal
  higher_basal <- valid &
    results$adj.P.Val < fdr_cutoff &
    results$logFC >= effect_cutoff
  
  # Higher in LumA
  higher_luma <- valid &
    results$adj.P.Val < fdr_cutoff &
    results$logFC <= -effect_cutoff
  
  results$status[higher_basal] <- "Higher in Basal"
  results$status[higher_luma] <- "Higher in LumA"
  results$status[!valid] <- "Not evaluated"
  
  results$status <- factor(
    results$status,
    levels = c(
      "Higher in LumA",
      "Below criteria",
      "Higher in Basal",
      "Not evaluated"
    )
  )
  
  return(results)
}

# Apply the same classification logic to both datasets
rna_results <- classify_results(
  rna_results,
  fdr_cutoff = fdr_threshold,
  effect_cutoff = rna_effect_threshold
)

protein_results <- classify_results(
  protein_results,
  fdr_cutoff = fdr_threshold,
  effect_cutoff = protein_effect_threshold
)

# Count genes in each category
table(rna_results$status)
table(protein_results$status)

# Inspect the highest-ranked genes meeting both criteria
head(
  subset(
    rna_results,
    status %in% c("Higher in Basal", "Higher in LumA")
  ),
  10
)

head(
  subset(
    protein_results,
    status %in% c("Higher in Basal", "Higher in LumA")
  ),
  10
)


library("ggplot2")
library("ggrepel")

make_volcano <- function(results, title, effect_cutoff,
                         fdr_cutoff = 0.05, top_n = 10) {
  
  # Keep rows with valid plotting values
  plot_data <- results[
    is.finite(results$logFC) &
      is.finite(results$adj.P.Val) &
      results$adj.P.Val >= 0 &
      results$adj.P.Val <= 1,
    ,
    drop = FALSE
  ]
  
  # Apply the supplied thresholds consistently
  plot_data <- classify_results(
    plot_data,
    fdr_cutoff = fdr_cutoff,
    effect_cutoff = effect_cutoff
  )
  
  # Avoid infinite plotting coordinates if a p-value is stored as zero
  plot_data$neg_log10_fdr <- -log10(
    pmax(plot_data$adj.P.Val, .Machine$double.xmin)
  )
  
  # Select genes for labels
  label_data <- plot_data[
    plot_data$status %in% c("Higher in LumA", "Higher in Basal"),
    ,
    drop = FALSE
  ]
  
  label_data <- label_data[
    order(label_data$adj.P.Val, label_data$gene),
    ,
    drop = FALSE
  ]
  
  label_data <- head(label_data, top_n)
  
  # Symmetric horizontal limits
  x_limit <- 1.1 * max(
    abs(plot_data$logFC),
    effect_cutoff
  )
  
  counts <- table(plot_data$status)
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = logFC, y = neg_log10_fdr, color = status)
  ) +
    ggplot2::geom_point(size = 1.3, alpha = 0.6) +
    ggplot2::geom_vline(
      xintercept = c(-effect_cutoff, effect_cutoff),
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_cutoff),
      linetype = "dashed",
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      seed = 42,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Higher in LumA" = "blue",
        "Below criteria" = "grey75",
        "Higher in Basal" = "red"
      ),
      drop = TRUE
    ) +
    ggplot2::coord_cartesian(xlim = c(-x_limit, x_limit)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.15))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf(
        "Higher in LumA: %d | Higher in Basal: %d | Below criteria: %d",
        counts["Higher in LumA"],
        counts["Higher in Basal"],
        counts["Below criteria"]
      ),
      x = "Log2 effect (Basal minus LumA)",
      y = expression(-log[10]("BH-adjusted p-value")),
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}


rna_volcano <- make_volcano(
  rna_results,
  title = "RNA: Basal versus LumA",
  effect_cutoff = rna_effect_threshold,
  fdr_cutoff = fdr_threshold
)

protein_volcano <- make_volcano(
  protein_results,
  title = "Protein: Basal versus LumA",
  effect_cutoff = protein_effect_threshold,
  fdr_cutoff = fdr_threshold
)

print(rna_volcano)
print(protein_volcano)


# Select the columns needed for integration
rna_effects <- rna_results[
  , c("gene", "logFC", "adj.P.Val", "status")
]

protein_effects <- protein_results[
  , c("gene", "logFC", "adj.P.Val", "status")
]

# Give each layer's columns distinct names
colnames(rna_effects) <- c(
  "gene", "rna_logFC", "rna_fdr", "rna_status"
)

colnames(protein_effects) <- c(
  "gene", "protein_logFC", "protein_fdr", "protein_status"
)

# Match results by gene identifier
integrated_results <- merge(
  rna_effects,
  protein_effects,
  by = "gene",
  all = FALSE
)

# Check the merged table
nrow(integrated_results)
head(integrated_results)

stopifnot(
  anyDuplicated(integrated_results$gene) == 0L
)


# Identify genes with finite effects in both layers
valid_effects <- with(
  integrated_results,
  is.finite(rna_logFC) & is.finite(protein_logFC)
)

# Compare effect sizes across genes
effect_correlations <- c(
  Pearson = cor(
    integrated_results$rna_logFC[valid_effects],
    integrated_results$protein_logFC[valid_effects],
    method = "pearson"
  ),
  Spearman = cor(
    integrated_results$rna_logFC[valid_effects],
    integrated_results$protein_logFC[valid_effects],
    method = "spearman"
  )
)

effect_correlations

# Cross-tabulate the classifications among shared tested genes
with(
  integrated_results,
  table(RNA = rna_status, Protein = protein_status)
)


#plot the effects
effect_plot <- ggplot2::ggplot(
  integrated_results[valid_effects, ],
  ggplot2::aes(x = rna_logFC, y = protein_logFC)
) +
  ggplot2::geom_point(
    alpha = 0.3,
    size = 1,
    color = "#336699"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  ggplot2::labs(
    title = "RNA–protein agreement in subtype effects",
    subtitle = sprintf(
      "Shared genes: %d | Pearson r = %.2f | Spearman rho = %.2f",
      sum(valid_effects),
      effect_correlations["Pearson"],
      effect_correlations["Spearman"]
    ),
    x = "RNA log2 effect (Basal minus LumA)",
    y = "Protein log2 effect (Basal minus LumA)"
  ) +
  ggplot2::theme_classic(base_size = 12)

print(effect_plot)


# Identify genes meeting our criteria in each layer
rna_pass <- integrated_results$rna_status %in%
  c("Higher in LumA", "Higher in Basal")

protein_pass <- integrated_results$protein_status %in%
  c("Higher in LumA", "Higher in Basal")

# Find opposite directions among genes passing in BOTH layers
opposite_effects <- with(
  integrated_results,
  sign(rna_logFC) != sign(protein_logFC)
)

discordant_genes <- integrated_results[
  rna_pass & protein_pass & opposite_effects,
  ,
  drop = FALSE
]

# Rank by the weaker statistical evidence of the two layers
discordant_genes <- discordant_genes[
  order(pmax(
    discordant_genes$rna_fdr,
    discordant_genes$protein_fdr
  )),
  ,
  drop = FALSE
]

# Display all 12 genes
discordant_genes[
  , c(
    "gene",
    "rna_logFC", "rna_fdr",
    "protein_logFC", "protein_fdr"
  )
]


# Verify that both matrices match the clinical sample order
stopifnot(
  identical(colnames(rna_de), clinical_compare$Sample.ID),
  identical(colnames(protein_de), clinical_compare$Sample.ID)
)

# Extract measurements for each discordant gene
discordant_values <- do.call(
  rbind,
  lapply(discordant_genes$gene, function(g) {
    
    rbind(
      data.frame(
        gene = g,
        sample_id = clinical_compare$Sample.ID,
        subtype = clinical_compare$subtype,
        layer = "RNA",
        value = as.numeric(rna_de[g, ])
      ),
      data.frame(
        gene = g,
        sample_id = clinical_compare$Sample.ID,
        subtype = clinical_compare$subtype,
        layer = "Protein",
        value = as.numeric(protein_de[g, ])
      )
    )
  })
)

# Summarize observed and missing measurements within each group
discordant_qc <- aggregate(
  cbind(
    observed = as.integer(!is.na(discordant_values$value)),
    missing = as.integer(is.na(discordant_values$value))
  ),
  by = discordant_values[c("gene", "layer", "subtype")],
  FUN = sum
)

discordant_qc$missing_percent <- with(
  discordant_qc,
  100 * missing / (observed + missing)
)

discordant_qc


#plot the first three genes in our existing evidence ranking: GRAMD1C, ABCA13, and SHC4
example_genes <- head(discordant_genes$gene, 3)

example_values <- subset(
  discordant_values,
  gene %in% example_genes & !is.na(value)
)

example_values$gene <- factor(
  example_values$gene,
  levels = example_genes
)

example_values$layer <- factor(
  example_values$layer,
  levels = c("RNA", "Protein")
)

discordant_plot <- ggplot2::ggplot(
  example_values,
  ggplot2::aes(x = subtype, y = value, fill = subtype)
) +
  ggplot2::geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    alpha = 0.5
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.12,
      height = 0,
      seed = 42
    ),
    size = 1.3,
    alpha = 0.7
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(gene, layer),
    ncol = 2,
    scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(
    values = c("LumA" = "#0072B2", "Basal" = "#D55E00")
  ) +
  ggplot2::labs(
    title = "Candidate discordant subtype associations",
    subtitle = "Observed measurements; missing values excluded from display",
    x = NULL,
    y = "Normalized log2 measurement"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(legend.position = "none")

print(discordant_plot)


subset(discordant_qc, gene %in% example_genes)


# Calculate subtype effects using the same patients in both layers
paired_sensitivity <- do.call(
  rbind,
  lapply(example_genes, function(g) {
    
    rna_values <- as.numeric(rna_de[g, ])
    protein_values <- as.numeric(protein_de[g, ])
    
    # Keep patients with both measurements
    paired <- is.finite(rna_values) &
      is.finite(protein_values)
    
    luma <- paired & clinical_compare$subtype == "LumA"
    basal <- paired & clinical_compare$subtype == "Basal"
    
    data.frame(
      gene = g,
      paired_LumA = sum(luma),
      paired_Basal = sum(basal),
      
      rna_mean_effect =
        mean(rna_values[basal]) - mean(rna_values[luma]),
      
      protein_mean_effect =
        mean(protein_values[basal]) - mean(protein_values[luma]),
      
      rna_median_effect =
        median(rna_values[basal]) - median(rna_values[luma]),
      
      protein_median_effect =
        median(protein_values[basal]) - median(protein_values[luma])
    )
  })
)

# Attach the original model estimates for comparison
paired_sensitivity <- merge(
  paired_sensitivity,
  integrated_results[
    , c("gene", "rna_logFC", "protein_logFC")
  ],
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

paired_sensitivity


# Select patients with both clinical covariates available
complete_covariates <- with(
  clinical_compare,
  is.finite(age_years) & is.finite(tumor_cellularity)
)

clinical_sensitivity <- clinical_compare[
  complete_covariates, , drop = FALSE
]

table(clinical_sensitivity$subtype)

# Create two models using the same patients
design_unadjusted <- model.matrix(
  ~ subtype,
  data = clinical_sensitivity
)

design_adjusted <- model.matrix(
  ~ subtype + age_years + tumor_cellularity,
  data = clinical_sensitivity
)

# Confirm that all coefficients can be estimated
stopifnot(
  qr(design_unadjusted)$rank == ncol(design_unadjusted),
  qr(design_adjusted)$rank == ncol(design_adjusted)
)

# Subset expression matrices to these same patients
rna_sensitivity <- rna_de[
  , complete_covariates, drop = FALSE
]

protein_sensitivity <- protein_de[
  , complete_covariates, drop = FALSE
]

stopifnot(
  identical(
    colnames(rna_sensitivity),
    clinical_sensitivity$Sample.ID
  ),
  identical(
    colnames(protein_sensitivity),
    clinical_sensitivity$Sample.ID
  )
)

# Fit unadjusted models on the reduced patient set
rna_fit_cc <- limma::eBayes(
  limma::lmFit(rna_sensitivity, design_unadjusted),
  trend = FALSE
)

protein_fit_cc <- limma::eBayes(
  limma::lmFit(protein_sensitivity, design_unadjusted),
  trend = FALSE
)

# Fit adjusted models on the identical patient set
rna_fit_adjusted <- limma::eBayes(
  limma::lmFit(rna_sensitivity, design_adjusted),
  trend = FALSE
)

protein_fit_adjusted <- limma::eBayes(
  limma::lmFit(protein_sensitivity, design_adjusted),
  trend = FALSE
)


# Compare subtype coefficients before and after adjustment
adjustment_check <- data.frame(
  gene = example_genes,
  
  rna_unadjusted_cc =
    rna_fit_cc$coefficients[example_genes, "subtypeBasal"],
  
  rna_adjusted =
    rna_fit_adjusted$coefficients[example_genes, "subtypeBasal"],
  
  protein_unadjusted_cc =
    protein_fit_cc$coefficients[example_genes, "subtypeBasal"],
  
  protein_adjusted =
    protein_fit_adjusted$coefficients[example_genes, "subtypeBasal"]
)

adjustment_check


# Extract all adjusted-model results
rna_adjusted_results <- limma::topTable(
  rna_fit_adjusted,
  coef = "subtypeBasal",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

protein_adjusted_results <- limma::topTable(
  protein_fit_adjusted,
  coef = "subtypeBasal",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

# Attach adjusted-model FDR values to our comparison table
adjustment_check$rna_adjusted_fdr <-
  rna_adjusted_results[
    adjustment_check$gene, "adj.P.Val"
  ]

adjustment_check$protein_adjusted_fdr <-
  protein_adjusted_results[
    adjustment_check$gene, "adj.P.Val"
  ]

# Check whether opposite directions and both criteria persist
adjustment_check$passes_adjusted_criteria <- with(
  adjustment_check,
  !is.na(rna_adjusted_fdr) &
    !is.na(protein_adjusted_fdr) &
    rna_adjusted_fdr < fdr_threshold &
    protein_adjusted_fdr < fdr_threshold &
    abs(rna_adjusted) >= rna_effect_threshold &
    abs(protein_adjusted) >= protein_effect_threshold &
    sign(rna_adjusted) != sign(protein_adjusted)
)

# Display the key results
adjustment_check[
  , c(
    "gene",
    "rna_adjusted", "rna_adjusted_fdr",
    "protein_adjusted", "protein_adjusted_fdr",
    "passes_adjusted_criteria"
  )
]
