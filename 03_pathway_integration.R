#==========================================================
#Pathway level integration
#==========================================================

## Install packages if needed
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", ask = FALSE, update = FALSE)
}

# Retrieve human Hallmark gene sets
hallmark_data <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)

# Convert the table into a named list of gene sets
hallmark_pathways <- lapply(
  split(hallmark_data$gene_symbol, hallmark_data$gs_name),
  unique
)

# Record the gene-set database version
unique(hallmark_data$db_version)
# gene-set database versions" = 2026.1.Hs"

# Check the number and sizes of gene sets
length(hallmark_pathways)
summary(lengths(hallmark_pathways))






# rank all tested genes from primary unadjusted analysis

make_gene_ranks <- function(results) {
  
  keep <- is.finite(results$t) &
    !is.na(results$gene) &
    nzchar(results$gene)
  
  rank_data <- results[keep, c("gene", "t")]
  
  stopifnot(anyDuplicated(rank_data$gene) == 0L)
  
  ranks <- rank_data$t
  names(ranks) <- rank_data$gene
  
  sort(ranks, decreasing = TRUE)
}

rna_ranks <- make_gene_ranks(rna_results)
protein_ranks <- make_gene_ranks(protein_results)

# Inspect the rankings
head(rna_ranks)
head(protein_ranks)

# Check how many ranked genes map to Hallmark gene sets
hallmark_genes <- unique(unlist(
  hallmark_pathways,
  use.names = FALSE
))

c(
  RNA_ranked = length(rna_ranks),
  RNA_in_Hallmark = sum(names(rna_ranks) %in% hallmark_genes),
  Protein_ranked = length(protein_ranks),
  Protein_in_Hallmark = sum(names(protein_ranks) %in% hallmark_genes)
)





# Confirm valid, uniquely named ranking statistics
stopifnot(
  all(is.finite(rna_ranks)),
  all(is.finite(protein_ranks)),
  anyDuplicated(names(rna_ranks)) == 0L,
  anyDuplicated(names(protein_ranks)) == 0L
)

# Run preranked enrichment for RNA
set.seed(42)

rna_gsea <- fgsea::fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = rna_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0,
  nproc = 1
)

# Run the same analysis for protein
set.seed(42)

protein_gsea <- fgsea::fgseaMultilevel(
  pathways = hallmark_pathways,
  stats = protein_ranks,
  minSize = 15,
  maxSize = 500,
  eps = 0,
  nproc = 1
)

# Convert results to ordinary data frames
rna_gsea <- as.data.frame(rna_gsea)
protein_gsea <- as.data.frame(protein_gsea)

# Sort by adjusted p-value
rna_gsea <- rna_gsea[
  order(rna_gsea$padj), , drop = FALSE
]

protein_gsea <- protein_gsea[
  order(protein_gsea$padj), , drop = FALSE
]

# Inspect the ten highest-ranked pathways
columns_to_show <- c("pathway", "size", "NES", "pval", "padj")

head(rna_gsea[, columns_to_show], 10)
head(protein_gsea[, columns_to_show], 10)



pathway_summary <- data.frame(
  layer = c("RNA", "Protein"),
  
  pathways_tested = c(
    nrow(rna_gsea),
    nrow(protein_gsea)
  ),
  
  enriched_toward_Basal = c(
    sum(rna_gsea$padj < 0.05 & rna_gsea$NES > 0, na.rm = TRUE),
    sum(protein_gsea$padj < 0.05 & protein_gsea$NES > 0, na.rm = TRUE)
  ),
  
  enriched_toward_LumA = c(
    sum(rna_gsea$padj < 0.05 & rna_gsea$NES < 0, na.rm = TRUE),
    sum(protein_gsea$padj < 0.05 & protein_gsea$NES < 0, na.rm = TRUE)
  ),
  
  missing_adjusted_pvalues = c(
    sum(is.na(rna_gsea$padj)),
    sum(is.na(protein_gsea$padj))
  )
)

pathway_summary




# Select pathway statistics from each layer
rna_pathways <- rna_gsea[
  , c("pathway", "size", "NES", "padj")
]

protein_pathways <- protein_gsea[
  , c("pathway", "size", "NES", "padj")
]

# Give the columns layer-specific names
colnames(rna_pathways) <- c(
  "pathway", "rna_size", "rna_NES", "rna_fdr"
)

colnames(protein_pathways) <- c(
  "pathway", "protein_size", "protein_NES", "protein_fdr"
)

# Retain pathways tested in either layer
pathway_comparison <- merge(
  rna_pathways,
  protein_pathways,
  by = "pathway",
  all = TRUE
)

# Identify pathways with usable results in both layers
valid <- with(
  pathway_comparison,
  is.finite(rna_NES) &
    is.finite(protein_NES) &
    !is.na(rna_fdr) &
    !is.na(protein_fdr)
)

rna_significant <- with(
  pathway_comparison,
  valid & rna_fdr < 0.05
)

protein_significant <- with(
  pathway_comparison,
  valid & protein_fdr < 0.05
)

same_direction <- with(
  pathway_comparison,
  sign(rna_NES) == sign(protein_NES)
)

# Assign evidence categories
pathway_comparison$category <- "Not evaluated in both"

pathway_comparison$category[valid] <- "Below criteria in both"

pathway_comparison$category[
  rna_significant & !protein_significant
] <- "RNA criteria only"

pathway_comparison$category[
  !rna_significant & protein_significant
] <- "Protein criteria only"

pathway_comparison$category[
  which(rna_significant & protein_significant & same_direction)
] <- "Concordant in both"

pathway_comparison$category[
  which(rna_significant & protein_significant & !same_direction)
] <- "Opposite in both"

table(pathway_comparison$category)

# Inspect pathways meeting criteria in both layers
both_layers <- subset(
  pathway_comparison,
  category %in% c("Concordant in both", "Opposite in both")
)

# Rank using the less significant of the two adjusted p-values
both_layers <- both_layers[
  order(pmax(both_layers$rna_fdr, both_layers$protein_fdr)),
  ,
  drop = FALSE
]

both_layers[
  , c(
    "pathway",
    "rna_NES", "rna_fdr",
    "protein_NES", "protein_fdr",
    "category"
  )
]





# Order pathways by average NES across the two layers
plot_order <- pathway_comparison$pathway[
  order(
    rowMeans(
      pathway_comparison[, c("rna_NES", "protein_NES")],
      na.rm = TRUE
    )
  )
]

# Stack RNA and protein results into one plotting table
pathway_plot_data <- rbind(
  data.frame(
    pathway = pathway_comparison$pathway,
    layer = "RNA",
    NES = pathway_comparison$rna_NES,
    fdr = pathway_comparison$rna_fdr
  ),
  data.frame(
    pathway = pathway_comparison$pathway,
    layer = "Protein",
    NES = pathway_comparison$protein_NES,
    fdr = pathway_comparison$protein_fdr
  )
)

pathway_plot_data$pathway <- factor(
  pathway_plot_data$pathway,
  levels = plot_order
)

pathway_plot_data$layer <- factor(
  pathway_plot_data$layer,
  levels = c("RNA", "Protein")
)

# Mark evidence using shape; keep color for enrichment direction
pathway_plot_data$evidence <- ifelse(
  is.na(pathway_plot_data$fdr),
  "Not evaluated",
  ifelse(pathway_plot_data$fdr < 0.05, "FDR < 0.05", "FDR >= 0.05")
)

# Symmetric NES color scale
nes_limit <- max(abs(pathway_plot_data$NES), na.rm = TRUE)

pathway_plot <- ggplot2::ggplot(
  pathway_plot_data,
  ggplot2::aes(x = layer, y = pathway)
) +
  ggplot2::geom_point(
    ggplot2::aes(color = NES, shape = evidence),
    size = 3
  ) +
  ggplot2::scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-nes_limit, nes_limit)
  ) +
  ggplot2::scale_shape_manual(
    values = c(
      "FDR < 0.05" = 16,
      "FDR >= 0.05" = 1,
      "Not evaluated" = 4
    )
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) {
      gsub("_", " ", sub("^HALLMARK_", "", x))
    }
  ) +
  ggplot2::labs(
    title = "Hallmark enrichment across RNA and protein",
    subtitle = "Positive NES: Basal-associated | Negative NES: LumA-associated",
    x = NULL,
    y = NULL,
    color = "NES",
    shape = NULL
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    legend.position = "right",
    axis.ticks.y = ggplot2::element_blank()
  )

print(pathway_plot)



ggplot2::ggsave(
  "results/figures/hallmark_rna_protein_comparison.png",
  plot = pathway_plot,
  width = 10,
  height = 14,
  units = "in",
  dpi = 300,
  bg = "white"
)






# Find genes available in both rankings
shared_rank_genes <- intersect(
  names(rna_ranks),
  names(protein_ranks)
)

length(shared_rank_genes)

# Keep identical gene membership, then rank each layer independently
rna_ranks_shared <- sort(
  rna_ranks[shared_rank_genes],
  decreasing = TRUE
)

protein_ranks_shared <- sort(
  protein_ranks[shared_rank_genes],
  decreasing = TRUE
)

stopifnot(
  setequal(
    names(rna_ranks_shared),
    names(protein_ranks_shared)
  )
)

# Repeat enrichment with the same settings
set.seed(42)

rna_gsea_shared <- as.data.frame(
  fgsea::fgseaMultilevel(
    pathways = hallmark_pathways,
    stats = rna_ranks_shared,
    minSize = 15,
    maxSize = 500,
    eps = 0,
    nproc = 1
  )
)

set.seed(42)

protein_gsea_shared <- as.data.frame(
  fgsea::fgseaMultilevel(
    pathways = hallmark_pathways,
    stats = protein_ranks_shared,
    minSize = 15,
    maxSize = 500,
    eps = 0,
    nproc = 1
  )
)




compare_coverage <- function(original, shared) {
  
  comparison <- merge(
    original[, c("pathway", "NES", "padj", "size")],
    shared[, c("pathway", "NES", "padj", "size")],
    by = "pathway",
    suffixes = c("_original", "_shared"),
    all = TRUE
  )
  
  comparison$NES_change <-
    comparison$NES_shared - comparison$NES_original
  
  comparison$direction_same <- with(
    comparison,
    sign(NES_original) == sign(NES_shared)
  )
  
  comparison$significant_original <-
    comparison$padj_original < 0.05
  
  comparison$significant_shared <-
    comparison$padj_shared < 0.05
  
  comparison
}

rna_coverage_check <- compare_coverage(
  rna_gsea, rna_gsea_shared
)

protein_coverage_check <- compare_coverage(
  protein_gsea, protein_gsea_shared
)

# Count changes in significance status
with(
  rna_coverage_check,
  table(
    Original = significant_original,
    Shared = significant_shared,
    useNA = "ifany"
  )
)

with(
  protein_coverage_check,
  table(
    Original = significant_original,
    Shared = significant_shared,
    useNA = "ifany"
  )
)

# Inspect the largest changes in NES
head(
  rna_coverage_check[
    order(-abs(rna_coverage_check$NES_change)),
  ],
  10
)

head(
  protein_coverage_check[
    order(-abs(protein_coverage_check$NES_change)),
  ],
  10
)





# Inspect pathways with missing shared-analysis results
subset(
  rna_coverage_check,
  is.na(NES_shared) | is.na(padj_shared)
)

subset(
  protein_coverage_check,
  is.na(NES_shared) | is.na(padj_shared)
)

# List pathways whose significance status changed
changed_columns <- c(
  "pathway",
  "NES_original", "padj_original", "size_original",
  "NES_shared", "padj_shared", "size_shared"
)

subset(
  rna_coverage_check,
  !is.na(significant_original) &
    !is.na(significant_shared) &
    significant_original != significant_shared,
  select = changed_columns
)

subset(
  protein_coverage_check,
  !is.na(significant_original) &
    !is.na(significant_shared) &
    significant_original != significant_shared,
  select = changed_columns
)



# Count shared ranked genes in every Hallmark set
shared_pathway_sizes <- data.frame(
  pathway = names(hallmark_pathways),
  shared_genes = vapply(
    hallmark_pathways,
    function(genes) sum(genes %in% shared_rank_genes),
    integer(1)
  )
)

subset(shared_pathway_sizes, shared_genes < 15)





# Join shared-gene enrichment results across layers
shared_pathway_comparison <- merge(
  rna_gsea_shared[, c("pathway", "NES", "padj")],
  protein_gsea_shared[, c("pathway", "NES", "padj")],
  by = "pathway",
  suffixes = c("_rna", "_protein"),
  all = TRUE
)

# Identify concordant enrichment supported in both layers
shared_pathway_comparison$concordant_significant <- with(
  shared_pathway_comparison,
  is.finite(NES_rna) &
    is.finite(NES_protein) &
    !is.na(padj_rna) &
    !is.na(padj_protein) &
    padj_rna < 0.05 &
    padj_protein < 0.05 &
    sign(NES_rna) == sign(NES_protein)
)

# Extract the original 21 concordant pathways
original_concordant <- subset(
  pathway_comparison,
  category == "Concordant in both"
)$pathway

# Compare their status after restricting gene coverage
concordance_sensitivity <- shared_pathway_comparison[
  match(original_concordant, shared_pathway_comparison$pathway),
  ,
  drop = FALSE
]

stopifnot(!anyNA(concordance_sensitivity$pathway))

table(concordance_sensitivity$concordant_significant)

# Show any original concordant pathways that no longer meet both criteria
subset(
  concordance_sensitivity,
  !concordant_significant
)

# Count all concordant pathways in the shared-gene analysis
sum(shared_pathway_comparison$concordant_significant)





# Identify newly concordant pathways under shared-gene coverage
subset(
  shared_pathway_comparison,
  concordant_significant &
    !(pathway %in% original_concordant)
)




# dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Convert fgsea list columns to text so they can be exported to CSV
export_fgsea <- function(x, file) {
  out <- as.data.frame(x)
  
  list_columns <- vapply(out, is.list, logical(1))
  
  out[list_columns] <- lapply(
    out[list_columns],
    function(column) {
      vapply(column, paste, collapse = ";", FUN.VALUE = character(1))
    }
  )
  
  write.csv(out, file, row.names = FALSE)
}

export_fgsea(
  rna_gsea,
  "results/tables/rna_hallmark_gsea.csv"
)

export_fgsea(
  protein_gsea,
  "results/tables/protein_hallmark_gsea.csv"
)

export_fgsea(
  rna_gsea_shared,
  "results/tables/rna_hallmark_gsea_shared_genes.csv"
)

export_fgsea(
  protein_gsea_shared,
  "results/tables/protein_hallmark_gsea_shared_genes.csv"
)

tables_to_export <- list(
  hallmark_pathway_comparison = pathway_comparison,
  rna_pathway_coverage_sensitivity = rna_coverage_check,
  protein_pathway_coverage_sensitivity = protein_coverage_check,
  shared_gene_pathway_comparison = shared_pathway_comparison,
  original_concordant_pathway_sensitivity = concordance_sensitivity,
  hallmark_shared_gene_sizes = shared_pathway_sizes
)

for (table_name in names(tables_to_export)) {
  write.csv(
    tables_to_export[[table_name]],
    file.path("results/tables", paste0(table_name, ".csv")),
    row.names = FALSE
  )
}

# Refresh the package and R version record
writeLines(
  capture.output(sessionInfo()),
  "docs/sessionInfo.txt"
)










