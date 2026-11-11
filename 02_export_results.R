# Run after 01_import_data.R in the same R session, from the project root.
required_objects <- c("rna_volcano", "protein_volcano", "effect_plot", "discordant_plot",
                      "rna_results", "protein_results", "integrated_results",
                      "discordant_genes", "paired_sensitivity", "adjustment_check")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects)) {
  stop("Run 01_import_data.R first. Missing objects: ", paste(missing_objects, collapse = ", "))
}
for (folder in c("results/figures", "results/tables", "docs")) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}
plots_to_save <- list(
  rna_volcano = rna_volcano,
  protein_volcano = protein_volcano,
  rna_protein_effect_agreement = effect_plot,
  discordant_gene_examples = discordant_plot
)
for (plot_name in names(plots_to_save)) {
  ggplot2::ggsave(
    filename = file.path("results/figures", paste0(plot_name, ".png")),
    plot = plots_to_save[[plot_name]], width = 10, height = 8,
    units = "in", dpi = 300, bg = "white"
  )
}
tables_to_save <- list(
  rna_differential_results = rna_results,
  protein_differential_results = protein_results,
  integrated_rna_protein_results = integrated_results,
  discordant_genes = discordant_genes,
  paired_sample_sensitivity = paired_sensitivity,
  clinical_adjustment_check = adjustment_check,
  missingness_summary = missing_summary,
  filtering_summary = filter_summary,
  discordant_gene_missingness = discordant_qc
)
for (table_name in names(tables_to_save)) {
  write.csv(tables_to_save[[table_name]],
            file.path("results/tables", paste0(table_name, ".csv")),
            row.names = FALSE)
}
writeLines(capture.output(sessionInfo()), "docs/sessionInfo.txt")
message("Exported figures, tables, and session information.")
