# CPTAC Breast Cancer Multi-omics Analysis

An ongoing R portfolio project by **Enoch Larrey**, integrating published breast cancer transcriptomics, proteomics, and clinical annotations to examine agreement and discordance between molecular layers.

## Research question

Which molecular features distinguish Basal and Luminal A (LumA) breast tumors, and how consistently do RNA and protein measurements reflect those subtype associations?

## Data and scope

Processed data are from the CPTAC prospective breast cancer cohort associated with Krug et al. (2020), *Proteogenomic Landscape of Breast Cancer Tumorigenesis and Targeted Therapy*, Cell.

- [NCI publication and study page](https://gdc.cancer.gov/about-data/publications/CPTAC-3_2020_1)
- [LinkedOmics data downloads](https://www.linkedomics.org/data_download/CPTAC-BRCA/)
- [Exact input files and access instructions](docs/data_access.md)

The cohort contains 122 tumors. The primary differential comparison uses 57 LumA and 29 Basal tumors. This is an exploratory reanalysis of processed matrices, rather than a raw-data processing pipeline or a clinical prediction tool.

## Analysis completed

1. Import tab-delimited RNA, protein, and clinical data; preserve clinical variable annotations.
2. Verify sample identities and align all layers.
3. Inspect missingness and sample distributions; retain genes observed in at least 80% of the full cohort.
4. Explore each layer with centered, unscaled PCA using complete, variable genes.
5. Require at least 80% gene coverage within each comparison group, then fit separate limma models (`~ subtype`, LumA reference; `eBayes(trend = FALSE)`). Missing measurements remain NA.
6. Prioritize genes using BH-adjusted p-value < 0.05 and absolute log2 effect >= 1 in each layer. BH correction is performed separately by layer.
7. Compare subtype effects across shared tested genes and inspect candidate opposite-direction associations.
8. Check three ranked examples using matched observed patients, medians, and age/cellularity-adjusted models. Adjusted and unadjusted sensitivity models use the same 70 patients (51 LumA, 19 Basal).

9. Run preranked Hallmark enrichment using moderated t statistics and compare RNA/protein enrichment directions.
10. Repeat enrichment using the 8,133 shared ranked genes to assess coverage sensitivity.

## Current findings

Results below were observed during the interactive analysis; a fresh-session reproduction check is still pending.

| Metric | RNA | Protein |
| --- | ---: | ---: |
| Features tested in the primary comparison | 15,743 | 8,777 |
| Higher in LumA, meeting both criteria | 1,830 | 1,375 |
| Higher in Basal, meeting both criteria | 1,493 | 739 |

Across 8,133 shared tested genes, subtype effects have Pearson correlation **0.739** and Spearman correlation **0.697**. Among genes meeting both criteria in both layers, 548 are higher in LumA in both, 320 are higher in Basal in both, and 12 have opposite directions.

The three highest-ranked discordant examples—GRAMD1C, ABCA13, and SHC4—retain opposite directions in paired-patient descriptive checks and meet the chosen criteria in age/cellularity-adjusted models. These are candidate associations, not independently validated mechanisms.

## Pathway checkpoint

Hallmark enrichment used MSigDB **2026.1.Hs** and `fgseaMultilevel` with `minSize = 15`, `maxSize = 500`, `eps = 0`, `nproc = 1`, and seed 42 for each run. All tested genes are retained in each ranking; rankings are not restricted to significant genes or Hallmark members. FDR < 0.05 is applied separately to each enrichment run.

| Primary enrichment result | RNA | Protein |
| --- | ---: | ---: |
| Hallmark pathways tested | 50 | 50 |
| Basal-associated (positive NES), FDR < 0.05 | 21 | 15 |
| LumA-associated (negative NES), FDR < 0.05 | 6 | 11 |

There are **21 pathways significant in both layers with matching directions**, no significant opposite-direction pathways, six meeting RNA criteria only, five meeting protein criteria only, and 18 below criteria in both. Basal-associated examples include E2F targets, G2M checkpoint, and MYC targets. Early and late estrogen-response gene sets are LumA-associated in both layers.

Restricting both rankings to shared genes retains **20 of the original 21** concordant pathways. TNFA signaling via NFκB keeps a positive NES in both layers but its RNA FDR rises to 0.076 (protein FDR 0.009). The shared-gene analysis contains 21 concordant pathways overall, with changed membership. PANCREAS BETA CELLS is excluded from the shared-gene runs because only 11 genes remain, below the minimum size of 15.

The pathway script completed on the author's machine with RNA ranking-tie warnings. Duplicate nonzero statistics account for 603/15,743 RNA ranking entries (3.83%) and 55/8,133 shared RNA entries (0.68%); neither protein ranking has ties. Original ranking values and ordering are retained. No jitter or alternative tie-order analysis is applied; sensitivity to tie ordering has not been assessed.

This checkpoint includes the uploaded pathway script and documented results from the author's R outputs. The newly generated pathway figure, CSV exports, and refreshed environment record remain to be synchronized from the author's local project. The full workflow has not been independently rerun here.

## Run the project

1. Clone or download this repository and open an R session with its root as the working directory. An RStudio Project can be created in this existing directory.
2. Download the three inputs described in [data access](docs/data_access.md).
3. Install dependencies once:

```r
install.packages(c("ggplot2", "ggrepel", "msigdbr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("limma", "fgsea"), ask = FALSE, update = FALSE)
```

4. In a fresh R session, run these in order:

```r
source("01_import_data.R", encoding = "UTF-8")
source("02_export_results.R", encoding = "UTF-8")
source("03_pathway_integration.R", encoding = "UTF-8")
```

The first script creates analysis objects and plots. The second writes PNG figures, CSV tables, and `docs/sessionInfo.txt`, replacing prior exports with the same names. QC, PCA, and variance diagnostics are currently drawn during analysis but are not exported by the export script. The pathway script requires the RNA/protein result objects from the first script and writes its own figure and tables. MSigDB may change on later downloads; verify the retrieved database version against this checkpoint. Dataset files are not included in Git.

## Repository contents

| Path | Purpose |
| --- | --- |
| `01_import_data.R` | Current sequential analysis, from import through sensitivity checks |
| `02_export_results.R` | Export current results and record the R environment |
| `03_pathway_integration.R` | Hallmark enrichment, cross-layer comparison, shared-gene coverage checks, and exports |
| `docs/data_access.md` | Source links and exact input filenames |
| `results/figures/` | Exported analysis figures |
| `results/tables/` | Exported differential, integration, and sensitivity tables |

## Interpretation and limitations

- Positive effects mean Basal minus LumA on each layer's supplied log2 scale. These normalized measures should not be treated as interchangeable absolute abundances.
- PAM50 is expression-derived; RNA separation by PAM50 is not independent validation.
- Missingness may be nonrandom. ABCA13 RNA is missing in 9 of 57 LumA tumors and none of the 29 Basal tumors.
- Initial PCA uses complete genes and may underrepresent features with subtype-specific missingness.
- The full-cohort coverage filter can omit subtype-restricted features. Filtering sensitivity checks remain planned.
- Age/cellularity adjustment uses a smaller, selected clinical subset and does not address all potential confounding or batch effects. Cellularity adjustment changes the biological estimand.
- Meeting criteria in one layer only does not prove that effects differ between layers. The effect-size filter is not a formal test of a minimum true effect.
- Opposite associations do not establish post-transcriptional regulation. Protein quantification, isoforms, tissue composition, and other explanations require further investigation.
- Pathway NES indicates enrichment toward a subtype, not direct evidence of pathway activation. RNA and protein gene coverage differs, and shared-gene restriction changes the ranked background and tested gene sets. Separate FDR thresholds do not constitute a formal test of cross-layer differences.
- RNA ranking ties are retained and may affect enrichment or leading-edge membership; tie-order robustness is untested.
- Effect correlations compare subtype coefficients across genes, not RNA–protein correlations across patients within a gene.

## Planned extensions

- Coverage and filtering sensitivity analyses.
- Within-gene RNA–protein correlations, accounting for subtype structure.
- Unsupervised multi-omics integration and a reproducible report.

Later analyses will be added through commits in this repository. Pathway analysis is included in this checkpoint; its local generated exports are awaiting synchronization.
