# Data access

This project reanalyzes published processed data; it does not process raw sequencing reads or mass-spectrometry files.

Source: [LinkedOmics CPTAC breast cancer prospective cohort](https://www.linkedomics.org/data_download/CPTAC-BRCA/).

Download these three files without changing their names and place them in a local `data/` folder at the repository root:

| Portal row | Filename | Measurements |
| --- | --- | --- |
| Clinical | `HS_CPTAC_BRCA_2018_CLI.tsi` | Clinical and molecular annotations |
| RNAseq (Gene level, Tumor) | `HS_CPTAC_BRCA_2018_RNA_GENE.cct` | Log2 FPKM normalized by gene median |
| Proteome (Gene level, TMT Log Ratio, Tumor) | `HS_CPTAC_BRCA_2018_Proteome_Ratio_Norm_gene_Median.cct` | Normalized TMT log2 ratios |

Filenames are case-sensitive on some systems. The matrices are tab-delimited text despite their `.cct` and `.tsi` extensions. The clinical file's first data row describes variable types; the script preserves it separately and removes it from patient records.

The analyzed downloads contain 122 samples, 23,121 RNA features, and 12,621 protein features. The script verifies sample identifiers, order, numeric expression columns, and gene identifiers. Source matrices are excluded from Git tracking. Download access and source terms remain governed by their original providers.

Associated publication: Krug et al. (2020), *Proteogenomic Landscape of Breast Cancer Tumorigenesis and Targeted Therapy*, Cell. [NCI study page](https://gdc.cancer.gov/about-data/publications/CPTAC-3_2020_1).

For reproducibility, record the date of your own download and preserve the original inputs. The `2018` in these filenames is part of the source naming convention, not the publication year.
