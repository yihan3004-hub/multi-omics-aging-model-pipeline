# Multi-Omics Age Prediction (Elastic Net)

This repository contains a reproducible core pipeline for multi-omics age prediction using Elastic Net regression in R.

## Main Script

- `train_lasso_multimodal_core.R`  
  (The filename is kept for compatibility, but the model is now Elastic Net.)

## Expected Project Structure

```text
project/
  train_lasso_multimodal_core.R
  README.md
  .gitignore
  LICENSE
  sessionInfo.txt
  data/
    Mvalue_trans.rds
    cytokine_rank_normalized.csv
    olink_rank_normalized.csv
    metabolite_rank_normalized.csv
    cellcounts_name_replace.csv
    microbiome_rank_normalized.csv
    Age_group_basicPhenos.csv
    split_info/
      iteration_1_fdr_info.rds
      iteration_2_fdr_info.rds
      ...
  results/
```

## Data Requirements

1. All omics matrices must use sample IDs as row names.
2. `Age_group_basicPhenos.csv` must contain:
   - `ID_500fg` (sample ID)
   - `Age` (target value)
3. Each split file in `data/split_info/` must contain:
   - `split_indices1`
   - `selected_features`
   - optional `fdr`

## Install Dependencies

In R:

```r
install.packages(c("caret", "dplyr", "glmnet"))
```

## Run

From the project directory:

```bash
Rscript train_lasso_multimodal_core.R
```

## Outputs

- `results/elasticnet_multi_omics_results.rds`  
  Per-iteration model results.
- `results/elasticnet_metrics_summary.csv`  
  Mean/SD summary for RMSE, R2, and selected feature count.

## Reproducibility

Before publishing results, refresh environment information:

```r
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
```

