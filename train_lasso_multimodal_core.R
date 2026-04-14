library(caret)
library(dplyr)
library(glmnet)

# -----------------------------
# Config: update these paths
# -----------------------------
cfg <- list(
  methylation_rds = "data/Mvalue_trans.rds",
  cytokine_csv = "data/cytokine_rank_normalized.csv",
  olink_csv = "data/olink_rank_normalized.csv",
  metabolite_csv = "data/metabolite_rank_normalized.csv",
  cellcounts_csv = "data/cellcounts_name_replace.csv",
  microbiome_csv = "data/microbiome_rank_normalized.csv",
  basic_phenos_csv = "data/Age_group_basicPhenos.csv",
  split_info_dir = "data/split_info",
  n_iterations = 100,
  output_rds = "results/elasticnet_multi_omics_results.rds",
  output_metrics_csv = "results/elasticnet_metrics_summary.csv"
)

read_inputs <- function(cfg) {
  list(
    Mvalue = readRDS(cfg$methylation_rds),
    cytokine = read.csv(cfg$cytokine_csv, row.names = 1, check.names = FALSE),
    olink = read.csv(cfg$olink_csv, row.names = 1, check.names = FALSE),
    metabolite = read.csv(cfg$metabolite_csv, row.names = 1, check.names = FALSE),
    cellcounts = read.csv(cfg$cellcounts_csv, row.names = 1, check.names = FALSE),
    microbiome = read.csv(cfg$microbiome_csv, row.names = 1, check.names = FALSE),
    basicPhenos = read.csv(cfg$basic_phenos_csv, row.names = 1, check.names = FALSE)
  )
}

build_combined_age <- function(dat) {
  common_samples <- Reduce(intersect, list(
    rownames(dat$cytokine),
    rownames(dat$olink),
    rownames(dat$metabolite),
    rownames(dat$cellcounts),
    rownames(dat$microbiome),
    rownames(dat$Mvalue),
    dat$basicPhenos$ID_500fg
  ))
  common_samples <- sort(common_samples)

  cytokine_common <- dat$cytokine[common_samples, , drop = FALSE]
  olink_common <- dat$olink[common_samples, , drop = FALSE]
  metabolite_common <- dat$metabolite[common_samples, , drop = FALSE]
  cellcounts_common <- dat$cellcounts[common_samples, , drop = FALSE]
  microbiome_common <- dat$microbiome[common_samples, , drop = FALSE]
  Mvalue_common <- dat$Mvalue[common_samples, , drop = FALSE]

  basicPhenos_common <- dat$basicPhenos[dat$basicPhenos$ID_500fg %in% common_samples, , drop = FALSE]
  basicPhenos_common <- basicPhenos_common[order(match(basicPhenos_common$ID_500fg, common_samples)), ]

  combined_table <- cbind(
    cytokine_common, olink_common, metabolite_common,
    cellcounts_common, microbiome_common, Mvalue_common
  )
  rownames(combined_table) <- common_samples
  combined_age <- cbind(combined_table, Age = basicPhenos_common$Age)

  list(
    combined_age = combined_age,
    feature_cols = list(
      cytokine = colnames(dat$cytokine),
      olink = colnames(dat$olink),
      metabolite = colnames(dat$metabolite),
      cellcounts = colnames(dat$cellcounts),
      microbiome = colnames(dat$microbiome)
    )
  )
}

train_elasticnet_once <- function(training_data, validation_data) {
  net_grid <- expand.grid(.alpha = seq(0, 1, by = 0.1), .lambda = seq(0, 1, by = 0.01))
  ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)

  fit <- train(
    Age ~ ., data = training_data, method = "glmnet", metric = "RMSE",
    tuneGrid = net_grid, trControl = ctrl, preProcess = "scale"
  )

  train_pred <- predict(fit, newdata = training_data)
  val_pred <- predict(fit, newdata = validation_data)

  model_coefs <- as.matrix(coef(fit$finalModel, s = fit$bestTune$.lambda))
  selected_features <- rownames(model_coefs)[model_coefs != 0]
  selected_features <- setdiff(selected_features, "(Intercept)")

  list(
    train_rmse = sqrt(mean((training_data$Age - train_pred)^2)),
    val_rmse = sqrt(mean((validation_data$Age - val_pred)^2)),
    train_r2 = cor(training_data$Age, train_pred)^2,
    val_r2 = cor(validation_data$Age, val_pred)^2,
    best_alpha = fit$bestTune$.alpha,
    best_lambda = fit$bestTune$.lambda,
    selected_features = selected_features
  )
}

run_experiment <- function(cfg, combined_age, feature_cols) {
  results <- vector("list", cfg$n_iterations)

  for (i in seq_len(cfg$n_iterations)) {
    message("Iteration: ", i)
    split_path <- file.path(cfg$split_info_dir, paste0("iteration_", i, "_fdr_info.rds"))

    results[[i]] <- tryCatch({
      info <- readRDS(split_path)
      split_indices <- as.integer(unlist(info$split_indices1))

      mvalue_features <- info$selected_features
      if (length(mvalue_features) > 20000 && !is.null(info$fdr)) {
        names(info$fdr) <- info$selected_features
        mvalue_features <- names(sort(info$fdr, decreasing = FALSE))[1:10000]
      }

      all_features <- unique(c(
        feature_cols$cytokine,
        feature_cols$olink,
        feature_cols$metabolite,
        feature_cols$cellcounts,
        feature_cols$microbiome,
        mvalue_features
      ))
      all_features <- setdiff(all_features, "Age")
      available_features <- intersect(all_features, colnames(combined_age))

      if (length(available_features) == 0) {
        stop("No available features in combined_age.")
      }

      training <- combined_age[split_indices, c(available_features, "Age"), drop = FALSE]
      validation <- combined_age[-split_indices, c(available_features, "Age"), drop = FALSE]
      train_elasticnet_once(training, validation)
    }, error = function(e) {
      message("Iteration ", i, " failed: ", e$message)
      NULL
    })

    gc(verbose = FALSE)
  }

  results[!sapply(results, is.null)]
}

summarize_metrics <- function(results) {
  data.frame(
    metric = c("train_rmse", "val_rmse", "train_r2", "val_r2", "n_selected_features"),
    mean = c(
      mean(sapply(results, `[[`, "train_rmse")),
      mean(sapply(results, `[[`, "val_rmse")),
      mean(sapply(results, `[[`, "train_r2")),
      mean(sapply(results, `[[`, "val_r2")),
      mean(sapply(results, function(x) length(x$selected_features)))
    ),
    sd = c(
      sd(sapply(results, `[[`, "train_rmse")),
      sd(sapply(results, `[[`, "val_rmse")),
      sd(sapply(results, `[[`, "train_r2")),
      sd(sapply(results, `[[`, "val_r2")),
      sd(sapply(results, function(x) length(x$selected_features)))
    )
  )
}

main <- function(cfg) {
  dir.create(dirname(cfg$output_rds), recursive = TRUE, showWarnings = FALSE)

  dat <- read_inputs(cfg)
  prepared <- build_combined_age(dat)

  results <- run_experiment(cfg, prepared$combined_age, prepared$feature_cols)
  saveRDS(results, cfg$output_rds)

  metrics <- summarize_metrics(results)
  write.csv(metrics, cfg$output_metrics_csv, row.names = FALSE)
  print(metrics)
}

main(cfg)
