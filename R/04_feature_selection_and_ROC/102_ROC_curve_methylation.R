# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data"
root.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan"
base.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation"
####
####  Load libraries
####
library(GSVA)
library(pROC)
library(randomForest)
library(randomForestSRC)
library(ggplot2)

library(foreach)
library(doParallel)
library(progressr)

library(caret)
library(ranger)
library(reshape2)

library(glmnet)
library(stabs)

base_seed <- 22
set.seed(base_seed)

# methylation.data <- read.csv(file = file.path(base.dir, "002 Differential Methylation Analysis DMA/Correct for PC1-3",
#                                               "sig_AsthmaVsHealthy_ATLANTIS_Mvalues_dasen.csv"),
#                              row.names = 1)
#
# pheno.data <- read.csv(file = file.path(base.dir, "100 Characteristics tables",
#                                         "ATLANTIS_linkingtable_with_group_Th.csv"),
#                        row.names = 1)
# pheno.eQTM <- read.csv(file = file.path(base.dir, "100 Characteristics tables",
#                                         "pheno_eQTM_samples.csv"),
#                        row.names = 1)
# DMA.results <- read.csv(file = file.path(base.dir, "002 Differential Methylation Analysis DMA/Correct for PC1-3",
#                                          "Tt2_significant_CpG_asthma_vs_control_BH.csv"),
#                         row.names = 1)

# Load data files
load(file = file.path(base.dir, "06_DMA/initialData.RData")) # Loads 'methylationData' object
pheno.data <- read.csv(
  file.path(base.dir, "06_DMA/ATLANTIS_linkingtable_with_group_Th.csv"),
  row.names = 1
)
pheno.eQTM <- read.csv(
  file.path(root.dir, "ATLANTIS_eQTM/01_data_prep/pheno_eQTM_samples.csv"),
  row.names = 1
)
DMA.results <- read.csv(
  file.path(base.dir, "06_DMA/Tt2_significant_CpG_asthma_vs_control_BH.csv"),
  row.names = 1
)

# Prepare methylation data
methylation.data <- methylationData
rm(methylationData)
rownames(methylation.data) <- methylation.data[, 1]
methylation.data <- methylation.data[, -1]
colnames(methylation.data) <- gsub("X", "", colnames(methylation.data))

# Filter pheno data and methylation data for relevant samples
# If no asthma status has been called, omit
# This is only in the linking table from Chiesi, not clinical data
pheno.data.filt <- pheno.data[!is.na(pheno.data$ASTHEA), ]

# Filter methylation data for relevant samples and transpose
methylation.data.filt <- methylation.data[, pheno.eQTM$meth_file_id]
# methylation.data.filt <- methylation.data[, colnames(methylation.data) %in% pheno.data.filt$meth_file_id]

methylation.data.filt <- as.data.frame(t(methylation.data.filt))

# Add asthma status to the filtered methylation data
methylation.data.filt$asthma <- pheno.data.filt[
  match(rownames(methylation.data.filt), pheno.eQTM$meth_file_id),
  "ASTHEA"
]
methylation.data.filt$asthma <- gsub(
  pattern = "A",
  replacement = 1,
  x = methylation.data.filt$asthma
)
methylation.data.filt$asthma <- gsub(
  pattern = "H",
  replacement = 0,
  x = methylation.data.filt$asthma
)
methylation.data.filt$asthma <- as.factor(methylation.data.filt$asthma)

# Prepare data matrix for GSVA and feature list
# data_matrix <- as.matrix(t(methylation.data.filt[, -ncol(methylation.data.filt)]))
all_features <- rownames(DMA.results)
# all_features <- sample(rownames(data_matrix), nrow(data_matrix)-1)
cat("Memory usage before RF: ", sum(gc()[, 2]), " MB\n")


cat("Performing Random Forest on 235 methylation sites\n")
# Choosing to random sample 70% of the total for training - 30% internally test
rf_data <- methylation.data.filt[, c(all_features, "asthma")]

################################################################################
#### Testing how many trees is the optimal tree count for random forest testing
####  Evaluated through the OOB (Out-of-Bag) method
# train_indices <- sample(seq_len(nrow(rf_data)), size = 0.7 * nrow(rf_data))
# train_data <- rf_data[train_indices, ]
# test_data <- rf_data[-train_indices, ]
#
#
# # Evaluate OOB Error for Different ntree Values
# oob_errors <- sapply(seq(100, 1000, by = 100), function(ntree) {
#   rf_model <- randomForest(asthma ~ ., data = train_data, ntree = ntree, importance = TRUE)
#   tail(rf_model$err.rate[, "OOB"], 1)  # Extract final OOB error rate
# })
#
# # Plot OOB error vs. Number of Trees
# tiff(filename = file.path(base.dir, "09_ROC_curve/OOB_error_Ntrees_features2.tiff"),
#      width = 3000, height = 3000, res = 300)
# plot(seq(100, 1000, by = 100), oob_errors, type = "b", col = "blue", pch = 19,
#      xlab = "Number of Trees", ylab = "OOB Error Rate",
#      main = "OOB Error vs. Number of Trees")
# dev.off()
# # End of testing the ntrees
#
# # Testing whether the feature importance is still stable even tho switching seeds
# set.seed(1)
# rf_model1 <- randomForest(asthma ~ ., data = rf_data, ntree = 500, importance = TRUE)
# set.seed(2)
# rf_model2 <- randomForest(asthma ~ ., data = rf_data, ntree = 500, importance = TRUE)
#
# top_features1 <- rownames(importance(rf_model1)[order(-importance(rf_model1)[, "MeanDecreaseGini"]), ])[1:20]
# top_features2 <- rownames(importance(rf_model2)[order(-importance(rf_model2)[, "MeanDecreaseGini"]), ])[1:20]
#
# cat("Overlap of top features:", intersect(top_features1, top_features2), "\n")
################################################################################

final_importance_file <- file.path(
  base.dir,
  "09_ROC_curve",
  "fold_level_importance_results.RDS"
)
use_random_forest <- TRUE

### Define Repeated K-Fold CV Parameters
num_folds <- 10 # Number of CV folds
repeats <- 5 # Number of repetitions for performance estimation

# Define the hyperparameter grid and train control
# These remain global since the grid and control methods are consistent across folds
rf_grid <- expand.grid(
  mtry = seq(2, 16, 2),
  splitrule = "gini",
  min.node.size = c(1, 5, 10, 20, 50, 100) # Includes smaller values for finer splits and larger values to prevent overfitting
)
# For the inner tuning on the training fold:
# 'method = "cv", number = X' here refers to the inner CV for hyperparameter tuning.
# Choose an appropriate number (e.g., 5 or 10) for inner CV folds.
# Adjust as desired.
inner_cv_control <- caret::trainControl(
  method = "cv",
  number = 10,
  allowParallel = TRUE
)

# Initialize a list to store importance results from each fold
# After running all repeats and folds, this will contain all fold-level rankings
all_fold_importance <- list()

### Global Hyperparameter Tuning and Feature Importance Extraction
if (file.exists(final_importance_file)) {
  cat("Loading previously saved fold-level importance results...\n")
  all_fold_importance_df <- readRDS(final_importance_file)
} else {
  # Example outer loop for repeats (assuming you have 'repeats' defined)
  for (r in seq_len(repeats)) {
    # Set a unique seed for each repeat
    set.seed(base_seed + r)
    cat(sprintf("Repeat %d: Seed set to %d\n", r, base_seed + r))

    # Example creation of folds for this repetition
    # 'cv_folds' should be a list of indices defining each fold’s test set
    cv_folds <- caret::createFolds(
      methylation.data.filt$asthma,
      k = num_folds,
      list = TRUE
    )

    for (fold_idx in seq_along(cv_folds)) {
      cat(sprintf(
        "Repeat %d, Fold %d: Preparing training and test data...\n",
        r,
        fold_idx
      ))

      # Define train/test indices for this fold
      test_indices <- cv_folds[[fold_idx]]
      train_indices <- setdiff(
        seq_len(nrow(methylation.data.filt)),
        test_indices
      )

      train_data <- methylation.data.filt[train_indices, ]
      test_data <- methylation.data.filt[test_indices, ]

      if (use_random_forest) {
        # Prepare the data frame for Random Forest tuning:
        rf_train_data <- train_data[, c(all_features, "asthma")]

        cat(sprintf(
          "Repeat %d, Fold %d: Running hyperparameter tuning on training data only...\n",
          r,
          fold_idx
        ))

        # Run caret::train for hyperparameter tuning using only the training data of this fold
        rf_cv_fold <- caret::train(
          asthma ~ .,
          data = rf_train_data,
          method = "ranger",
          tuneGrid = rf_grid,
          trControl = inner_cv_control,
          num.trees = 500
        )

        # Extract the best parameters for this fold
        fold_best_params <- rf_cv_fold$bestTune

        cat(sprintf(
          "Repeat %d, Fold %d: Best Params found: mtry=%d, splitrule=%s, min.node.size=%d\n",
          r,
          fold_idx,
          fold_best_params$mtry,
          fold_best_params$splitrule,
          fold_best_params$min.node.size
        ))

        # Train the final model on the entire training fold using best parameters
        final_tuneGrid <- expand.grid(
          mtry = fold_best_params$mtry,
          splitrule = fold_best_params$splitrule,
          min.node.size = fold_best_params$min.node.size
        )

        cat(sprintf(
          "Repeat %d, Fold %d: Training final model with best parameters to extract importance...\n",
          r,
          fold_idx
        ))

        final_model_fold <- caret::train(
          asthma ~ .,
          data = rf_train_data,
          method = "ranger",
          trControl = caret::trainControl(
            method = "none",
            allowParallel = TRUE
          ), # No CV here, just fit once on the full training set
          tuneGrid = final_tuneGrid,
          num.trees = 500,
          importance = "permutation"
        )

        # Extract feature importance for this fold
        fold_importance <- final_model_fold$finalModel$variable.importance

        # Sort features by importance
        fold_importance <- sort(fold_importance, decreasing = TRUE)

        # Create a rank vector: the first feature in fold_importance gets rank 1, next gets rank 2, etc.
        feature_ranks <- seq_along(fold_importance)

        # Predict on the test data
        test_predictions <- predict(final_model_fold, newdata = test_data)
        # Calculate overall accuracy on the test set
        test_accuracy <- mean(test_predictions == test_data$asthma)

        # Store fold-specific importance in a data frame with fold/repeat info
        fold_importance_df <- data.frame(
          Feature = names(fold_importance),
          Importance = as.numeric(fold_importance),
          Rank = feature_ranks,
          Test_Accuracy = test_accuracy,
          Repeat = r,
          Fold = fold_idx,
          stringsAsFactors = FALSE
        )

        # Append to the global list of all importance results
        all_fold_importance[[paste0(
          "Repeat_",
          r,
          "_Fold_",
          fold_idx
        )]] <- fold_importance_df

        cat(sprintf(
          "Repeat %d, Fold %d: Importance extracted and stored.\n",
          r,
          fold_idx
        ))
      }
    } # end of fold loop
  } # end of repeats loop
  cat("All repeats and folds completed. Aggregating importance results...\n")

  # Combine all fold-level importance data into one data frame
  all_fold_importance_df <- do.call(rbind, all_fold_importance)

  # Save the aggregated fold-level importance table if desired
  saveRDS(
    all_fold_importance_df,
    file = file.path(
      base.dir,
      "09_ROC_curve",
      "fold_level_importance_results.RDS"
    )
  )

  cat(
    "Fold-level importance extraction complete. You can now proceed with frequency-based consensus or other methods to select top stable features.\n"
  )
}

# # Define file paths for saving results
# rf_cv_file <- file.path(base.dir, "09_ROC_curve/rf_cv_results.RData")
# importance_cv_file <- file.path(base.dir, "09_ROC_curve/importance_cv_results.RData")

# # Check if files exist
# if (file.exists(rf_cv_file) && file.exists(importance_cv_file)) {
#   # Load previously saved results
#   cat("Loading precomputed Random Forest cross-validation results...\n")
#   load(rf_cv_file)
#   load(importance_cv_file)
# } else {
#   # Set up parallel backend
#   num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")) - 1
#   if (num_cores < 1) num_cores <- 1
#   cl <- parallel::makeCluster(num_cores, type = "PSOCK")
#   doParallel::registerDoParallel(cl)
#
#   # Perform computation-intensive process
#   cat("Running Random Forest cross-validation in parallel...\n")
#
#   # Cross-validating the top features
#   # K-fold Cross-validation Random Forest testing
#   control <- caret::trainControl(method = "cv", number = 10, allowParallel = TRUE)
#   rf_cv <- caret::train(asthma ~ ., data = rf_data, method = "rf", trControl = control, ntree = 500)
#
#   cat("Finished Random Forest cross-validation.\n")
#
#   # Extract feature importance from cross-validation
#   # This also takes into account the mean decrease in Accuracy (MDA)
#   #   and the mean decrease in Gini index (MDG)
#   importance_cv <- varImp(rf_cv)
#
#   # Save results for later use
#   cat("Saving Random Forest cross-validation results...\n")
#   save(rf_cv, file = rf_cv_file)
#   save(importance_cv, file = importance_cv_file)
#
#   # Clean up parallel backend
#   stopCluster(cl)
#   doParallel::stopImplicitCluster()
#
#   cat("Random Forest cross-validation and feature importance extraction complete.\n")
#   cat("Finished Random Forest\n")
#   cat("Memory usage after Random Forest: ", sum(gc()[, 2]), " MB\n")
# }

### Function to Select Top Features
get_top_features <- function(feature_table, num_features) {
  # Check if requested number of features is within range
  max_features <- length(unique(feature_table$Feature))
  if (num_features > max_features) {
    message(paste(
      "Requested number of features exceeds max (",
      max_features,
      "). Using max.",
      sep = ""
    ))
    num_features <- max_features
  } else if (num_features < 1) {
    message("Requested < 1 feature. Using 1 feature.")
    num_features <- 1
  }

  # Calculate the frequency of each feature appearing in the top 20
  feature_frequency <- aggregate(
    Rank ~ Feature,
    data = feature_table[feature_table$Rank <= 20, ],
    FUN = length
  )
  colnames(feature_frequency)[2] <- "FrequencyInTop20"

  # Calculate the median rank for each feature
  median_ranks <- aggregate(Rank ~ Feature, data = feature_table, FUN = median)
  colnames(median_ranks)[2] <- "MedianRank"

  # Merge frequency and median rank into a single data frame
  feature_stats <- merge(feature_frequency, median_ranks, by = "Feature")

  # Sort by frequency (descending) and then by median rank (ascending)
  feature_stats <- feature_stats[
    order(-feature_stats$FrequencyInTop20, feature_stats$MedianRank),
  ]

  # Select the top features based on the combined ranking
  top_features <- head(feature_stats$Feature, num_features)

  return(top_features)
}

num_features_to_use <- 20 # Specify how many features you want
top_features_cv <- get_top_features(all_fold_importance_df, num_features_to_use)

print(top_features_cv)

# Garbage collecting, due to memory intensive randomForest modelling
gc()

################################################################################
####  Visualising the top features
visualise_importance <- FALSE
if (visualise_importance == TRUE) {
  # Convert caret::varImp to a data frame for plotting
  importance_df <- as.data.frame(all_fold_importance_df)
  # importance_df$Feature <- rownames(importance_df)

  # Sort features by importance
  importance_df <- importance_df[order(-importance_df$Importance), ]

  # 1. Barplot: Visualize Top Features by Importance
  p <- ggplot(
    importance_df,
    aes(x = reorder(Feature, -Importance), y = Importance)
  ) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(
      title = "Feature Importance (Caret - Barplot)",
      x = "Features",
      y = "Normalized Importance"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(
    plot = p,
    filename = file.path(
      base.dir,
      "09_ROC_curve/barplot_important_features.png"
    ),
    dpi = 300,
    width = 4,
    height = 3
  )

  # 2. Scatterplot: Highlight Importance of Each Feature
  p <- ggplot(
    importance_df,
    aes(x = Importance, y = reorder(Feature, -Importance))
  ) +
    geom_point(color = "blue", size = 3, alpha = 0.7) +
    geom_text(aes(label = Feature), hjust = -0.2, vjust = 0, size = 3) +
    labs(
      title = "Feature Importance (Caret - Scatterplot)",
      x = "Normalized Importance",
      y = "Features"
    ) +
    theme_minimal()
  ggsave(
    plot = p,
    filename = file.path(
      base.dir,
      "09_ROC_curve/scatterplot_important_features.png"
    ),
    dpi = 300,
    width = 4,
    height = 3
  )

  # 3. Heatmap: Overview of Feature Importance
  importance_long <- melt(importance_df, id.vars = "Feature")

  p <- ggplot(importance_long, aes(x = Feature, y = variable, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue",
      high = "red",
      mid = "white",
      midpoint = mean(importance_long$value)
    ) +
    labs(
      title = "Feature Importance Heatmap (Caret)",
      x = "Features",
      y = "Metric"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  ggsave(
    plot = p,
    filename = file.path(
      base.dir,
      "09_ROC_curve/heatmap_important_features.png"
    ),
    dpi = 300,
    width = 4,
    height = 3
  )
}
################################################################################
# Initialize an empty list to store the methylation sets
methylation_sets <- list()

# Generate methylation sets
for (i in seq_along(top_features_cv)) {
  # Adding the next important feature into the set
  # methylation_sets[[paste0("set_", i)]] <- top_features_cv[1:i]

  # Important features individually as a set
  methylation_sets[[paste0("set_", i, "_single")]] <- top_features_cv[i]
}
# Example: Inspect the first few methylation sets
print(methylation_sets[2])

# To log memory usage during parallel processing
log_memory_usage <- function(step_name) {
  mem_usage <- sum(gc()[, 2]) # Sum of memory used (in MB)
  cat(sprintf("[%s] Memory usage: %.2f MB\n", step_name, mem_usage))
  return(mem_usage)
}

# Set up parallel backend
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")) - 1
if (num_cores < 1) num_cores <- 1
cl <- parallel::makeCluster(num_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

# Parameters
num_folds <- 10 # Number of folds for k-fold CV
repeats <- 5 # Number of repetitions for CV
nested_cv_results <- list() # Store results for each methylation set

results_path <- file.path(base.dir, "09_ROC_curve/confusion_matrices")
if (!dir.exists(results_path)) dir.create(results_path)

# Iterate over all methylation sets
for (set_name in names(methylation_sets)) {
  start_time <- Sys.time()
  cat("Processing methylation set:", set_name, "\n")

  # Current selected features for this set
  selected_features <- methylation_sets[[set_name]]

  # The selected features form a single "gene set" for GSVA
  methylation_set_list <- list(current_set = selected_features)

  # Initialize storage for AUCs and confusion matrices
  repeat_aucs <- numeric()
  confusion_matrices <- list()

  # Perform repeated k-fold CV
  results <- foreach(
    i = 1:repeats,
    .packages = c("caret", "GSVA", "pROC", "stats", "base")
  ) %dopar%
    {
      # Set a unique seed for each repeat derived from the base seed
      set.seed(base_seed + i)

      # Create stratified folds for this repetition (based on asthma factor)
      cv_folds <- caret::createFolds(
        methylation.data.filt$asthma,
        k = num_folds,
        list = TRUE
      )
      fold_results <- list()

      for (fold in seq_along(cv_folds)) {
        # Split data into training and test sets
        test_indices <- cv_folds[[fold]]
        train_indices <- setdiff(
          seq_len(nrow(methylation.data.filt)),
          test_indices
        )

        train_data <- methylation.data.filt[train_indices, ]
        test_data <- methylation.data.filt[test_indices, ]

        # Step 1: GSVA on training data (transform methylation features into a composite score)
        gsva_train <- gsva(
          expr = as.matrix(t(train_data[, !colnames(train_data) == "asthma"])),
          gset.idx.list = methylation_set_list,
          method = "gsva",
          min.sz = 1
        )

        # Step 2: GSVA on test data
        gsva_test <- gsva(
          expr = as.matrix(t(test_data[, !colnames(test_data) == "asthma"])),
          gset.idx.list = methylation_set_list,
          method = "gsva",
          min.sz = 1
        )

        # Prepare data for GLM
        glm_train_data <- data.frame(
          Score = gsva_train[1, ],
          Group = train_data$asthma
        )
        glm_test_data <- data.frame(
          Score = gsva_test[1, ],
          Group = test_data$asthma
        )

        # Step 3: Train a logistic GLM model
        # No hyperparameters tuned here; GLM is straightforward
        glm_model <- glm(
          Group ~ Score,
          data = glm_train_data,
          family = binomial
        )

        # Step 4: Predict probabilities on test set
        test_probs <- predict(
          glm_model,
          newdata = glm_test_data,
          type = "response"
        )

        # Step 5: Calculate ROC and determine optimal threshold
        roc_obj <- roc(glm_test_data$Group, test_probs)
        optimal_threshold_coords <- coords(
          roc_obj,
          "best",
          ret = "threshold",
          best.method = "youden"
        )
        optimal_threshold <- optimal_threshold_coords$threshold

        # Step 6: Convert probabilities to predictions using optimal threshold
        predicted_classes <- factor(
          ifelse(test_probs >= optimal_threshold, 1, 0),
          levels = c(0, 1)
        )
        actual_classes <- factor(glm_test_data$Group, levels = c(0, 1))

        # Calculate confusion matrix for this fold
        confusion_matrix <- caret::confusionMatrix(
          predicted_classes,
          actual_classes
        )

        # Store fold-level results
        fold_results[[fold]] <- list(
          Confusion_Matrix = confusion_matrix,
          AUC = auc(roc_obj),
          Fold = fold,
          Optimal_Threshold = optimal_threshold
        )
        gc()
      }
      # Return the following object from the loop and store into 'results'
      fold_results
    }

  # Consolidate results from all repeats and folds
  for (i in seq_along(results)) {
    for (fold in seq_along(results[[i]])) {
      confusion_matrices[[paste0("Repeat_", i, "_Fold_", fold)]] <- results[[
        i
      ]][[fold]]
      repeat_aucs <- c(repeat_aucs, results[[i]][[fold]]$AUC)
    }
  }

  # log_memory_usage("Before Saving Confusion Matrices")
  # Save confusion matrices for this methylation set
  # saveRDS(confusion_matrices, file = file.path(results_path, paste0("confusion_matrices_", set_name, ".RDS")))
  saveRDS(
    confusion_matrices,
    file = file.path(
      results_path,
      paste0("confusion_matrices_", set_name, "_singles.RDS")
    )
  )
  # log_memory_usage("After Saving Confusion Matrices")

  # Store average AUC for this methylation set
  nested_cv_results[[set_name]] <- list(
    Mean_AUC = mean(repeat_aucs),
    AUCs = repeat_aucs
  )

  # Record end time of the iteration
  end_time <- Sys.time()

  # Calculate elapsed time in hours for this iteration
  elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Convert to hh:mm:ss format
  hours <- floor((elapsed_seconds %% (24 * 3600)) / 3600)
  minutes <- floor((elapsed_seconds %% 3600) / 60)
  seconds <- round(elapsed_seconds %% 60)

  # Format the elapsed time as a string
  formatted_time <- sprintf("%02d:%02d:%02d", hours, minutes, seconds)

  # Print or log the formatted elapsed time
  cat(sprintf(
    "Methylation set: %s took %s (hh:mm:ss)\n",
    set_name,
    formatted_time
  ))
}

# Save the overall nested CV results
# saveRDS(nested_cv_results, file = file.path(base.dir, "09_ROC_curve/nested_cv_results.RDS"))
saveRDS(
  nested_cv_results,
  file = file.path(base.dir, "09_ROC_curve/nested_cv_results_singles.RDS")
)

# Stop parallel backend
stopCluster(cl)

cat("Cross-validation with confusion matrix calculation complete.\n")


summary_results <- data.frame(
  Methylation_Set = names(nested_cv_results),
  Mean_AUC = sapply(nested_cv_results, function(x) x$Mean_AUC),
  Fold_AUCs = sapply(
    nested_cv_results,
    function(x) paste(round(x$AUCs, 3), collapse = ", ")
  )
)

# Create a bar plot of mean AUCs
p <- ggplot(
  summary_results,
  aes(x = reorder(Methylation_Set, -Mean_AUC), y = Mean_AUC)
) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Mean AUC by Methylation Set",
    x = "Methylation Set",
    y = "Mean AUC"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# ggsave(plot = p, filename = file.path(base.dir, "09_ROC_curve/methylation_sets_performance.png"),
#        dpi = 300, width = 4, height = 3)
ggsave(
  plot = p,
  filename = file.path(
    base.dir,
    "09_ROC_curve/methylation_sets_performance_singles.png"
  ),
  dpi = 300,
  width = 4,
  height = 3
)

quit()
