# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data"
root.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan"
base.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation"

####
#### Script for Global Feature Ranking, Fixed Hyperparameter Selection, and
#### Repeated K-Fold Cross-Validation for Asthma Prediction using
#### Methylation Data and GSVA-Transformed Features
####
#### This script:
#### 1. Loads methylation and phenotype data.
#### 2. Performs a global random forest-based feature importance ranking with
####    hyperparameter tuning to identify the top features predictive of asthma.
#### 3. Selects top features and constructs methylation sets incrementally
####    (1 feature, top 2 features, top 3 features, etc.).
#### 4. Uses a fixed set of globally tuned hyperparameters from the random
####    forest step (for feature importance) and applies a simple GLM model
####    with GSVA-transformed methylation data to estimate classification
####    performance.
#### 5. Runs a repeated k-fold cross-validation to evaluate the predictive
####    performance of each methylation set and compares their AUC scores.
####
#### The final output includes confusion matrices and AUCs for each methylation set,
#### enabling the selection of the best performing methylation set.
####

### Load Required Packages
library(GSVA) # For gene set variation analysis (GSVA)
library(pROC) # For ROC curve and AUC calculations
library(randomForest) # For legacy random forest (if needed)
library(randomForestSRC) # RF-SRC (if needed)
library(ggplot2) # For data visualization
library(foreach) # For parallel loops
library(doParallel) # For parallel backend registration
library(caret) # For ML modeling and cross-validation
library(reshape2) # For data reshaping
library(ranger) # For Random Forest with tuning
library(glmnet)
library(stabs)

### Set Base Seed for Reproducibility
base_seed <- 22
set.seed(base_seed)

# methylation.data <- read.csv(file = file.path(base.dir, "006 Differential Methylated Region Analysis DMR",
#                                               "DMR_mval_AsthmaVsHealthy_allSamples.csv"),
#                              row.names = 1)
#
# pheno.data <- read.csv(file = file.path(base.dir, "100 Characteristics tables",
#                                         "ATLANTIS_linkingtable_with_group_Th.csv"),
#                        row.names = 1)
# pheno.eQTM <- read.csv(file = file.path(base.dir, "100 Characteristics tables",
#                                         "pheno_eQTM_samples.csv"),
#                        row.names = 1)
# DMR.results <- read.csv(file = file.path(base.dir, "006 Differential Methylated Region Analysis DMR",
#                                          "DMR_results_asthma_vs_control_750gap.csv"),
#                         row.names = 1)

### Assumes 'base.dir' and 'root.dir' are defined externally
### Load Input Data
# The user must ensure these files exist and are correctly formatted.
methylation.data <- read.csv(
  file = file.path(base.dir, "08_DMR/DMR_mval_AsthmaVsHealthy_allSamples.csv"),
  row.names = 1
)
pheno.data <- read.csv(
  file.path(base.dir, "06_DMA/ATLANTIS_linkingtable_with_group_Th.csv"),
  row.names = 1
)
pheno.eQTM <- read.csv(
  file.path(root.dir, "ATLANTIS_eQTM/01_data_prep/pheno_eQTM_samples.csv"),
  row.names = 1
)
DMR.results <- read.csv(
  file.path(base.dir, "08_DMR/DMR_results_asthma_vs_control_750gap.csv"),
  row.names = 1
)

### Data Preprocessing
# methylation.data <- methylationData
# rm(methylationData)
# rownames(methylation.data) <- methylation.data[, 1]
# methylation.data <- methylation.data[, -1]
colnames(methylation.data) <- gsub("X", "", colnames(methylation.data))

# Filter phenotype data to include only samples with defined asthma status
pheno.data.filt <- pheno.data[!is.na(pheno.data$ASTHEA), ]

# Filter methylation data for the subset of samples that have corresponding
# phenotype and gene expression data (as defined in pheno.eQTM)
methylation.data.filt <- methylation.data[, pheno.eQTM$meth_file_id]
# methylation.data.filt <- methylation.data[, colnames(methylation.data) %in% pheno.data.filt$meth_file_id]

# Transpose methylation data so rows = samples, columns = features
methylation.data.filt <- as.data.frame(t(methylation.data.filt))

# Characters like : and - gives problems downstream
# Make column names syntactically valid for R
colnames(methylation.data.filt) <- make.names(colnames(methylation.data.filt))

# Add asthma status (0 = Healthy, 1 = Asthma) to the data frame
methylation.data.filt$asthma <- pheno.data.filt[
  match(rownames(methylation.data.filt), pheno.eQTM$meth_file_id),
  "ASTHEA"
]
methylation.data.filt$asthma <- gsub("A", "1", methylation.data.filt$asthma)
methylation.data.filt$asthma <- gsub("H", "0", methylation.data.filt$asthma)
methylation.data.filt$asthma <- as.factor(methylation.data.filt$asthma)

### Filter DMR Results
# Include only regions meeting certain criteria (e.g. p.adjust < 0.05, n > 1)
DMR.results <- DMR.results[DMR.results$n > 1 & DMR.results$p.adjust < 0.05, ]
rownames(DMR.results) <- paste0(
  DMR.results$chr,
  ":",
  DMR.results$start,
  "-",
  DMR.results$end
)
rownames(DMR.results) <- make.names(rownames(DMR.results))
all_features <- rownames(DMR.results)

head(all_features)
cat(
  "Memory usage before Random Forest model fitting: ",
  sum(gc()[, 2]),
  " MB\n"
)

### Prepare Data for Random Forest-Based Global Feature Ranking
# Here we use the globally optimal hyperparameters from RF to rank features
rf_data <- methylation.data.filt[, c(all_features, "asthma")]

### Define File Paths for Saving/Loading RF Results
# rf_cv_file <- file.path(base.dir, "09_ROC_curve/rf_cv_results_DMR.RData")
# importance_cv_file <- file.path(base.dir, "09_ROC_curve/importance_cv_results_DMR.RData")
# best_params_file <- file.path(base.dir, "09_ROC_curve/best_feature_params_DMR.RData")
final_importance_file <- file.path(
  base.dir,
  "09_ROC_curve",
  "fold_level_importance_results_DMR.RDS"
)
use_random_forest <- TRUE

### Define Repeated K-Fold CV Parameters
num_folds <- 10 # Number of CV folds
repeats <- 5 # Number of repetitions for performance estimation
# Define the hyperparameter grid and train control
# These remain global since the grid and control methods are consistent across folds
rf_grid <- expand.grid(
  mtry = seq(2, 34, 2),
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
      } else {
        #### Might want to change this into Elastic Net regression
        ####   Which combines the Lasso L1 and Ridge L2 methods
        ####   Shrinking variables for variable selection
        ####   But also accounting for possible multicollinearity (multiple correlations between features)
        ####    Which is highly possible with methylation sites
        # Ridge implementation
        x_train <- as.matrix(train_data[, all_features])
        y_train <- as.numeric(as.character(train_data$asthma))

        cat(sprintf("Repeat %d, Fold %d: Running Ridge...\n", r, fold_idx))
        ridge_model <- cv.glmnet(
          x_train,
          y_train,
          family = "binomial", # For binary outcomes
          alpha = 0, # Alpha = 1 for Lasso, Alpha = 0 for Ridge
          nfolds = 10 # Inner cross-validation folds for lambda tuning
        )

        # Extract the best lambda
        best_lambda <- ridge_model$lambda.min

        # print(best_lambda)
        # Refit Lasso model on full training data using best lambda
        final_ridge_model <- glmnet(
          x_train,
          y_train,
          family = "binomial",
          alpha = 0,
          lambda = best_lambda
        )

        # Extract non-zero coefficients and their ranks
        ridge_coefs <- as.numeric(as.matrix(coef(final_ridge_model)))
        # print(lasso_coefs)
        names(ridge_coefs) <- rownames(coef(final_ridge_model))

        # Removing the Intercept coefficient
        ridge_coefs <- ridge_coefs[
          names(ridge_coefs) != "(Intercept)",
          drop = FALSE
        ]
        ranked_features <- names(ridge_coefs)[order(
          ridge_coefs,
          decreasing = TRUE
        )]

        # Predicting and assessing performance on test_data
        predictions <- predict(
          final_lasso_model,
          newx = as.matrix(test_data[, -which(names(test_data) == "asthma")])
        )
        actual_values <- as.numeric(test_data[["asthma"]])

        ss_total <- sum((actual_values - mean(actual_values))^2) # Total sum of squares
        ss_residual <- sum((actual_values - predictions)^2) # Residual sum of squares
        r_squared <- 1 - (ss_residual / ss_total) # R-squared
        rmse <- sqrt(mean((actual_values - predictions)^2)) # RMSE

        # Save all interesting parameters in 1 table
        fold_lasso_ranks <- data.frame(
          Feature = ranked_features,
          Coefficient = ridge_coefs,
          Rank = rank(-ridge_coefs),
          Repeat = r,
          Fold = fold_idx,
          R_squared = r_squared,
          RMSE = rmse,
          stringsAsFactors = FALSE
        )

        all_fold_importance[[paste0(
          "Repeat_",
          r,
          "_Fold_",
          fold_idx
        )]] <- fold_lasso_ranks
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
      "fold_level_importance_results_DMR.RDS"
    )
  )

  cat(
    "Fold-level importance extraction complete. You can now proceed with frequency-based consensus or other methods to select top stable features.\n"
  )
}

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

### Define the Number of Top Features to Investigate (Clinical Constraint)
num_features_to_use <- 20
top_features_cv <- get_top_features(all_fold_importance_df, num_features_to_use)
print(top_features_cv)

# Garbage Collection after heavy computations
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
  # print(p)
  ggsave(
    plot = p,
    filename = file.path(
      base.dir,
      "09_ROC_curve/barplot_important_features_DMR.png"
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
    # geom_text(aes(label = Feature), hjust = -0.2, vjust = 0, size = 3) +
    labs(
      title = "Feature Importance (Caret - Scatterplot)",
      x = "Normalized Importance",
      y = "Features"
    ) +
    theme_minimal()
  # print(p)
  ggsave(
    plot = p,
    filename = file.path(
      base.dir,
      "09_ROC_curve/scatterplot_important_features_DMR.png"
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
  # print(p)
  ggsave(
    plot = p,
    filename = file.path(
      base.dir,
      "09_ROC_curve/heatmap_important_features_DMR.png"
    ),
    dpi = 300,
    width = 4,
    height = 3
  )
}
################################################################################
### Create Methylation Sets from Top Features
# set_1 contains the top 1 feature, set_2 the top 2 features, etc.
methylation_sets <- list()
for (i in seq_along(top_features_cv)) {
  # Appending the next important feature into the set
  # methylation_sets[[paste0("set_", i)]] <- top_features_cv[1:i]

  # Important features individually as a set
  methylation_sets[[paste0("set_", i, "_single")]] <- top_features_cv[i]
}
print(methylation_sets[2]) # Example inspection

### Set up Parallel Backend for CV
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")) - 1
if (num_cores < 1) num_cores <- 1
cl <- parallel::makeCluster(num_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

### Define Repeated K-Fold CV Parameters
num_folds <- 10 # Number of CV folds
repeats <- 5 # Number of repetitions for performance estimation

# We will store performance results for each methylation set here
cv_results <- list()

# Directory for confusion matrices and results
results_path <- file.path(base.dir, "09_ROC_curve/confusion_matrices")
if (!dir.exists(results_path)) dir.create(results_path)

### Evaluate Each Methylation Set Using Repeated K-Fold CV
# Note: At this point, hyperparameters and feature ranking are fixed globally.
#       Each CV run is just to estimate how well these features perform in a GLM model.

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

  # Save confusion matrices for this methylation set
  # saveRDS(confusion_matrices, file = file.path(results_path, paste0("confusion_matrices_DMR", set_name, ".RDS")))
  saveRDS(
    confusion_matrices,
    file = file.path(
      results_path,
      paste0("confusion_matrices_DMR", set_name, "_singles.RDS")
    )
  )

  # Store mean AUC across all folds/repeats
  cv_results[[set_name]] <- list(
    Mean_AUC = mean(repeat_aucs),
    AUCs = repeat_aucs
  )

  # Timing information
  end_time <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(end_time, start_time, units = "secs"))
  hours <- floor((elapsed_seconds %% (24 * 3600)) / 3600)
  minutes <- floor((elapsed_seconds %% 3600) / 60)
  seconds <- round(elapsed_seconds %% 60)
  formatted_time <- sprintf("%02d:%02d:%02d", hours, minutes, seconds)

  cat(sprintf(
    "Methylation set: %s took %s (hh:mm:ss)\n",
    set_name,
    formatted_time
  ))
}

# Save the overall nested CV results
# Save the summarized results for all sets
# saveRDS(cv_results, file = file.path(base.dir, "09_ROC_curve/nested_cv_results_DMR.RDS"))
saveRDS(
  cv_results,
  file = file.path(base.dir, "09_ROC_curve/nested_cv_results_singles_DMR.RDS")
)

# Stop parallel backend
stopCluster(cl)

cat("Repeated k-fold CV evaluation complete.\n")


### Summarize Results
summary_results <- data.frame(
  Methylation_Set = names(cv_results),
  Mean_AUC = sapply(cv_results, function(x) x$Mean_AUC),
  Fold_AUCs = sapply(
    cv_results,
    function(x) paste(round(x$AUCs, 3), collapse = ", ")
  )
)

### Visualize Methylation Set Performance
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
# ggsave(plot = p, filename = file.path(base.dir, "09_ROC_curve/methylation_sets_performance_DMR.png"),
#        dpi = 300, width = 4, height = 3)
ggsave(
  plot = p,
  filename = file.path(
    base.dir,
    "09_ROC_curve/methylation_sets_performance_singles_DMR.png"
  ),
  dpi = 300,
  width = 4,
  height = 3
)

cat("Analysis complete. Results and figures have been saved.\n")

quit()

# Run locally to inspect/examine the results
all_fold_importance_df <- readRDS(
  "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/009 ROC curves/fold_level_importance_results_DMR.RDS"
)

test <- readRDS(
  "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/009 ROC curves/confusion_matrices/confusion_matrices_DMRset_1.RDS"
)
nested_cv_results <- readRDS(
  "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/009 ROC curves/nested_cv_results.RDS"
)
nested_cv_results <- readRDS(
  "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/009 ROC curves/nested_cv_results_singles.RDS"
)
nested_cv_results <- readRDS(
  "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/009 ROC curves/nested_cv_results_DMR.RDS"
)

# Convert nested_cv_results to a data frame
nested_results_df <- do.call(
  rbind,
  lapply(names(nested_cv_results), function(set_name) {
    # Extract information for this methylation set
    set_results <- nested_cv_results[[set_name]]

    # Flatten results into a data frame
    data.frame(
      Methylation_Set = set_name,
      Mean_AUC = set_results$Mean_AUC,
      Fold_AUCs = paste(round(set_results$AUCs, 3), collapse = ", ") # Combine fold AUCs into a string
    )
  })
)

nested_results_long <- do.call(
  rbind,
  lapply(1:nrow(nested_results_df), function(i) {
    aucs <- as.numeric(strsplit(nested_results_df$Fold_AUCs[i], ",\\s*")[[1]])
    data.frame(
      Methylation_Set = nested_results_df$Methylation_Set[i],
      Fold_AUC = aucs
    )
  })
)

nested_results_long$Methylation_Set <- factor(
  nested_results_long$Methylation_Set,
  levels = nested_results_df[
    order(-nested_results_df$Mean_AUC),
    "Methylation_Set"
  ]
)


ggplot(nested_results_long, aes(x = Methylation_Set, y = Fold_AUC)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 19,
    size = 2,
    color = "red"
  ) +
  labs(
    title = "Distribution of Fold AUCs Across Methylation Sets",
    x = "Methylation Set",
    y = "AUC"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


auc_summary <- data.frame(
  Set = names(nested_cv_results),
  Mean_AUC = sapply(nested_cv_results, function(x) x$Mean_AUC)
)

ggplot(auc_summary, aes(x = reorder(Set, -Mean_AUC), y = Mean_AUC)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(
    title = "Mean AUC Across Methylation Sets",
    x = "Methylation Set",
    y = "Mean AUC"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


load(
  "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/009 ROC curves/1. Methylation/rf_cv_results.RData"
)
load(
  "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/009 ROC curves/1. Methylation/importance_cv_results.RData"
)
