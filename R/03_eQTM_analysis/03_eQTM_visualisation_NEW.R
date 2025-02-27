# -----------------------------------------------------------------------------
# File: 03_eQTM_visualisation.R
# Section 1: Library Loading, Data Loading, and Preprocessing
# -----------------------------------------------------------------------------

# Define paths to libraries and directories dynamically
set_paths <- function(is_hpc = FALSE) {
  if (!is_hpc) {
    library_path <- "C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages"
    .libPaths(library_path)
  }
  
  list(
    root_dir = if (is_hpc) "/path/to/HPC/root" else "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS",
    base_dir = if (is_hpc) "/path/to/HPC/base" else "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/005 eQTM analysis",
    input_dir = if (is_hpc) "/path/to/HPC/input" else "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/002 Differential Methylation Analysis DMA/",
    gene_location_file = if (is_hpc) "/path/to/HPC/geneLocation_63k.csv" else "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/005 eQTM analysis/geneLocation_63k.csv"
  )
}

# Load required libraries
load_required_libraries <- function() {
  library("ggplot2")
  library("dplyr")
  library("limma")
  library("ggpubr")
  library("openintro")
  library("patchwork")
  library("scales")
  library("PupillometryR")
  library("lubridate")
  # Source external functions
  source("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Documentation/Github/004 Bioinformatics/Scripts/M2beta.R")
  
}

# Load RData file and return objects
load_rdata <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found:", file_path)
  }
  
  # Create a custom environment to isolate loaded objects
  custom_env <- new.env()
  
  # Load the .rdata file into the custom environment
  load(file_path, envir = custom_env)
  
  # Retrieve the names and values of loaded objects
  loaded_objects <- as.list(custom_env)
  
  # Print loaded object names for debugging
  cat("Successfully loaded:", file_path, "\n")
  cat("Loaded objects:", names(loaded_objects), "\n")
  
  return(loaded_objects)
}

# Load all necessary data files
load_data <- function(paths) {
  # Load RData files using `load_rdata`
  eqtm_data <- load_rdata(file.path(paths$base_dir, "toEQTM_asthma_vs_control.rdata"))
  
  # Read CSV files
  meth_data <- read.csv(file = file.path(paths$base_dir, "methylationData_eQTM_samples_asthma_vs_control.csv"),
                        header = TRUE, row.names = 1)
  # Transform methylation data from M-values -> beta values
  # For easier interpretation of figures)
  meth_data <- M2beta(meth_data)
  
  expr_data <- read.csv(file = file.path(paths$base_dir, "normalized_geneExpressionData_eQTM_asthma_vs_control_samples.csv"),
                        header = TRUE, row.names = 1)
  DMA_result <- read.csv(file = file.path(paths$input_dir, "Correct for PC1-3/Tt2_significant_CpG_asthma_vs_control_BH.csv"),
                         header = TRUE, row.names = 1)
  # DMA_result <- read.csv(file = file.path(paths$input_dir, "Correct for PC1-3 - Th2 groups/DMA_Tt2_significant_CpG_Th2HighVsHealthy_BH.csv"),
  #                        header = TRUE, row.names = 1)
  Th2_groups <- read.csv(file = file.path(paths$base_dir, "covariates_eQTM_samples_IDs.csv"),
                         header = TRUE, row.names = 1)
  
  list(
    eqtm_data = eqtm_data,  # Contains objects loaded from the RData file
    meth_data = meth_data,
    expr_data = expr_data,
    DMA_result = DMA_result,
    covariates = Th2_groups,
    Th2_groups = Th2_groups
  )
}

# Preprocess data
preprocess_data <- function(raw_data, base_dir) {
  # Clean column names for methylation and expression data
  colnames(raw_data$meth_data) <- gsub(pattern = "X", replacement = "", x = colnames(raw_data$meth_data))
  colnames(raw_data$expr_data) <- gsub(pattern = "X", replacement = "", x = colnames(raw_data$expr_data))
  
  # Filter covariates for matching patients
  raw_data$covariates <- raw_data$covariates[order(raw_data$covariates$PT), ]
  raw_data$covariates <- raw_data$covariates[raw_data$covariates$meth_file_id %in% colnames(raw_data$meth_data), ]
  
  # Process covariates
  covariates <- raw_data$covariates
  
  # Modify the ASTHEA column
  if ("ASTHEA" %in% colnames(covariates)) {
    covariates$ASTHEA <- gsub("A", "asthma", covariates$ASTHEA)
    covariates$ASTHEA <- gsub("H", "healthy", covariates$ASTHEA)
  }
  
  # Modify the group_th column
  if ("group_th" %in% colnames(covariates)) {
    covariates$group_th <- gsub("high", "Th2-high", covariates$group_th)
    covariates$group_th <- gsub("low", "Th2-low", covariates$group_th)
    # Leave "healthy" and "undeterm" unchanged
  }
  
  # Filter methylation and expression data based on covariates
  meth_data_filtered <- raw_data$meth_data[, raw_data$covariates$meth_file_id]
  expr_data_filtered <- raw_data$expr_data[, raw_data$covariates$GenomeScan_ID]
  
  # Load the gene conversion table
  gene_conversion_table <- read.csv(file = paths$gene_location_file,
                                    header = TRUE, row.names = 1, na.strings = "")
  gene_conversion_table <- gene_conversion_table[complete.cases(gene_conversion_table), ]
  
  # Create a conversion vector from ENSEMBL IDs to gene symbols
  conversion_vector <- setNames(gene_conversion_table$hgnc_symbol, gene_conversion_table$ensembl_gene_id)
  
  # Convert rownames of expression data
  new_rownames <- ifelse(
    rownames(expr_data_filtered) %in% names(conversion_vector),
    conversion_vector[rownames(expr_data_filtered)],
    rownames(expr_data_filtered)
  )
  new_rownames <- make.unique(new_rownames, sep = ".")
  
  rownames(expr_data_filtered) <- new_rownames
  
  
  list(
    meth_data = meth_data_filtered,
    expr_data = expr_data_filtered,
    covariates = covariates,
    conversion_vector = conversion_vector,
    gene_conversion_table = gene_conversion_table
  )
}

# -----------------------------------------------------------------------------
# Section 2: CpG Selection and Processing
# Purpose: Select and process specific CpG sites for further analysis
# -----------------------------------------------------------------------------

# Function to determine the number of CpGs to process
select_cpg_count <- function(DMA_result, top_n = NULL) {
  if (is.null(top_n)) {
    # Use all CpGs if no specific number is provided
    nrow(DMA_result)
  } else {
    # Use the top N CpGs, ensuring it doesn't exceed the total available
    min(top_n, nrow(DMA_result))
  }
}

# Filter and process CpG data
process_cpg_data <- function(DMA_result, finaltable, expr_data, conversion_vector, use_all_genes = FALSE, top_n = NULL) {
  # Determine the number of CpGs to process
  no.CpGs <- select_cpg_count(DMA_result, top_n)
  
  # Select CpG sites from DMA result
  filtered_cpg_table <- finaltable[finaltable$snps %in% rownames(DMA_result)[1:no.CpGs], ]
  
  # Order the filtered eQTM list with the order of significance from DMA result
  filtered_cpg_table <- filtered_cpg_table[match(rownames(DMA_result)[1:no.CpGs], filtered_cpg_table$snps), ]
  
  # Remove entries with NA (CpG sites not nearby a gene in the DNA)
  filtered_cpg_table <- na.omit(filtered_cpg_table)
  
  # Filter the gene expression data for the remaining genes
  eQTM_GE_data <- expr_data[rownames(expr_data) %in% filtered_cpg_table$gene, ]
  
  # Optionally save all expression data to eQTM_GE_data for future use
  if (use_all_genes) {
    eQTM_GE_data <- expr_data
  }
  
  # Convert ENSEMBL IDs in filtered CpG table to gene symbols
  symbol_vector <- character(length(filtered_cpg_table$gene))
  multiple_matches <- list()
  
  for (i in seq_along(filtered_cpg_table$gene)) {
    gene <- filtered_cpg_table$gene[i]
    found_index <- grep(pattern = gene, x = names(conversion_vector))
    
    if (length(found_index) == 1) {
      # Use the gene symbol if exactly one match is found
      symbol <- conversion_vector[gene]
      if (is.na(symbol)) {
        symbol <- gene
        cat("This ENSEMBL ID has no gene symbol\n")
      }
    } else if (length(found_index) > 1) {
      # Save multiple matches for review
      multiple_matches[[gene]] <- conversion_vector[found_index]
      symbol <- gene
    } else {
      # Retain original gene ID if no match is found
      symbol <- gene
    }
    
    symbol_vector[i] <- symbol
  }
  
  # Add the symbols to the filtered CpG table
  filtered_cpg_table$symbol <- symbol_vector
  
  valid_pairs <- data.frame(
    CpG = filtered_cpg_table$snps,
    Gene = filtered_cpg_table$symbol
  )
  
  list(
    filtered_cpg_table = filtered_cpg_table,     # Table with filtered CpGs and their gene associations
    gene_expression_data = eQTM_GE_data,         # Filtered expression data
    multiple_matches = multiple_matches,         # Debugging info for ENSEMBL-to-symbol mapping
    valid_pairs = valid_pairs                    # Table of CpG-Gene pairs for visualization
  )
}

# Analyze CpG-Gene comparisons dynamically
analyze_comparison <- function(
    filtered_cpg_table, meth_data, covariates, DMA_result, 
    group_of_interest, comparison_type = "asthma_vs_healthy", 
    threshold = NULL, root_dir = ""
) {
  # Validate group_of_interest
  valid_groups <- c("asthma", "Th2-high", "Th2-low")
  if (!group_of_interest %in% valid_groups) {
    stop("Invalid group_of_interest. Must be one of: 'asthma', 'Th2-high', 'Th2-low'.")
  }
  
  # Initialize results data frame
  results <- data.frame(
    CpG = character(),
    Mean_Diff_1 = numeric(),
    Mean_Diff_2 = numeric(),
    Combined_Diff = numeric(),
    adj.P.Val = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each CpG
  for (i in seq_len(nrow(filtered_cpg_table))) {
    cpg_name <- filtered_cpg_table[i, "snps"]
    
    # Extract methylation data
    df <- t(meth_data[cpg_name, , drop = FALSE])
    df <- as.data.frame(df)
    
    # Dynamically select the appropriate column for group assignment
    if (group_of_interest == "asthma") {
      # Use the "ASTHEA" column for asthma comparisons
      df$group <- covariates[covariates$meth_file_id %in% rownames(df), "ASTHEA"]
    } else if (group_of_interest %in% c("Th2-high", "Th2-low")) {
      # Use the "group_th" column for endotype comparisons
      df$group <- covariates[covariates$meth_file_id %in% rownames(df), "group_th"]
    } else {
      stop("Invalid group_of_interest. Must be 'asthma', 'Th2-high', or 'Th2-low'.")
    }
    
    # Validate presence of groups
    if (!all(c("healthy", group_of_interest) %in% df$group)) {
      next  # Skip if required groups are missing
    }
    
    # Calculate group means
    group_means <- tapply(df[, 1], df$group, mean, na.rm = TRUE)
    
    # Comparison logic
    if (comparison_type == "asthma_vs_healthy" && group_of_interest == "asthma") {
      # Asthma vs Healthy: Single comparison
      mean_diff <- group_means["asthma"] - group_means["healthy"]
      results <- rbind(results, data.frame(
        CpG = cpg_name,
        Mean_Diff_1 = mean_diff,
        Mean_Diff_2 = NA,
        Combined_Diff = mean_diff,  # Only one difference in this case
        adj.P.Val = DMA_result[cpg_name, "adj.P.Val"],
        stringsAsFactors = FALSE
      ))
    } else if (comparison_type == "endotype") {
      # Endotype comparisons
      if (!all(c("Th2-low", "Th2-high") %in% df$group)) {
        next  # Skip if required endotype groups are missing
      }
      mean_diff_1 <- group_means[group_of_interest] - group_means["healthy"]
      mean_diff_2 <- group_means[group_of_interest] - group_means[setdiff(c("Th2-high", "Th2-low"), group_of_interest)]
      combined_diff <- abs(mean_diff_1) - abs(mean_diff_2)
      results <- rbind(results, data.frame(
        CpG = cpg_name,
        Mean_Diff_1 = mean_diff_1,
        Mean_Diff_2 = mean_diff_2,
        Combined_Diff = combined_diff,
        adj.P.Val = DMA_result[cpg_name, "adj.P.Val"],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Apply threshold filtering if specified
  if (!is.null(threshold)) {
    results <- results[abs(results$Combined_Diff) >= threshold, ]
  }
  
  # Sort results
  if (comparison_type == "asthma_vs_healthy") {
    results <- results[order(-abs(results$Mean_Diff_1)), ]  # Largest difference first
  } else if (comparison_type == "endotype") {
    results <- results[order(abs(results$Combined_Diff)), ]  # Closest to zero first
  }
  
  # Save results
  # write.table(
  #   results, 
  #   file = file.path(root_dir, paste0("/Results/", group_of_interest, "_Comparison_Results.txt")),
  #   quote = FALSE, row.names = FALSE, sep = "\t"
  # )
  
  return(results)
  list(
    results = results,
    valid_pairs = eQTM_GE_data,
    multiple_matches = multiple_matches  # Useful for debugging or review
  )
}

# -----------------------------------------------------------------------------
# Function: visualize_cpg_gene_pairs
# Purpose: Generate visualizations for CpG-Gene pairs with dynamic grouping
# -----------------------------------------------------------------------------
# Functions adapted from https://stackoverflow.com/a/39611375
#   used for the nearest decimal to edit the axes limits
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

# Visualize CpG-Gene Pairs
visualize_cpg_gene_pairs <- function(valid_pairs, meth_data, gene_expression_data, covariates, 
                                     plot_type = "patchwork", group_of_interest = "asthma", 
                                     add_stats = TRUE, plots_per_figure = 9) {
  if (!plot_type %in% c("patchwork", "individual")) {
    stop("Invalid plot_type. Choose either 'patchwork' or 'individual'.")
  }
  
  # Validate group_of_interest
  valid_groups <- c("asthma", "Th2-high", "Th2-low")
  if (!group_of_interest %in% valid_groups) {
    stop("Invalid group_of_interest. Must be one of: 'asthma', 'Th2-high', or 'Th2-low'.")
  }
  
  # Filter the valid_pairs table based on the CpG column in valid_pairs and asthma_results
  # valid_pairs <- valid_pairs[valid_pairs$CpG %in% asthma_results$CpG, ]
  
  
  # Initialize the plot list
  plot_list <- list()
  
  for (i in seq_len(nrow(valid_pairs))) {
    cpg_name <- valid_pairs$CpG[i]
    gene_name <- valid_pairs$Gene[i]
    
    # Extract methylation and gene expression data
    df <- t(meth_data[cpg_name, , drop = FALSE])
    df <- cbind(df, t(gene_expression_data[gene_name, , drop = FALSE]))
    df <- as.data.frame(df)
    colnames(df) <- c("cpg", "gene")
    
    # Dynamically map the group column
    if (group_of_interest == "asthma") {
      df$group <- covariates[covariates$meth_file_id %in% rownames(df), "ASTHEA"]
      df$group <- gsub("healthy", "Healthy", df$group)
      df$group <- gsub("asthma", "Asthma", df$group)
      color_palette <- c("Healthy" = "#1f77b4", "Asthma" = "#ff7f0e")
      df$group <- factor(df$group, levels = c("Healthy", "Asthma"))
    } else {
      df$group <- covariates[covariates$meth_file_id %in% rownames(df), "group_th"]
      df$group <- gsub("healthy", "Healthy", df$group)
      df$group <- gsub("Th2-high", "Th2-high", df$group)
      df$group <- gsub("Th2-low", "Th2-low", df$group)
      df$group <- gsub("undeterm", "Undetermined subtype", df$group)
      color_palette <- c("Healthy" = "#1f77b4", "Th2-low" = "#ffb74d", "Th2-high" = "#e6550d", "Undetermined subtype" = "black")
      df$group <- factor(df$group, levels = c("Healthy", "Th2-low", "Undetermined subtype", "Th2-high"))
    }
    
    # Set axis limits
    max_cpg <- ceiling_dec(max(df$cpg), 1)
    min_cpg <- floor_dec(min(df$cpg), 1)
    max_gene <- ceiling_dec(max(df$gene), 1)
    min_gene <- floor_dec(min(df$gene), 1)
    
    
    # Create eQTM plot
    p_eQTM <- ggplot(df, aes(x = gene, y = cpg)) +
      geom_point(aes(colour = group), size = 1.3, alpha = 0.9, stroke = 0.7) +
      geom_smooth(linetype = "dashed", se = FALSE, color = 'black', formula = 'y ~ x', method = "lm") +
      # stat_cor(r.accuracy = 0.01, p.accuracy = 0.001, label.x.npc = 0.70, label.y.npc = 0.88, size = 5) +
      # stat_regline_equation(label.x.npc = 0.70, label.y.npc = 0.95, size = 5) +
      scale_color_manual(values = color_palette) +
      scale_x_continuous(limits = c(min_gene, max_gene), breaks = pretty_breaks(5)) +
      scale_y_continuous(limits = c(min_cpg, max_cpg), breaks = pretty_breaks(5)) +
      guides(colour = guide_legend(title="Status")) +
      ## Commented these out, nicer looking for final figure
      # labs(x = paste0(gene_name, " Expression (normalized logCPM)"),
      #      y = paste0(cpg_name, " Methylation (Beta)"),
      #      color = "Group") +
      labs(x = paste0(gene_name),
           y = paste0(cpg_name),
           color = "Group") +
      theme_bw() + 
      theme(legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text.x = element_text(size = 14),
            axis.title = element_text(size = 14), ## Changed the size, nicer looking for final figure
            axis.text.y = element_text(size = 14))
    
    # Conditionally add stat_cor and stat_regline_equation
    if (add_stats) {
      p_eQTM <- p_eQTM +
        stat_cor(r.accuracy = 0.01, p.accuracy = 0.001, label.x.npc = 0.70, label.y.npc = 0.88, size = 5) +
        stat_regline_equation(label.x.npc = 0.70, label.y.npc = 0.95, size = 5)
    }
    
    
    # Create CpG plot
    p_CpG <- ggplot(df, aes(x = group, y = cpg, fill = group)) +
      geom_flat_violin(position = position_nudge(x = 0.20), alpha = 0.8) +
      geom_point(aes(colour = group), position = position_jitter(width = 0.15), alpha = 0.5, size = 1.3, show.legend = F) +
      geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
      scale_color_manual(values = color_palette) +
      scale_fill_manual(values = color_palette) +
      scale_y_continuous(limits = c(min_cpg, max_cpg), breaks = pretty_breaks(5)) +
      labs(x = "Status") +
      theme_bw() + 
      theme(legend.position = "none",
            axis.ticks.x=element_blank(), 
            axis.title.x = element_blank(),
            axis.text.x = element_blank(), #
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    
    
    # Create gene expression plot
    p_gene <- ggplot(df, aes(x = group, y = gene, fill = group)) +
      geom_flat_violin(position = position_nudge(x = 0.20), alpha = 0.8) +
      geom_point(aes(colour = group), position = position_jitter(width = 0.15), alpha = 0.5, size = 1.3, show.legend = F) +
      geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
      scale_fill_manual(values = color_palette, guide = NULL) +
      scale_color_manual(values = color_palette) +
      scale_y_continuous(limits = c(min_gene, max_gene), breaks = pretty_breaks(5)) +
      labs(x = "Status") +
      coord_flip() + 
      theme_bw() + 
      theme(legend.position = "none", 
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(), 
            axis.title.y = element_blank(),
            axis.text.y = element_blank(), #
            axis.ticks.y=element_blank(),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
      )
      
    
    # Combine plots for individual option
    if (plot_type == "individual") {
      layout <- 
      "
      222222222##
      222222222##
      33333333311
      33333333311
      33333333311
      33333333311
      33333333311
      33333333311
      33333333311
      "
      p_all <- p_CpG + p_gene + p_eQTM + plot_layout(design = layout)
      plot_list[[i]] <- p_all
      
    } else if (plot_type == "patchwork") {
      plot_list[[i]] <- p_eQTM
    }
  }
  
  if (plot_type == "patchwork") {
    # Split plots into chunks of `plots_per_figure`
    plot_chunks <- split(plot_list, ceiling(seq_along(plot_list) / plots_per_figure))
    
    # Combine each chunk into a separate patchwork figure
    combined_figures <- lapply(plot_chunks, function(chunk) {
      wrap_plots(chunk, ncol = 2, nrow = 5)
    })
    
    return(combined_figures) 
  } else if (plot_type == "individual") {
    
    return(plot_list)
  }
}

# -----------------------------------------------------------------------------
# Main Execution Workflow
# -----------------------------------------------------------------------------

# Step 1: Set Paths
# Dynamically set paths based on whether running on HPC or local
paths <- set_paths(is_hpc = FALSE)

# Step 2: Load Required Libraries
# Load all necessary R libraries for the script
load_required_libraries()

# Step 3: Load Data
# Load the raw data from specified paths
raw_data <- load_data(paths)

# Step 4: Preprocess Data
# Perform preprocessing tasks such as ensuring unique row names
processed_data <- preprocess_data(raw_data, paths)

# Step 5: CpG Selection and Processing
# Use the DMA results and final table to filter and select CpG sites
top_n <- NULL  # Specify the number of CpGs to select
cpg_data <- process_cpg_data(
  DMA_result = raw_data$DMA_result,
  finaltable = raw_data$eqtm_data$finaltable,
  expr_data = processed_data$expr_data,
  conversion_vector = processed_data$conversion_vector,
  use_all_genes = FALSE,  # Use only filtered genes
  top_n = top_n  # Specify the number of CpGs to select
)


# Step 6: Analyze CpG-Gene Comparisons
# Perform analysis based on the group of interest (e.g., "asthma" or "Th2-high")
asthma_results <- analyze_comparison(
  filtered_cpg_table = cpg_data$filtered_cpg_table,
  meth_data = processed_data$meth_data,
  covariates = processed_data$covariates,
  DMA_result = raw_data$DMA_result,
  group_of_interest = "asthma",  # Comparison group
  comparison_type = "asthma_vs_healthy",  # "asthma_vs_healthy" or "endotype"
  threshold = 0.05,  # Filter for mean differences >= 0.01
  root_dir = paths$base_dir
)

## Change the pairs here of own interest
# cpg_data$valid_pairs # Table with "CpG" and "Gene" as colnames
cpg_data$valid_pairs <- cpg_data$valid_pairs[match(c("cg01482588", "cg21501207", 
                                                     "cg01482377", "cg02333649", 
                                                     "cg26554722", "cg24224501", 
                                                     "cg07010930", "cg19398575", 
                                                     "cg22689016", "cg25710507"), 
                                                   cpg_data$valid_pairs$CpG), ]

cpg_data$valid_pairs <- cpg_data$valid_pairs[match(c("cg01482588", "cg21501207"), 
                                                   cpg_data$valid_pairs$CpG), ]

cpg_data$valid_pairs <- data.frame(CpG = "cg22719314", Gene = "FLG")

# When combining the correlation plots together
# Generate and save the visualization
combined_plot <- visualize_cpg_gene_pairs(
  valid_pairs = cpg_data$valid_pairs,
  meth_data = processed_data$meth_data,
  gene_expression_data = processed_data$expr_data,
  covariates = processed_data$covariates,
  plot_type = "patchwork",
  group_of_interest = "asthma",
  add_stats = TRUE,
  plots_per_figure = 10
)
print(combined_plot[[1]])

# Display or save each figure
for (i in seq_along(combined_plot)) {
  # print(combined_plot[[i]])  # Display in R
  ggsave(plot = combined_plot[[i]],
         filename = file.path(paths$root_dir, "Results/005 eQTM analysis",
                              paste0("patchwork_figure_", i, ".tiff"))
  )
}

# Save plot to a file
tiff(file.path(paths$root_dir, "Results/005 eQTM analysis/CpG_Gene_Plots2.tiff"), 
     width = 8.5, height = 17, units = "in", res = 300)
print(combined_plot[[1]])
dev.off()



# When to extract the individual eQTM figures (correlation + marginal)
individual_plots <- visualize_cpg_gene_pairs(
  valid_pairs = cpg_data$valid_pairs,
  meth_data = processed_data$meth_data,
  gene_expression_data = processed_data$expr_data,
  covariates = processed_data$covariates,
  plot_type = "individual",
  group_of_interest = "asthma",
  add_stats = TRUE,
  plots_per_figure = 10
)
print(individual_plots[[6]])

# Display or save each figure
for (i in seq_along(individual_plots)) {
  # print(individual_plots[[i]])  # Display in R
  ggsave(plot = individual_plots[[i]],
         filename = file.path(paths$root_dir, "Results/005 eQTM analysis",
                              paste0("individual_figure_", i, ".tiff"))
  )
  dev.off()
}
# Save plot to a file
tiff(file.path(paths$root_dir, "Results/005 eQTM analysis/Gene_plot2.tiff"), 
     width = 1250, height = 1250, res = 300)
print(individual_plots[[2]])
dev.off()

####



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
# output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/"

##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("ggplot2")
library("dplyr")
library("limma")

library("ggpubr")

library("openintro")
library("patchwork")
library("scales")
library("PupillometryR")

library("lubridate")
##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
# Set input folder
root.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/005 eQTM analysis"

# Load in eQTM RData
load(file = file.path(base.dir, "toEQTM_Th2HighVsHealthy_AllSamples.rdata")) # Contains eQTM results
# load(file = file.path(root.dir, "/Data/007 DMR eQTM analysis/toDMR_EQTM_AsthmaVsHealthy_allSamples.rdata")) # Contains eQTM results

# test.df <- finaltable[finaltable$snps %in% top.cpgs, ]
# test.df <- test.df[test.df$finalFDR < 0.05, ]
# write.csv(test.df, "C:/Users/24976197/OneDrive - UTS/Desktop/test_df.csv")

# Set input folder
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/005 eQTM analysis"
input.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/002 Differential Methylation Analysis DMA/"

# write.csv(x = finaltable, file = file.path(base.dir, "EQTM_results_all.csv"))

# dir <- "figures"
# !dir.exists(file.path(input.dir, dir)) && dir.create(file.path(input.dir, dir), recursive = TRUE)

meth.data <- read.csv(file = file.path(base.dir, "methylationData_eQTM_AllSamples_Th2HighVsHealthy.csv"),
                      header = TRUE,
                      row.names = 1)
# meth.data <- read.csv(file = file.path(root.dir, "Data/007 DMR eQTM analysis/DMR_eQTM_mval_AsthmaVsHealthy.csv"),
#                       header = TRUE,
#                       row.names = 1)

source("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Documentation/Github/004 Bioinformatics/Scripts/M2beta.R")
meth.data <- M2beta(meth.data)
expr.data <- read.csv(file = file.path(base.dir, "normalized_geneExpressionData_eQTM_asthma_vs_control_samples.csv"),
                      header = TRUE,
                      row.names = 1)
DMA.result <- read.csv(file = file.path(input.dir, "Correct for PC1-3 - Th2 groups/DMA_Tt2_significant_CpG_Th2HighVsHealthy_BH.csv"),
                       header = TRUE,
                       row.names = 1)
DMA.result <- read.csv(file = file.path(root.dir, "Data/006 Differential Methylated Region Analysis DMR/DMR_results_asthma_vs_control_750gap.csv"),
                       header = TRUE,
                       row.names = 1)
covariates <- read.csv(file = file.path(base.dir, "covariates_eQTM_samples_original.csv"),
                       header = TRUE,
                       row.names = 1)
# covariates <- as.data.frame(t(covariates))

Th2.groups <- read.csv(file = file.path(base.dir, "covariates_eQTM_samples_IDs.csv"),
                       header = TRUE,
                       row.names = 1)

colnames(meth.data) <- gsub(pattern = "X", replacement = "", x = colnames(meth.data))

colnames(expr.data) <- gsub(pattern = "X", replacement = "", x = colnames(expr.data))

# Check the patients that are being used
covariates <- covariates[order(covariates$PT),]
covariates <- Th2.groups
covariates <- covariates[order(covariates$PT),]
covariates <- covariates[covariates$meth_file_id %in% colnames(meth.data), ]

# Th2.groups <- Th2.groups[Th2.groups$meth_file_id %in% colnames(meth.data), ]

# Filter out the patients
meth.data <- meth.data[ ,covariates$meth_file_id]
expr.data <- expr.data[ ,covariates$GenomeScan_ID]

##############################################################################
####################  ###########################
##############################################################################

##############################################################################
# Select the significant CpG sites from the DMA analysis
#   and search what genes they may influence
#   Don't care about eQTM on its own, this way gives more biological relevance
rownames(DMA.result)
# DMA.result <- DMA.result[-which(DMA.result$n==1),]
# DMA.result <- DMA.result[DMA.result$p.adjust < 0.05, ]
# DMA.result$region <- paste0(DMA.result$chr,":", DMA.result$start, "-", DMA.result$end)
# rownames(DMA.result) <- DMA.result$region

# Search and try to get the first 100 significant DMP from the DMA analysis
# Trial and error by increasing the number of rownames
no.CpGs <- nrow(DMA.result)
DMA.finaltable.list.100 <- finaltable[finaltable$snps %in% rownames(DMA.result)[1:no.CpGs], ]

# Order the filtered eQTM list with the order of significance from the DMA result
DMA.finaltable.list.100 <- DMA.finaltable.list.100[match(rownames(DMA.result)[1:no.CpGs], DMA.finaltable.list.100$snps), ]

# If a CpG couldn't be found in the original eQTM list, it gave an NA
#   This is because the methylation site was not nearby (50k)
#     a gene in the DNA. Therefore, removing the entries with an NA
DMA.finaltable.list.100 <- na.omit(DMA.finaltable.list.100)

eQTM.GE.data <- expr.data[rownames(expr.data) %in% DMA.finaltable.list.100$gene,]
####
eQTM.GE.data <- expr.data
####
gene.conversion.table <- read.csv(file = file.path(base.dir, "geneLocation_63k.csv"),
                                  header = TRUE,
                                  row.names = 1,
                                  na.strings = "")
gene.conversion.table <- gene.conversion.table[complete.cases(gene.conversion.table), ]

# Converting the gene expression ensembl ID to gene symbol
conversion_vector <- setNames(gene.conversion.table$hgnc_symbol, gene.conversion.table$ensembl_gene_id)

new_rownames <- ifelse(
  rownames(eQTM.GE.data) %in% names(conversion_vector),
  conversion_vector[rownames(eQTM.GE.data)],
  rownames(eQTM.GE.data)
)

rownames(eQTM.GE.data) <- new_rownames

# Converting the ensembl ID from the eQTM finaltable list
# Initialize the symbol.vector with the same length as the genes in DMA.finaltable.list.100
symbol.vector <- character(length(DMA.finaltable.list.100$gene))
multiple_matches <- list()

# Loop through each gene to find and assign the corresponding symbol
for (i in seq_along(DMA.finaltable.list.100$gene)) {
  gene <- DMA.finaltable.list.100$gene[i]
  found_index <- grep(pattern = gene, x = gene.conversion.table$ensembl_gene_id)
  
  if (length(found_index) == 1) {
    # If exactly one match is found, use the symbol
    symbol <- gene.conversion.table[found_index, "hgnc_symbol"]
    # If the symbol is NA, use the Ensembl gene ID instead
    if (is.na(symbol)) {
      symbol <- gene
      cat("This ENSEMBL ID has no gene symbol\n")
    }
  } else if (length(found_index) > 1) {
    # Save entries with multiple matches for later review
    multiple_matches[[gene]] <- gene.conversion.table[found_index, "hgnc_symbol"]
    # Retain the original gene ID for now
    symbol <- gene
  } else {
    # If no match is found, retain the original gene ID
    symbol <- gene
  }
  
  # Assign the symbol (or gene ID) to the corresponding position in symbol.vector
  symbol.vector[i] <- symbol
}

DMA.finaltable.list.100$symbol <- symbol.vector

DMA.finaltable.list.100.sign <- DMA.finaltable.list.100[DMA.finaltable.list.100$finalFDR < 0.05, ]

finaltable_sig_20 <- DMA.finaltable.list.100.sign

# Initialize data frames to store valid and invalid pairs
valid_pairs <- data.frame(
  CpG = character(),
  Gene = character(),
  Mean_Th2_Low = numeric(),
  Mean_Healthy = numeric(),
  Mean_Th2_High = numeric(),
  Mean_Undetermined = numeric(),
  Condition = character(),
  stringsAsFactors = FALSE
)

invalid_pairs <- data.frame(
  CpG = character(),
  Gene = character(),
  Mean_Th2_Low = numeric(),
  Mean_Healthy = numeric(),
  Mean_Th2_High = numeric(),
  Mean_Undetermined = numeric(),
  Reason = character(),
  stringsAsFactors = FALSE
)

# Loop through each CpG-Gene pair and perform the checks
for (i in 1:nrow(finaltable_sig_20)) {
  cpg.name <- finaltable_sig_20[i, "snps"]
  gene.name <- finaltable_sig_20[i, "symbol"]
  
  # Extract the methylation and gene expression data
  df <- t(meth.data[cpg.name, ])
  df <- cbind(df, t(eQTM.GE.data[gene.name, ]))
  df <- as.data.frame(df)
  colnames(df) <- c("cpg", "gene")
  
  # Map groups (asthma status)
  df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "group_th"]
  df$asthma <- gsub(pattern = "healthy", replacement = "Healthy", x = df$asthma)
  df$asthma <- gsub(pattern = "high", replacement = "Th2-High", x = df$asthma)
  df$asthma <- gsub(pattern = "low", replacement = "Th2-Low", x = df$asthma)
  df$asthma <- gsub(pattern = "undeterm", replacement = "Undetermined subtype", x = df$asthma)
  df$asthma <- factor(df$asthma, levels = c("Healthy", "Th2-Low", "Undetermined subtype", "Th2-High"))
  
  # Calculate the mean values for each group, including "Undetermined subtype"
  mean_cpg <- tapply(df$cpg, df$asthma, mean, na.rm = TRUE)
  
  # Perform the conditions check, focusing on Th2-High first
  if (!is.na(mean_cpg["Th2-High"]) && !is.na(mean_cpg["Healthy"]) && !is.na(mean_cpg["Th2-Low"])) {
    if (mean_cpg["Th2-High"] > mean_cpg["Healthy"] && mean_cpg["Th2-Low"] < mean_cpg["Healthy"]) {
      condition <- "Th2-High > Healthy, Th2-Low < Healthy"
      # Add to valid_pairs if the condition is met
      valid_pairs <- rbind(valid_pairs, data.frame(
        CpG = cpg.name,
        Gene = gene.name,
        Mean_Th2_Low = mean_cpg["Th2-Low"],
        Mean_Healthy = mean_cpg["Healthy"],
        Mean_Th2_High = mean_cpg["Th2-High"],
        Mean_Undetermined = mean_cpg["Undetermined subtype"],
        Condition = condition,
        stringsAsFactors = FALSE
      ))
    } else if (mean_cpg["Th2-High"] < mean_cpg["Healthy"] && mean_cpg["Th2-Low"] > mean_cpg["Healthy"]) {
      condition <- "Th2-High < Healthy, Th2-Low > Healthy"
      # Add to valid_pairs if the condition is met
      valid_pairs <- rbind(valid_pairs, data.frame(
        CpG = cpg.name,
        Gene = gene.name,
        Mean_Th2_Low = mean_cpg["Th2-Low"],
        Mean_Healthy = mean_cpg["Healthy"],
        Mean_Th2_High = mean_cpg["Th2-High"],
        Mean_Undetermined = mean_cpg["Undetermined subtype"],
        Condition = condition,
        stringsAsFactors = FALSE
      ))
    } else {
      # Condition failed: Capture the reason
      reason <- ifelse(mean_cpg["Th2-High"] > mean_cpg["Healthy"], 
                       "Th2-Low not less than Healthy", 
                       "Th2-Low not greater than Healthy")
      invalid_pairs <- rbind(invalid_pairs, data.frame(
        CpG = cpg.name,
        Gene = gene.name,
        Mean_Th2_Low = mean_cpg["Th2-Low"],
        Mean_Healthy = mean_cpg["Healthy"],
        Mean_Th2_High = mean_cpg["Th2-High"],
        Mean_Undetermined = mean_cpg["Undetermined subtype"],
        Reason = reason,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    # Capture the reason for missing data in invalid_pairs
    reason <- "Missing data for one or more groups"
    invalid_pairs <- rbind(invalid_pairs, data.frame(
      CpG = cpg.name,
      Gene = gene.name,
      Mean_Th2_Low = mean_cpg["Th2-Low"],
      Mean_Healthy = mean_cpg["Healthy"],
      Mean_Th2_High = mean_cpg["Th2-High"],
      Mean_Undetermined = mean_cpg["Undetermined subtype"],
      Reason = reason,
      stringsAsFactors = FALSE
    ))
  }
}

# Filter the valid_pairs CpG, meaning that the Th2-high and Th2-low 
#   going in different direction of methylation
DMA.result.filt <- DMA.result[valid_pairs$CpG, ]
# Re-order on significance 
DMA.result.filt <- DMA.result.filt[order(DMA.result.filt$adj.P.Val), ]
# Change rownames so we order the valid-pairs on significance too
rownames(valid_pairs) <- valid_pairs$CpG
valid_pairs <- valid_pairs[rownames(DMA.result.filt), ]

# Reorder the eQTM final result table of only valid-pairs and on significance
#   This is needed for presentation slide/table
rownames(finaltable_sig_20) <- finaltable_sig_20$snps
finaltable_sig_20 <- finaltable_sig_20[rownames(valid_pairs), ]
write.table(finaltable_sig_20, file = file.path(root.dir, "/Results/005 eQTM analysis/Th2-High_unique_eQTM.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")


COLOR_TEMP = c("#d5896f","#dab785","#70a288", "#6A9FB5")

# Asthma  Healthy
COLOR_TEMP = c("#ff7f0e","#1f77b4")
# Healthy   Th2-high
COLOR_TEMP = c("#1f77b4", "#e6550d")
# Healthy   Th2-Low
COLOR_TEMP = c("#1f77b4", "#ffb74d")
# Th2-low   Th2-high 
COLOR_TEMP = c("#ffb74d", "#e6550d")
# Healthy   Th2-high  Th2-low   Undetermined
COLOR_TEMP = c("#1f77b4", "#ffb74d", "black", "#e6550d")

# Functions adapted from https://stackoverflow.com/a/39611375
#   used for the nearest decimal to edit the axes limits
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

# pdf(file = file.path(root.dir, "/Results/005 eQTM analysis/TOP_eQTM_Th2HighvsHealthy_NEW.pdf"),
#     width = 10, height = 8.5)
# pdf(file = file.path(root.dir, "/Results/005 eQTM analysis/Most_important_CpG_feature_ROC.pdf"),
#     width = 10, height = 8.5)

# DMPs
valid_pairs2 <- data.frame(
  CpG = c("cg21501207", "cg02333649", "cg01482588", "cg24224501", "cg19398575", 
           "cg25710507", "cg20505547", "cg16911809", "cg01482377", "cg26554722"),
  Gene = c("SH2D1B", "MRPL40", "NTRK1", "LRRC17", "SAMSN1", 
                  "NTRK2", "DENND5A", "CCL26", "NTRK1", "TLE3")
)

# DMRs
# valid_pairs2 <- data.frame(
#   CpG = c("chr17:28717304-28718284",
#           "chr7:79452350-79454722",
#           "chr6:43645019-43645658",
#           "chr3:6861137-6861640",
#           "chr1:156862644-156864002",
#           "chr2:191845500-191848140",
#           "chr11:57646984-57649073",
#           "chr8:94949850-94950235",
#           "chr9:84672079-84672092"),
#   Gene = c("RPL23A",
#            "MAGI2",
#            "GTPBP2",
#            "GRM7",
#            "NTRK1",
#            "CAVIN2-AS1",
#            "SERPING1",
#            "NDUFAF6",
#            "NTRK2"
#            )
# )


# Create a list to store the individual plots
plot_list <- list()

for (i in 1:nrow(valid_pairs2)) {
  print(i)
  
  cpg.name <- valid_pairs2[i, "CpG"]
  gene.name <- valid_pairs2[i, "Gene"]
   
  # cpg.name <- "chr22:41367089-41367508"
  # gene.name <- "ENSG00000167074"
  
  # Extract the methylation and gene expression data
  df <- t(meth.data[cpg.name, ])
  df <- cbind(df, t(eQTM.GE.data[gene.name, ]))
  df <- as.data.frame(df)
  colnames(df) <- c("cpg", "gene")
  
  incl.endotypes <- FALSE
  
  # Map groups (asthma status)
  if (incl.endotypes == FALSE) {
    COLOR_TEMP = c("#1f77b4", "#ff7f0e")
    # Process using 'ASTHEA'
    df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "ASTHEA"]
    df$asthma <- gsub(pattern = "H", replacement = "Healthy", x = df$asthma)
    df$asthma <- gsub(pattern = "A", replacement = "Asthma", x = df$asthma)
    df$asthma <- factor(df$asthma, levels = c("Healthy", "Asthma"))
    
  } else {
    COLOR_TEMP = c("#1f77b4", "#ffb74d", "black", "#e6550d")
    # Process using 'group_th'
    df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "group_th"]
    df$asthma <- gsub(pattern = "healthy", replacement = "Healthy", x = df$asthma)
    df$asthma <- gsub(pattern = "high", replacement = "Th2-High", x = df$asthma)
    df$asthma <- gsub(pattern = "low", replacement = "Th2-Low", x = df$asthma)
    df$asthma <- gsub(pattern = "undeterm", replacement = "Undetermined subtype", x = df$asthma)
    df$asthma <- factor(df$asthma, levels = c("Healthy", "Th2-Low", "Undetermined subtype", "Th2-High"))
  }
  
  max.cpg.limit <- ceiling_dec(max(df$cpg), 1)
  min.cpg.limit <- floor_dec(min(df$cpg), 1)
  max.gene.limit <- ceiling_dec(max(df$gene), 1)
  min.gene.limit <- floor_dec(min(df$gene), 1)
  
  # mean_cpg_rounded <- round(c(
  #   Healthy = valid_pairs$Mean_Healthy[i],
  #   `Th2-Low` = valid_pairs$Mean_Th2_Low[i],
  #   `Undetermined subtype` = valid_pairs$Mean_Undetermined[i],
  #   `Th2-High` = valid_pairs$Mean_Th2_High[i]
  #   ), digits = 3)
  
  # eQTM plot
  ggplot(data = df) +
    aes(y = cpg, 
        x = gene) +
    geom_point(aes(colour = asthma),
               size = 1.3,
               alpha = 0.9,
               stroke = 0.7) +
    geom_smooth(linetype = "dashed",
                se = F,
                color = 'black',
                formula = 'y ~ x',
                method = "lm") +
    # geom_smooth(aes(colour = asthma),
    #             formula = 'y ~ x',
    #             method = "lm") +
    ## Stats within the correlation plot
    stat_cor(r.accuracy = 0.01, p.accuracy = 0.001,
             label.x.npc = 0.75, label.y.npc = 0.88, size = 5) +
    stat_regline_equation(label.x.npc = 0.75, label.y.npc = 0.95, size = 5) +
    
    scale_y_continuous(limits = c(min.cpg.limit,max.cpg.limit), breaks = pretty_breaks(5)) +
    scale_x_continuous(limits = c(min.gene.limit, max.gene.limit), breaks = pretty_breaks(5)) +
    # scale_fill_manual(values = COLOR_TEMP) +
    scale_color_manual(values = COLOR_TEMP) +
    guides(colour = guide_legend(title="Status")) +
    labs(x = paste(gene.name, "expression [normalized]"),
         y = paste(cpg.name, "methylation [Beta]")
    ) +
    theme_bw() + theme(legend.position = "none",
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       axis.text.x = element_text(size = 14),
                       axis.text.y = element_text(size = 14)) -> p_eQTM
  # facet_wrap(~ asthma)-> p_eQTM
  
  # CpG plot
  ggplot(data = df) +
    aes(y = cpg, 
        x = asthma,
        fill = asthma) + #
    geom_flat_violin(position = position_nudge(x = .20), 
                     alpha = .8) +
    geom_point(aes(color = asthma), 
               position = position_jitter(width = .10),
               size = .3, 
               alpha = .5,
               show.legend = F) +
    
    geom_boxplot(width = .3, 
                 outlier.shape = NA,
                 alpha = .5) +
    labs(x = "Status") +
    scale_fill_manual(values = COLOR_TEMP, guide = NULL) +
    scale_color_manual(values = COLOR_TEMP) +
    scale_y_continuous(limits = c(min.cpg.limit, max.cpg.limit), breaks = pretty_breaks(5)) +
    # expand_limits(x = 3, y = 10) + 
    theme_bw() + 
    theme(axis.ticks.x=element_blank(), 
          axis.title.x = element_blank(), 
          axis.text.x = element_blank(), 
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) -> p_CpG
    # annotate("text", x = 1, y = min.cpg.limit,
    #          label = paste("Healthy:", mean_cpg_rounded["Healthy"]),
    #          size = 3, hjust = 0.5) +
    # annotate("text", x = 2, y = min.cpg.limit,
    #          label = paste("Th2-Low:", mean_cpg_rounded["Th2-Low"]),
    #          size = 3, hjust = 0.5) +
    # annotate("text", x = 3, y = min.cpg.limit,
    #          label = paste("Th2-Intermediate:", mean_cpg_rounded["Undetermined subtype"]),
    #          size = 3, hjust = 0.5) +
    # annotate("text", x = 4, y = min.cpg.limit,
    #          label = paste("Th2-High:", mean_cpg_rounded["Th2-High"]),
    #          size = 3, hjust = 0.5) -> p_CpG
  
  #gene expression plot
  ggplot(data = df) +
    aes(y = gene, 
        x = asthma,
        fill = asthma) + #
    geom_flat_violin(position = position_nudge(x = .2), 
                     alpha = .8) +
    geom_point(aes(color = asthma), 
               position = position_jitter(width = .15),
               size = .3, 
               alpha = .5,
               show.legend = F) +
    geom_boxplot(width = .3, 
                 outlier.shape = NA,
                 alpha = .5) +
    # coord_flip() +
    labs(x = "Status", fill = NULL) +
    scale_fill_manual(values = COLOR_TEMP, guide = NULL) +
    scale_color_manual(values = COLOR_TEMP) +
    scale_y_continuous(limits = c(min.gene.limit, max.gene.limit), breaks = pretty_breaks(5)) +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y=element_blank(),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
    ) + coord_flip() -> p_gene
  
  layout = "222222222##
            222222222##
            33333333311
            33333333311
            33333333311
            33333333311
            33333333311
            33333333311
            33333333311"
  
  # print(p_CpG)
  # print(p_gene)
  print(p_eQTM)
  
  p_all <- p_CpG + p_gene + p_eQTM + 
    plot_layout(design = layout)
  
  
  # print(p_all)
  # Add combined plot to list
  plot_list[[i]] <- p_eQTM
  
  # print(p_all)
}

combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 3)

tiff(file.path("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS",
               "/Results/005 eQTM analysis/top_sig_CpGs_Th2_High_Healthy_FINAL.tiff"), width = 3000, height = 3000, res = 300)
# tiff(file.path("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS",
#                "/Results/007 DMR eQTM analysis/top_sig_DMRs_Asthma_Healthy_FINAL.tiff"), width = 3000, height = 3000, res = 300)
print(combined_plot)
dev.off()
##############################################################################

# 
# 
# # Select the significant eQTM
# finaltable_sig2 <- finaltable_sig[finaltable_sig$gene %in% rownames(expr.data), ]
# # Select the significant CpG sites from the DMA analysis
# finaltable_sig2 <- rbind(finaltable_sig2, finaltable[finaltable$snps %in% rownames(DMA.result), ])
# finaltable_sig2 <- finaltable_sig2[!duplicated(finaltable_sig2), ]
# 
# eQTM.GE.data <- expr.data[rownames(expr.data) %in% finaltable$gene,]
# 
# finaltable_sig_20 <- head(finaltable_sig2, n =20)
# 
# 
# gene.conversion.table <- read.csv(file = file.path(base.dir, "geneLocation_63k.csv"),
#                                   header = TRUE,
#                                   row.names = 1,
#                                   na.strings = "")
# 
# symbol.vector <- c()
# for (gene in rownames(eQTM.GE.data)) {
#   found_index <- grep(pattern = gene, x = gene.conversion.table$ensembl_gene_id)
#   symbol <- gene.conversion.table[found_index, "hgnc_symbol"]
#   if(is.na(symbol)) {
#     print("This ENSEMBL ID has no gene symbol")
#     symbol <- gene.conversion.table[found_index, "ensembl_gene_id"]
#   }
#   symbol.vector <- c(symbol.vector, symbol)
# }
# rownames(eQTM.GE.data) <- symbol.vector
# 
# symbol.vector <- c()
# for (gene in finaltable_sig_20$gene) {
#   found_index <- grep(pattern = gene, x = gene.conversion.table$ensembl_gene_id)
#   symbol <- gene.conversion.table[found_index, "hgnc_symbol"]
#   if(is.na(symbol)) {
#     print("This ENSEMBL ID has no gene symbol")
#     symbol <- gene.conversion.table[found_index, "ensembl_gene_id"]
#   }
#   symbol.vector <- c(symbol.vector, symbol)
# }
# finaltable_sig_20$gene <- symbol.vector
# 
# 
# all.DMR.eQTMS <- finaltable[order(finaltable$finalFDR), ]
# symbol.vector <- c()
# for (gene in all.DMR.eQTMS$gene) {
#   found_index <- grep(pattern = gene, x = gene.conversion.table$ensembl_gene_id)
#   if(length(found_index) == 0) {next}
#   if(length(found_index) > 1) {next}
#   symbol <- gene.conversion.table[found_index, "hgnc_symbol"]
#   if(is.na(symbol)) {
#     cat("This ENSEMBL ID has no gene symbol\n")
#     symbol <- gene.conversion.table[found_index, "ensembl_gene_id"]
#   }
#   symbol.vector <- c(symbol.vector, symbol)
# }
# all.DMR.eQTMS$gene <- symbol.vector
# finaltable_sig_20 <- all.DMR.eQTMS
# 
# finaltable_sig_20 <- finaltable_sig_20[finaltable_sig_20$gene == "ENSG00000197558", ]
# finaltable_sig_20$gene <- "SSPOP"
# finaltable_sig_20 <- finaltable_sig_20[grep(pattern = "0655", x = finaltable_sig_20$snps), ]
# 
# 
# 
# 
# View(finaltable_sig_20)
# # cg16964439_BC21 ENSG00000148053
# # 
# df <- data.frame()
# df <- t(meth.data["cg24881657_BC11", ])
# df <- cbind(df, eQTM.GE.data["ENSG00000105639", ])
# colnames(df) <- c("cpg", "gene")
# df <- as.data.frame(df)
# df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "ASTHEA"]
# 
# ggplot(data = df,
#        aes(x = gene,
#            y = cpg)) +
#   geom_point() +
#   facet_wrap(~ asthma)
# # 
# # 
# # pdf(file = file.path(root.dir, "/Results/005 eQTM analysis/Top20eQTM.pdf"))
# # 
# # for (i in 1:nrow(finaltable_sig_20)) {
# #   df <- data.frame()
# #   df <- t(meth.data[i, ])
# #   df <- cbind(df, eQTM.GE.data[i, ])
# #   df <- as.data.frame(df)
# #   colnames(df) <- c("cpg", "gene")
# #   df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "ASTHEA"]
# #   df$asthma <- gsub(pattern = "0", replacement = "Healthy", x = df$asthma)
# #   df$asthma <- gsub(pattern = "1", replacement = "Asthma", x = df$asthma)
# #   
# #   plot <- ggplot(data = df,
# #          aes(x = gene,
# #              y = cpg,
# #              colour = asthma)) +
# #     scale_colour_manual(values = c("red","green3")) +
# #     # expand_limits(x = 3, y = 10) +
# # 
# #     # scale_x_continuous(limits = c(-3, 10),
# #     #                    breaks = seq(-3, 10, 1)) +
# #     # scale_y_continuous(limits = c(0,1),
# #     #                    breaks = seq(0,1, 0.1)) +
# #     geom_point() +
# #     facet_wrap(~ asthma) +
# #     geom_smooth(method='lm', formula= y~x) +
# #     labs(x = rownames(eQTM.GE.data)[i],
# #        y = rownames(meth.data[i, ]),
# #        colour = "Status",
# #        title = "") +
# #     guides(fill = guide_legend(nrow=1,
# #                                byrow=TRUE))+
# #     theme_bw() +
# #     theme(legend.position = "bottom",
# #           legend.margin = margin(-5, 0, 0, 0))
# #   
# #   plot2 <- ggpubr::ggscatter(df, x = "gene", y = "cpg", size = 0.3, 
# #                        rug = TRUE,                                # Add marginal rug
# #                        color = "asthma", palette = "jco") +
# #     ggpubr::stat_cor(aes(color = asthma), method = "spearman") + 
# #     xlab(label = rownames(eQTM.GE.data)[i]) + ylab(label = rownames(meth.data[i, ]))
# #   plot(plot2)
# #   # print(plot)
# # }
# # dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# COLOR_TEMP = c("#d5896f","#dab785","#70a288", "#6A9FB5")
# 
# # Asthma  Healthy
# COLOR_TEMP = c("#ff7f0e","#1f77b4")
# # Healthy   Th2-high
# COLOR_TEMP = c("#1f77b4", "#e6550d")
# # Healthy   Th2-Low
# COLOR_TEMP = c("#1f77b4", "#ffb74d")
# # Th2-low   Th2-high 
# COLOR_TEMP = c("#ffb74d", "#e6550d")
# # Healthy   Th2-high  Th2-low   Undetermined
# COLOR_TEMP = c("#1f77b4", "#ffb74d", "black", "#e6550d")
# 
# # Functions adapted from https://stackoverflow.com/a/39611375
# #   used for the nearest decimal to edit the axes limits
# floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
# ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
# 
# 
# # Create progression bar for for-loop
# # n_iter <- nrow(finaltable_sig_20)
# # progress.bar <- txtProgressBar(min = 0,
# #                                initial = 0, 
# #                                max = n_iter,
# #                                style = 3,
# #                                width = n_iter, # Needed to avoid multiple printings
# #                                char = "=") 
# # init <- numeric(n_iter)
# # end <- numeric(n_iter)
# 
# pdf(file = file.path(root.dir, "/Results/005 eQTM analysis/TOP_eQTM_Th2HighvsHealthy_NEW.pdf"),
#     width = 10, height = 8.5)
# for (i in 1:nrow(finaltable_sig_20)) {
#   print(i)
#   df <- data.frame()
#   cpg.name <- finaltable_sig_20[i, "snps"]
#   gene.name <- finaltable_sig_20[i, "symbol"]
#   df <- t(meth.data[cpg.name, ])
#   df <- cbind(df, t(eQTM.GE.data[gene.name, ]))
#   df <- as.data.frame(df)
#   colnames(df) <- c("cpg", "gene")
#   
#   ## Asthma vs Healthy
#   # df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "ASTHEA"]
#   # df$asthma <- gsub(pattern = "H", replacement = "Healthy", x = df$asthma)
#   # df$asthma <- gsub(pattern = "A", replacement = "Asthma", x = df$asthma)
#   
#   ## Th2-High vs Healthy | Th2-Low vs Healthy | Th2-High vs Th2-Low
#   df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "group_th"]
#   df$asthma <- gsub(pattern = "healthy", replacement = "Healthy", x = df$asthma)
#   df$asthma <- gsub(pattern = "high", replacement = "Th2-High", x = df$asthma)
#   df$asthma <- gsub(pattern = "low", replacement = "Th2-Low", x = df$asthma)
#   df$asthma <- gsub(pattern = "undeterm", replacement = "Undetermined subtype", x = df$asthma)
#   df$asthma <- factor(df$asthma, levels = c("Healthy", "Th2-Low", "Undetermined subtype", "Th2-High"))
#   
#   mean_cpg <- tapply(df$cpg, df$asthma, mean, na.rm = TRUE)
# 
#   print(mean_cpg)
#   if (!is.na(mean_cpg["Th2-Low"]) && !is.na(mean_cpg["Healthy"]) && !is.na(mean_cpg["Th2-High"])) {
#     if (mean_cpg["Th2-Low"] > mean_cpg["Healthy"] && mean_cpg["Th2-High"] < mean_cpg["Healthy"]) {
#       print("Condition passed: Th2-Low > Healthy, Th2-High < Healthy")
#     } else if (mean_cpg["Th2-Low"] < mean_cpg["Healthy"] && mean_cpg["Th2-High"] > mean_cpg["Healthy"]) {
#       print("Condition passed: Th2-Low < Healthy, Th2-High > Healthy")
#     } else {
#       print("Condition failed, skipping this CpG-Gene pair")
#       next
#     }
#   } else {
#     print("Not enough data for all groups, skipping this CpG-Gene pair")
#     next
#   }
#   
#   mean_cpg_rounded <- round(mean_cpg, 3)
#   
#   
#   max.cpg.limit <- ceiling_dec(max(df$cpg), 1)
#   min.cpg.limit <- floor_dec(min(df$cpg), 1)
#   max.gene.limit <- ceiling_dec(max(df$gene), 1)
#   min.gene.limit <- floor_dec(min(df$gene), 1)
#   
#   # eQTM plot
#   ggplot(data = df) +
#     aes(y = cpg, 
#         x = gene) +
#     geom_point(aes(colour = asthma),
#                size = 1.3,
#                alpha = 0.9,
#                stroke = 0.7) +
#     geom_smooth(linetype = "dashed",
#                 se = F,
#                 color = 'black',
#                 formula = 'y ~ x',
#                 method = "lm") +
#     # geom_smooth(aes(colour = asthma),
#     #             formula = 'y ~ x',
#     #             method = "lm") +
#     stat_cor(r.accuracy = 0.01, p.accuracy = 0.001, 
#              label.x.npc = 0.75, label.y.npc = 0.88, size = 5) +
#     stat_regline_equation(label.x.npc = 0.75, label.y.npc = 0.95, size = 5) + 
#     
#     scale_y_continuous(limits = c(min.cpg.limit,max.cpg.limit), breaks = pretty_breaks(5)) +
#     scale_x_continuous(limits = c(min.gene.limit, max.gene.limit), breaks = pretty_breaks(5)) +
#     # scale_fill_manual(values = COLOR_TEMP) +
#     scale_color_manual(values = COLOR_TEMP) +
#     guides(colour = guide_legend(title="Status")) +
#     labs(x = paste(gene.name, "expression [normalized]"),
#          y = paste(cpg.name, "methylation [Beta]")
#     ) +
#     theme_bw() + theme(legend.position = "none",
#                        panel.grid.major = element_blank(),
#                        panel.grid.minor = element_blank()) -> p_eQTM
#   # facet_wrap(~ asthma)-> p_eQTM
#   
#   # CpG plot
#   ggplot(data = df) +
#     aes(y = cpg, 
#         x = asthma,
#         fill = asthma) + #
#     geom_flat_violin(position = position_nudge(x = .20), 
#                      alpha = .8) +
#     geom_point(aes(color = asthma), 
#                position = position_jitter(width = .10),
#                size = .3, 
#                alpha = .5,
#                show.legend = F) +
#     
#     geom_boxplot(width = .3, 
#                  outlier.shape = NA,
#                  alpha = .5) +
#     labs(x = "Status") +
#     scale_fill_manual(values = COLOR_TEMP, guide = NULL) +
#     scale_color_manual(values = COLOR_TEMP) +
#     scale_y_continuous(limits = c(min.cpg.limit, max.cpg.limit), breaks = pretty_breaks(5)) +
#     # expand_limits(x = 3, y = 10) + 
#     theme_bw() + 
#     theme(axis.ticks.x=element_blank(), 
#           axis.title.x = element_blank(), 
#           # axis.title.y=element_blank(),
#           # axis.text.y=element_blank(), 
#           # axis.ticks.y=element_blank(),
#           # panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()) +
#     annotate("text", x = 1, y = min.cpg.limit, 
#              label = paste("Healthy:", mean_cpg_rounded["Healthy"]),
#              size = 3, hjust = 0.5) +
#     annotate("text", x = 2, y = min.cpg.limit, 
#              label = paste("Th2-Low:", mean_cpg_rounded["Th2-Low"]),
#              size = 3, hjust = 0.5) +
#     annotate("text", x = 3, y = min.cpg.limit, 
#              label = paste("Th2-Intermediate:", mean_cpg_rounded["Undetermined subtype"]),
#              size = 3, hjust = 0.5) +
#     annotate("text", x = 4, y = min.cpg.limit, 
#              label = paste("Th2-High:", mean_cpg_rounded["Th2-High"]),
#              size = 3, hjust = 0.5) -> p_CpG
#   
#   #gene expression plot
#   ggplot(data = df) +
#     aes(y = gene, 
#         x = asthma,
#         fill = asthma) + #
#     geom_flat_violin(position = position_nudge(x = .2), 
#                      alpha = .8) +
#     geom_point(aes(color = asthma), 
#                position = position_jitter(width = .15),
#                size = .3, 
#                alpha = .5,
#                show.legend = F) +
#     geom_boxplot(width = .3, 
#                  outlier.shape = NA,
#                  alpha = .5) +
#     # coord_flip() +
#     labs(x = "Status", fill = NULL) +
#     scale_fill_manual(values = COLOR_TEMP, guide = NULL) +
#     scale_color_manual(values = COLOR_TEMP) +
#     scale_y_continuous(limits = c(min.gene.limit, max.gene.limit), breaks = pretty_breaks(5)) +
#     theme_bw() +
#     theme(axis.title.x=element_blank(),
#           # axis.text.x=element_blank(),
#           axis.ticks.x=element_blank(), 
#           # axis.title.y = element_blank(), 
#           # axis.ticks.y=element_blank(),
#           # panel.border = element_blank(),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank()
#     ) -> p_gene
#   
#   
#   # layout = "11333333333
#   #           11333333333
#   #           11333333333
#   #           11333333333
#   #           11333333333   
#   #           11333333333
#   #           11333333333
#   #           ##222222222
#   #           ##222222222"
#   
#   
#   layout = "222222222##
#             222222222##
#             33333333311
#             33333333311
#             33333333311
#             33333333311
#             33333333311
#             33333333311
#             33333333311"
#   
#   # print(p_CpG)
#   # print(p_gene)
#   # print(p_eQTM)
#   
#   p_all <- p_CpG + p_gene + p_eQTM + 
#     plot_layout(design = layout) 
#   # + 
#   #   labs(title = "Gene expression vs. beta-values between asthma and healthy")
#   
#   print(p_all)
#   
# }
# dev.off()
# 
# 
# 
# 
# 
# # 
# # for (i in 1:nrow(finaltable_sig_20)) {
# #   # Initiate tracking process of loop
# #   init[i] <- Sys.time()
# #   #---------------------
# #   
# #   # Data manipulation to create figures
# #   df <- data.frame()
# #   cpg.name <- finaltable_sig_20[i, "snps"]
# #   gene.name <- finaltable_sig_20[i, "gene"]
# #   df <- t(meth.data[cpg.name, ])
# #   df <- cbind(df, t(eQTM.GE.data[gene.name, ]))
# #   df <- as.data.frame(df)
# #   colnames(df) <- c("cpg", "gene")
# #   
# #   ## Determine colours for the comparison - Asthma vs Healthy
# #   # df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "ASTHEA"]
# #   # df$asthma <- gsub(pattern = "0", replacement = "Healthy", x = df$asthma)
# #   # df$asthma <- gsub(pattern = "1", replacement = "Asthma", x = df$asthma)
# #   
# #   ## Determine colours for the comparison - Th2-High vs Healthy | Th2-Low vs Healthy | Th2-High vs Th2-Low
# #   df$asthma <- covariates[covariates$meth_file_id %in% rownames(df), "group_th"]
# #   df$asthma <- gsub(pattern = "healthy", replacement = "Healthy", x = df$asthma)
# #   df$asthma <- gsub(pattern = "high", replacement = "Th2-High", x = df$asthma)
# #   df$asthma <- gsub(pattern = "low", replacement = "Th2-Low", x = df$asthma)
# #   df$asthma <- gsub(pattern = "undeterm", replacement = "Undetermined subtype", x = df$asthma)
# #   
# #   # Setting the axes limits top and bottom
# #   max.cpg.limit <- ceiling_dec(max(df$cpg), 1)
# #   min.cpg.limit <- floor_dec(min(df$cpg), 1)
# #   max.gene.limit <- ceiling_dec(max(df$gene), 1)
# #   min.gene.limit <- floor_dec(min(df$gene), 1)
# #   
# #   # eQTM plot
# #   ggplot(data = df) +
# #     aes(y = cpg, 
# #         x = gene) +
# #     geom_point(aes(colour = asthma),
# #                size = 1.3,
# #                alpha = 0.9,
# #                stroke = 0.7) +
# #     geom_smooth(linetype = "dashed",
# #                 se = F,
# #                 color = 'black',
# #                 formula = 'y ~ x',
# #                 method = "lm") +
# #     # geom_smooth(aes(colour = asthma),
# #     #             formula = 'y ~ x',
# #     #             method = "lm") +
# #     stat_cor(r.accuracy = 0.01, p.accuracy = 0.001, 
# #              label.x.npc = 0.75, label.y.npc = 0.88, size = 5) +
# #     stat_regline_equation(label.x.npc = 0.75, label.y.npc = 0.95, size = 5) + 
# #     
# #     scale_y_continuous(limits = c(min.cpg.limit,max.cpg.limit), breaks = pretty_breaks(5)) +
# #     scale_x_continuous(limits = c(min.gene.limit, max.gene.limit), breaks = pretty_breaks(5)) +
# #     # scale_fill_manual(values = COLOR_TEMP) +
# #     scale_color_manual(values = COLOR_TEMP) +
# #     guides(colour = guide_legend(title="Status")) +
# #     labs(x = paste(gene.name, "expression [normalized]"),
# #          y = paste(cpg.name, "methylation [Beta]")
# #          ) +
# #     theme_bw() + theme(legend.position = "none") -> p_eQTM
# #     # facet_wrap(~ asthma)-> p_eQTM
# # 
# #   # CpG plot
# #   ggplot(data = df) +
# #     aes(y = cpg, 
# #         x = asthma,
# #         fill = asthma) + #
# #     geom_flat_violin(position = position_nudge(x = .20), 
# #                      alpha = .8) +
# #     geom_point(aes(color = asthma), 
# #                position = position_jitter(width = .10),
# #                size = .3, 
# #                alpha = .5,
# #                show.legend = F) +
# #     
# #     geom_boxplot(width = .3, 
# #                  outlier.shape = NA,
# #                  alpha = .5) +
# #     labs(x = "Status") +
# #     scale_fill_manual(values = COLOR_TEMP, guide = NULL) +
# #     scale_color_manual(values = COLOR_TEMP) +
# #     scale_y_continuous(limits = c(min.cpg.limit, max.cpg.limit), breaks = pretty_breaks(5)) +
# #     # expand_limits(x = 3, y = 10) + 
# #     theme_bw() + 
# #     theme(axis.ticks.x=element_blank(), 
# #           axis.title.x = element_blank(), 
# #           axis.title.y=element_blank(),
# #           axis.text.y=element_blank(), 
# #           axis.ticks.y=element_blank(),
# #           panel.border = element_blank(),
# #           panel.grid.major = element_blank(),
# #           panel.grid.minor = element_blank()
# #           ) -> p_CpG
# # 
# #   #gene expression plot
# #   ggplot(data = df) +
# #     aes(y = gene, 
# #         x = asthma,
# #         fill = asthma) + #
# #     geom_flat_violin(position = position_nudge(x = .2), 
# #                      alpha = .8) +
# #     geom_point(aes(color = asthma), 
# #                position = position_jitter(width = .15),
# #                size = .3, 
# #                alpha = .5,
# #                show.legend = F) +
# #     geom_boxplot(width = .3, 
# #                  outlier.shape = NA,
# #                  alpha = .5) +
# #     coord_flip() +
# #     labs(x = "Status", fill = NULL) +
# #     scale_fill_manual(values = COLOR_TEMP, guide = NULL) +
# #     scale_color_manual(values = COLOR_TEMP) +
# #     scale_y_continuous(limits = c(min.gene.limit, max.gene.limit), breaks = pretty_breaks(5)) +
# #     theme_bw() +
# #     theme(axis.title.x=element_blank(),
# #           axis.text.x=element_blank(),
# #           axis.ticks.x=element_blank(), 
# #           axis.title.y = element_blank(), 
# #           axis.ticks.y=element_blank(),
# #           panel.border = element_blank(),
# #           panel.grid.major = element_blank(),
# #           panel.grid.minor = element_blank()
# #           ) -> p_gene
# # 
# # 
# #   layout = "222222222##
# #             222222222##
# #             33333333311
# #             33333333311
# #             33333333311
# #             33333333311
# #             33333333311
# #             33333333311
# #             33333333311"
# #   
# #   
# #   p_all <- p_CpG + p_gene + p_eQTM + 
# #     plot_layout(design = layout) 
# #   
# #   print(p_all)
# # 
# #   #---------------------
# #   
# #   end[i] <- Sys.time()
# #   
# #   setTxtProgressBar(progress.bar, i)
# #   time <- round(seconds_to_period(sum(end - init)), 0)
# #   
# #   
# #   # Estimated remaining time based on the
# #   # mean time that took to run the previous iterations
# #   est <- n_iter * (mean(end[end != 0] - init[init != 0])) - time
# #   remaining <- round(seconds_to_period(est), 0)
# #   
# #   cat(paste(" // Execution time:", time,
# #             " // Estimated time remaining:", remaining), "")
# # 
# # }
# # 
# # dev.off()
# # close(progress.bar)
