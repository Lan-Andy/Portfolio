# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
# output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
# output.dir.plots <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results"

##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("dplyr")
library("ggplot2")
library("limma")
library("bacon")

##############################################################################
################## Set seed for reproducible results #########################
##############################################################################
set.seed(22)

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
input.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/05_pca"
output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/06_DMA"
output.dir.plots <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/06_DMA/figures"

linkingtable <- read.csv(
  "/groups/umcg-griac/tmp02/projects/TatianaKarp/ATL_methylation/methID_rnaID_clinID.csv"
)
pheno.data <- read.csv(
  "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_gene_expression/data/atlantis_patient_data.csv",
  header = TRUE,
  row.names = 1
)
# methylationData <- read.csv(file = "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/04_normalisation/ATLANTIS_Mvalues_dasen.csv")
PC.all.cp <- read.csv(
  "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/05_pca/ctrlProbes/my_pca_data.csv"
)


# linkingtable <- read.csv(file.path(output.dir, "methID_rnaID_clinID.csv"))
# pheno.data <- read.csv(file.path(output.dir, "atlantis_patient_data.csv"),
#                        header = TRUE, row.names = 1, na.strings = "")
# PC.all.cp <- read.csv(file = file.path(output.dir, "ATLANTIS QC/Data/PCA technical probes/my_pca_data.csv"),
#                       header = TRUE)

# save(methylationData, file = file.path(output.dir, "initialData.RData"))
##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
load(file = file.path(output.dir, "initialData.RData"))
# methylationData <- as.matrix(
#   read.table(file = file.path(output.dir, "Data/sliced_methylationData_Mvalues.csv"),
#              header = TRUE,
#              sep = ",",
#              row.names = 1
#              )
# )

# If the colnames do consist an X of make.names(), replace it with 'nothing'
colnames(methylationData) <- gsub(
  pattern = "X",
  replacement = "",
  colnames(methylationData)
)
row.names(methylationData) <- methylationData[, 1]
methylationData <- as.matrix(methylationData[, -1])

# Change the VISIT column to have the same values as the big clinical table
linkingtable$VISIT <- gsub(
  pattern = "Visit 1a",
  replacement = "VISIT 1",
  x = linkingtable$VISIT
)
linkingtable$VISIT <- gsub(
  pattern = "Visit 2",
  replacement = "VISIT 2",
  x = linkingtable$VISIT
)


##############################################################################
#################### Selection methylation samples ###########################
##############################################################################
# Select the correct methylation samples
# Some samples have visit 1 and 2, some samples only have visit 2
# We chose visit 1, and only visit 2 whether they have no visit 1.
## Trying intersect ##
# sum(is.na(linkingtable$PT)) # 53 methylation samples that do not have patient ID
# 448 methylation samples (send of for sequencing) - 53 methylation samples (no patient ID)
#   396 methylation samples leftover
# table(linkingtable[!is.na(linkingtable$PT), "VISIT"])
# From 396 samples leftover:
#   329 methylation samples were from VISIT 1
#   66 methylation samples were from VISIT 2
# Filter the methylation which have a patient id and are only from visit 1
linktable1 <- linkingtable[
  !is.na(linkingtable$PT) & linkingtable$VISIT == "VISIT 1",
]

# Check the colnames of methylation samples is in the link table
methylation.patients <- intersect(
  colnames(methylationData),
  linktable1$meth_file_id
)
# linktable 2 are the visit 1 samples (even when they have visit 2, visit 1 is chosen here)
linktable2 <- linktable1[linktable1$meth_file_id %in% methylation.patients, ]

# Identifying samples with only visit 2
unique.patient.ids <- linkingtable[!is.na(linkingtable$PT), ]
single.visit2.ids <- data.frame()
for (patient.id in unique.patient.ids$PT) {
  temp <- linkingtable[linkingtable$PT %in% patient.id, c("PT", "VISIT")]
  if ("VISIT 2" %in% temp) {
    if (!"VISIT 1" %in% temp) {
      single.visit2.ids <- rbind(single.visit2.ids, temp)
      # print("There is no VISIT 1")
    } else {
      # print("There is VISIT 1")
    }
  }
  rm(temp, patient.id)
}

# linktable3 consist of:
#   1 - Samples that either have visit 1 or visit 1 and 2
#   2 - Samples that only have visit 2 methylation data
linktable3 <- rbind(
  linktable2,
  linkingtable[linkingtable$PT %in% single.visit2.ids$PT, ]
)

# write.csv(linktable3, file = file.path(output.dir, "linktable_DMA.csv"))
# quit()

pd.methyl <- pheno.data[pheno.data$PT %in% linktable2$PT, ]
pd.methyl.v1 <- pd.methyl[pd.methyl$VISIT %in% "VISIT 1", ]

pd.methyl <- pheno.data[pheno.data$PT %in% single.visit2.ids$PT, ]
pd.methyl <- pd.methyl %>% tidyr::fill(SMOKE)
pd.methyl.v2 <- pd.methyl[pd.methyl$VISIT %in% "VISIT 2", ]

pd.methyl.sub <- rbind(pd.methyl.v1, pd.methyl.v2)

# Selecting the PCs only from the samples we are interested in (see pd.methyl.sub patient.ids)
PC.all.cp <- PC.all.cp[PC.all.cp$meth_file_id %in% linktable3$meth_file_id, ]

# Merge all the dataframes to 1 big dataframe
pd.methyl.sub <- base::merge(pd.methyl.sub, linktable3, by = "PT")
pd.methyl.sub <- base::merge(pd.methyl.sub, PC.all.cp, by = "meth_file_id")

pd.methyl.sub$SMOKE <- gsub(
  pattern = "-| ",
  replacement = "_",
  x = pd.methyl.sub$SMOKE
)
pd.methyl.sub$SMOKE <- gsub(
  pattern = "Non",
  replacement = "Never",
  x = pd.methyl.sub$SMOKE
)

methylationData.sub <- methylationData[, pd.methyl.sub$meth_file_id]

################################################################################
# Determine whether the patient belongs to clinical Th2-high or Th2-low
#   Made by Tatiana Karp
pheno.data.Th2 <- pd.methyl.sub %>%
  dplyr::filter(ASTHEA == 'A') %>%
  # dplyr::mutate(include = if_else(((SYS_COR == 'Yes') | (BIO == 'Yes')), 'NO', 'YES'))%>%
  dplyr::mutate(include = dplyr::if_else((SYS_COR == 'Yes'), 'NO', 'YES')) %>%
  dplyr::mutate(
    group_th = dplyr::case_when(
      ((LABEOSV > 0.3) & (FENRES > 25)) ~ 'high',
      ((LABEOSV < 0.15) & (FENRES < 25)) ~ 'low',
      ((LABEOSV > 0.3) & (is.na(FENRES))) ~ 'high',
      ((LABEOSV < 0.15) & (is.na(FENRES))) ~ 'low'
    )
  ) %>%
  dplyr::mutate(
    group_th = dplyr::if_else(is.na(group_th), 'undeterm', group_th)
  ) %>%
  dplyr::mutate(group_th = dplyr::if_else(include == 'NO', 'cor_bio', group_th))

# Adding the Th2 group into the main pheno table
#   to be able to compare the Th2-high/low to healthy participants
pd.methyl.sub$group_th <- "healthy"
# Assign the Th2 group towards the right sample
pd.methyl.sub[
  which(pd.methyl.sub$meth_file_id %in% pheno.data.Th2$meth_file_id),
  "group_th"
] <- pheno.data.Th2$group_th

################################################################################

##### DIFFERENTIAL METHYLATION ANALYSIS (UNBIASED)#####
# To account for 95% of the total technical variance, we want to add the
#   principal components that have a cumulative proportion of 95%
# The principal components 1-12 have a cumulative proportion of 95.1%
#   Decided that wer are losing too many significant hits and the 95%
#     is an arbitrary cut-off. We have chosen to use PC1-PC3, which
#     contribute to 86% of the total technical variance.

#Design matrix
#~0 for blocking
design <- model.matrix(
  ~ 0 + ASTHEA + AGE + SEX + SMOKE + PC1 + PC2 + PC3,
  data = pd.methyl.sub
)

# duplicateCorrelation() is not needed, we are not looking over time, or we do not have replicates
#   so we do not need to block for duplicates
# #Estimate the correlation between duplicate spots (regularly spaced replicate spots on the same array) or between technical replicates from a series of arrays.
# #Blocking is done with the purpose to balance the design with respect to a factor that is known or strongly suspected to have an influence but is not in itself of interest, and it is usually assumed that block factors do not interact with experimental factors.
# corfit <- duplicateCorrelation(
#   methylationData,
#   design
# )

## Load in this object for downstream differential methylation analysis + visualisation
# load(file.path(input.dir, "Workspaces/Contmatrix.rdata"))

# fit a linear model to each cpg site
fit <- limma::lmFit(methylationData.sub, design)


# Specify columns to compare
cont.matrix <- limma::makeContrasts(
  test = ASTHEAA - ASTHEAH,
  levels = design
)

# contrasts.fit converts the coefficients and standard errors to reflect the contrasts rather than the original design matrix, but does not compute t-statistics or p-values.
fit2 <- limma::contrasts.fit(fit, contrast = cont.matrix[, 'test'])


# eBayes is used to rank sites in order of evidence for differential methylation
# eBayes computes t-statistics and p-values from the coefficients and standard errors.
fit2 <- limma::eBayes(fit2)


# For future DMR analysis, calculating the standard error per CpG site
fit2.SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
write.csv(
  fit2.SE,
  file = file.path(output.dir, "fit2_standard_error_DMA_asthma_vs_control.csv")
)

# topTableF ranks cpg sites on the basis of moderated F-statistics for one or more coefficients. FDR calculated using Benjamin-Hochberg method
tT.BH <- limma::topTable(
  fit2,
  adjust = "BH",
  sort.by = "P",
  number = nrow(fit2)
)
tT.Bonf <- limma::topTable(
  fit2,
  adjust = "bonferroni",
  sort.by = "P",
  number = nrow(fit2)
)

tT <- tT.BH
tT$adj.P.BH <- tT.BH$adj.P.Val
tT$adj.P.Bonf <- tT.Bonf$adj.P.Val

################################################################################
# Using R package bacon to calculate the bias and genomic inflation factor
# Apply bacon correction to t-statistics
bacon_results <- bacon(tT$t)

# Extract bias and inflation factor
bias <- bacon_results@estimates[1, "mu.0"]
inflation <- bacon_results@estimates[1, "sigma.0"]

# Print the estimated bias and genomic inflation factor
cat(
  "For the comparison Asthma vs. Healthy:\n",
  "Estimated bias is:",
  bias,
  "and the genomic inflation factor is:",
  inflation,
  "for this analysis.\n"
)

# Comparison        Bias  Genomic Inflation Factor
# asthma vs control 0.01  1.08

# Correct t-statistics using bacon estimates
bacon_corrected_t <- (bacon_results@teststatistics - bias) / inflation
str(bacon_corrected_t)
class(bacon_corrected_t)


################################################################################
#### Calculate the new bias and inflation factor after correction
bacon_results <- bacon(as.vector(bacon_corrected_t))
bias <- bacon_results@estimates[1, "mu.0"]
inflation <- bacon_results@estimates[1, "sigma.0"]

# Print the estimated bias and genomic inflation factor
cat(
  "For the comparison Asthma vs. Healthy:\n",
  "New estimated bias is:",
  bias,
  "and the new genomic inflation factor is:",
  inflation,
  "for this analysis.\n"
)
####
################################################################################
# Add new columns with consistent and concise names
# Bacon-corrected t-statistics and two-sided p-value
tT$t.Bacon <- as.vector(as.vector(bacon_corrected_t))
tT$P.Value.Bacon <- 2 * pnorm(-abs(tT$t.Bacon))
# Multiple testing over bacon corrected p-values
tT$adj.P.Val.Bacon.BH <- p.adjust(tT$P.Value.Bacon, method = "BH")
tT$adj.P.Val.Bacon.Bonf <- p.adjust(tT$P.Value.Bacon, method = "bonferroni")

# Reorder columns for clarity
column_order <- c(
  "logFC",
  "AveExpr",
  "t",
  "P.Value",
  "adj.P.BH",
  "adj.P.Bonf",
  "t.Bacon",
  "P.Value.Bacon",
  "adj.P.Val.Bacon.BH",
  "adj.P.Val.Bacon.Bonf",
  "B"
)

# Ensure only existing columns are selected (prevents errors)
column_order <- base::intersect(column_order, colnames(tT))
tT <- tT[, column_order]

# Add logical columns indicating significance (TRUE/FALSE) for each correction
tT$Significant.Orig.BH <- tT$adj.P.BH < 0.05
tT$Significant.Orig.Bonf <- tT$adj.P.Bonf < 0.05
tT$Significant.Bacon.BH <- tT$adj.P.Val.Bacon.BH < 0.05
tT$Significant.Bacon.Bonf <- tT$adj.P.Val.Bacon.Bonf < 0.05

# Save the results
write.csv(
  tT,
  file.path(output.dir, "DMA_tT_CpG_asthma_vs_control_all_corrections.csv")
)

#### Changed the code below to the one here above, due to testing of the genomic inflation factor
# # topTableF ranks cpg sites on the basis of moderated F-statistics for one or more coefficients. FDR calculated using Benjamin-Hochberg method
# tT <- limma::topTable(fit2, adjust="BH", sort.by="P", number=nrow(fit2))
#
# # selecting sites that have adj.p.val less than 0.05
# selection <- which(tT$adj.P.Val<0.05)
#
# # Only includes sites that have significant FDR (Adj p val < 0.05)
# tT2 <- tT[selection,]
#
# # Save the results
# write.csv(tT2, file.path(output.dir, "Tt2_significant_CpG_asthma_vs_control_BH.csv"))
# write.csv(tT, file.path(output.dir, "Tt_CpG_asthma_vs_control_BH.csv"))
#
# # Redo the multiple testing correction with bonferroni
# tT <- limma::topTable(fit2, adjust="bonferroni", sort.by="P", number=nrow(fit2))
# selection <- which(tT$adj.P.Val<0.05)
#
# tT2 <- tT[selection,]
# write.csv(tT2, file.path(output.dir, "Tt2_significant_CpG_asthma_vs_control_bonferroni.csv"))
# write.csv(tT, file.path(output.dir, "Tt_CpG_asthma_vs_control_bonferroni.csv"))
################################################################################

################################################################################
# Performing multiple comparisons all in one go
# Th2 high vs Healthy
# Th2 low vs Healthy
# Th2 high vs Th2 low

design <- model.matrix(
  ~ 0 + group_th + AGE + SEX + SMOKE + PC1 + PC2 + PC3,
  data = pd.methyl.sub
)
# Remove the "group_th" prefix, for easier interpreting names
colnames(design) <- gsub(
  pattern = "group_th",
  replacement = "",
  x = colnames(design)
)

fit <- limma::lmFit(methylationData.sub, design)

cont.matrix <- limma::makeContrasts(
  Th2HighVsHealthy = high - healthy,
  Th2LowVsHealthy = low - healthy,
  Th2HighvsTh2Low = high - low,
  levels = colnames(design)
)

# Perform multiple comparison depending on the amount of contrasts that have been made
for (hypothesis in colnames(cont.matrix)) {
  print(hypothesis)
  fit2 <- limma::contrasts.fit(fit, contrast = cont.matrix[, hypothesis])
  fit2 <- limma::eBayes(fit2)
  # # For future DMR analysis, calculating the standard error per CpG site
  fit2.SE <- sqrt(fit2$s2.post) * fit2$stdev.unscaled
  write.csv(
    fit2.SE,
    file = file.path(
      output.dir,
      paste0("fit2_standard_error_DMA", hypothesis, ".csv")
    )
  )

  # Similarly to the code above whilst comparing asthma vs healthy
  tT.BH <- limma::topTable(
    fit2,
    adjust = "BH",
    sort.by = "P",
    number = nrow(fit2)
  )
  tT.Bonf <- limma::topTable(
    fit2,
    adjust = "bonferroni",
    sort.by = "P",
    number = nrow(fit2)
  )

  tT <- tT.BH
  tT$adj.P.BH <- tT.BH$adj.P.Val
  tT$adj.P.Bonf <- tT.Bonf$adj.P.Val

  ##############################################################################
  #### Genomic inflation factors were calculated for all comparisons
  #### Due to the low lambda from all comparisons, it is not needed to correct
  ####  for lambda. Therefore we skip bacon corrected t-statistic
  bacon_results <- bacon(tT$t)

  # Extract bias and inflation factor
  bias <- bacon_results@estimates[1, "mu.0"]
  inflation <- bacon_results@estimates[1, "sigma.0"]

  # Print the estimated bias and genomic inflation factor
  cat(
    "For the comparison ",
    hypothesis,
    ":\n",
    "Estimated bias is:",
    bias,
    "and the genomic inflation factor is:",
    inflation,
    "for this analysis.\n"
  )

  # Comparison        Bias  Genomic Inflation Factor
  # Th2HighVsHealthy  -0.08  1.31
  # Th2LowVsHealthy    0.03  1.03
  # Th2HighvsTh2Low   -0.07  1.07

  # Correct t-statistics using bacon estimates
  bacon_corrected_t <- (bacon_results@teststatistics - bias) / inflation
  ##############################################################################
  #### Calculate the new bias and inflation factor after correction ####
  bacon_results <- bacon(as.vector(bacon_corrected_t))
  bias <- bacon_results@estimates[1, "mu.0"]
  inflation <- bacon_results@estimates[1, "sigma.0"]

  # Print the estimated bias and genomic inflation factor
  cat(
    "For the comparison ",
    hypothesis,
    ":\n",
    "New estimated bias is:",
    bias,
    "and the new genomic inflation factor is:",
    inflation,
    "for this analysis.\n"
  )
  ####                                                              ####
  ##############################################################################
  # Add new columns with consistent and concise names
  # Bacon-corrected t-statistics and two-sided p-value
  tT$t.Bacon <- as.vector(bacon_corrected_t)
  tT$P.Value.Bacon <- 2 * pnorm(-abs(tT$t.Bacon))
  # Multiple testing over bacon corrected p-values
  tT$adj.P.Val.Bacon.BH <- p.adjust(tT$P.Value.Bacon, method = "BH")
  tT$adj.P.Val.Bacon.Bonf <- p.adjust(tT$P.Value.Bacon, method = "bonferroni")

  # Reorder columns for clarity
  column_order <- c(
    "logFC",
    "AveExpr",
    "t",
    "P.Value",
    "adj.P.BH",
    "adj.P.Bonf",
    "t.Bacon",
    "P.Value.Bacon",
    "adj.P.Val.Bacon.BH",
    "adj.P.Val.Bacon.Bonf",
    "B"
  )

  # Ensure only existing columns are selected (prevents errors)
  column_order <- base::intersect(column_order, colnames(tT))
  tT <- tT[, column_order]

  # Add logical columns indicating significance (TRUE/FALSE) for each correction
  tT$Significant.Orig.BH <- tT$adj.P.BH < 0.05
  tT$Significant.Orig.Bonf <- tT$adj.P.Bonf < 0.05
  tT$Significant.Bacon.BH <- tT$adj.P.Val.Bacon.BH < 0.05
  tT$Significant.Bacon.Bonf <- tT$adj.P.Val.Bacon.Bonf < 0.05

  # Save the results
  write.csv(
    tT,
    file.path(
      output.dir,
      paste0("DMA_tT_CpG_", hypothesis, "_all_corrections.csv")
    )
  )
  ##############################################################################

  ##############################################################################
  # #### Change the code below to the one here above, due to testing of the genomic inflation factor
  # #### Keeping the old code for archive sake
  # # Benjamini-Hochberg
  # tT <- limma::topTable(fit2, adjust="BH", sort.by = "P", number=nrow(fit2))
  # write.csv(tT, file.path(output.dir, paste0("DMA_Tt_CpG_",hypothesis,"_BH.csv")))
  #
  # selection <- which(tT$adj.P.Val<0.05)
  # tT2 <- tT[selection,]
  # write.csv(tT2, file.path(output.dir, paste0("DMA_Tt2_significant_CpG_",hypothesis,"_BH.csv")))
  #
  # # Bonferroni
  # tT <- limma::topTable(fit2, adjust="bonferroni", sort.by = "P", number=nrow(fit2))
  # write.csv(tT, file.path(output.dir, paste0("DMA_Tt_CpG_",hypothesis,"_bonferroni.csv")))
  #
  # selection <- which(tT$adj.P.Val<0.05)
  # tT2 <- tT[selection,]
  # write.csv(tT2, file.path(output.dir, paste0("DMA_Tt2_significant_CpG_",hypothesis,"_bonferroni.csv")))
  # ####
}


quit()
