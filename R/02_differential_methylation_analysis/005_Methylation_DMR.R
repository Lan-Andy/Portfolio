# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
# input.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
base.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/"
output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/08_DMR"

##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("dmrff")

# ##############################################################################
# #################### Load in raw data for analysis ###########################
# ##############################################################################
# Changing this and the standard error file, depending on the hypothesis/input files
# Asthma vs Healthy   = 06_DMA/Tt_CpG_asthma_vs_control_bonferroni.csv
#   Standard error    = 06_DMA/fit2_standard_error_DMA_asthma_vs_control.csv

# Th2-high vs Healthy = 06_DMA/DMA_Tt_CpG_Th2HighVsHealthy_bonferroni.csv
#   Standard error    = 06_DMA/fit2_standard_error_DMATh2HighVsHealthy.csv

# Th2-low vs Healthy  = 06_DMA/DMA_Tt_CpG_Th2LowVsHealthy_bonferroni.csv
#   Standard error    = 06_DMA/fit2_standard_error_DMATh2LowVsHealthy.csv

# Th2-high vs Th2-low = 06_DMA/DMA_Tt_CpG_Th2HighvsTh2Low_bonferroni.csv
#   Standard error    = 06_DMA/fit2_standard_error_DMATh2HighvsTh2Low.csv
DMA.results <- read.csv(file = file.path(base.dir, "06_DMA/Tt_CpG_asthma_vs_control_bonferroni.csv"),
                    header = TRUE,
                    row.names = 1)

DMA_standard_error <- read.csv(file = file.path(base.dir, "06_DMA/fit2_standard_error_DMA_asthma_vs_control.csv"),
                               header = TRUE,
                               row.names = 1)


CpG_annotation <- read.csv(file = file.path(base.dir, "07_anno/CpGLocation_854k.csv"),
                           header = TRUE,
                           row.names = 1)
# methylation <-  read.csv(file = file.path(base.dir, "04_normalisation/ATLANTIS_betas_dasen.csv"),
#                          header = TRUE,
#                          row.names = 1)
# save(methylation, file = file.path(output.dir, "initialData_ATLANTIS_betas.RData"))
load(file = file.path(output.dir, "initialData_ATLANTIS_betas.RData"))
# 
linktable <- read.csv(file = file.path(base.dir, "06_DMA/ATLANTIS_linkingtable_with_group_Th.csv"),
                      header = TRUE,
                      row.names = 1)

DMA.subjects <- read.csv(file = file.path(base.dir, "06_DMA/linktable_DMA.csv"))

print("Finished loading all the data for the DMR analysis")

##############################################################################
######################### Select and filter raw data #########################
##############################################################################
# linktable <- read.csv(file = file.path(input.dir, "Data/100 Characteristics tables/ATLANTIS_linkingtable_with_group_Th.csv"),
#                       header = TRUE,
#                       row.names = 1)
# Removing the samples where asthma status (= unknown clinical data)


linktable <- linktable[which(linktable$meth_file_id %in% DMA.subjects$meth_file_id), ]
linktable <- linktable[!is.na(linktable$ASTHEA), ]
dim(linktable)

asthma_samples <- linktable[linktable$ASTHEA == "A", "meth_file_id"]
healthy_samples <- linktable[linktable$ASTHEA == "H",  "meth_file_id"]
Th2.high_samples <- linktable[linktable$group_th == "high", "meth_file_id"]
Th2.low_samples <- linktable[linktable$group_th == "low",  "meth_file_id"]

# Change to samples of interest
select_samples <- c(asthma_samples, healthy_samples)
# select_samples <- c(Th2.high_samples, healthy_samples)
# select_samples <- c(Th2.low_samples, healthy_samples)
# select_samples <- c(Th2.high_samples, Th2.low_samples)

colnames(methylation) <- gsub(pattern = "X", replacement = "", x = colnames(methylation))
methylation <- methylation[, colnames(methylation) %in% select_samples]

print("Starting to convert the beta-values into M-values")
source(file.path(base.dir,"M2beta.R"))
methylation <- as.matrix(beta2M(methylation))
print("Done converting beta-values")


print(head(methylation))
print(str(methylation))
print(dim(methylation))

DMA.results$SE <- DMA_standard_error[rownames(DMA.results), ]
DMA.results$chr <- CpG_annotation[rownames(DMA.results), "chr"]
DMA.results$pos <- CpG_annotation[rownames(DMA.results), "pos"]

print(head(DMA.results))
print(dim(DMA.results))
##############################################################################
#################### Differential methylated region analysis #################
##############################################################################
# Perform the DMR analysis

# Interesting to look into the different gap lengths whether there is an optimal significance threshold
# Investigating what will happen with the amount of significant DMRs 
#   when we increase the maximum gap length between each CpG site
# 27/07/2024: Found that the amount of DMRs increases due to the amount of tests that are performed
#   Investigated the ratio of amount of DMRs to the amount of tests, flattens out at 750 gap 
# gap.lengths <- c(100, seq(from = 250, to = 2500, by = 250))
# 
# for (i in gap.lengths) {
#   dmrs <- dmrff(
#     estimate = DMA.results$logFC,
#     se = DMA.results$SE,
#     p.value = DMA.results$P.Value,
#     chr = DMA.results$chr,
#     pos = DMA.results$pos,
#     methylation = as.matrix(methylation),
#     maxgap = i
#   )
#   write.csv(x = dmrs, file = file.path(output.dir, paste0("DMR_results_", i, "gap.csv")))
#   print(paste0("DMR analysis with gap of ", i, " has been performed and saved."))
# }
# 
# quit()
print("Starting the DMR analysis")
dmrs <- dmrff(
  estimate = DMA.results$logFC,
  se = DMA.results$SE,
  p.value = DMA.results$P.Value,
  chr = DMA.results$chr,
  pos = DMA.results$pos,
  methylation = methylation,
  maxgap = 750
)

# save(dmrs, file = file.path(output.dir, "DMR_results_750gap_asthma_vs_control.RData"))
# load(file = file.path(output.dir, "DMR_results_500gap.RData"))

print("Saving the output of the DMR analysis")
write.csv(x = dmrs, file = file.path(output.dir, "DMR_results_asthma_vs_control_bonferroni_750gap.csv"))




