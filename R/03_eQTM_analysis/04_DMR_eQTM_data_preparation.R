# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"

##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("dplyr")
library("limma")
# library("biomaRt")
# library("ggplot2")
library("scales")
# library("PupillometryR")

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
input.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation"
input.dir2 <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_eQTM/01_data_prep"
output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_eQTM/04_DMR_eQTM_data_prep"
output.dir2 <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/08_DMR/"

# beta.values <- read.csv(file.path(input.dir,"08_DMR/DMR_all_CpGs_betas.csv"),
#                         header= TRUE, row.names = 1)
# colnames(beta.values) <- gsub(pattern = "X", replacement = "", colnames(beta.values))
# 
# source(file.path(input.dir, "M2beta.R"))
# 
# M.values <- beta2M(beta.values)


load("/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/06_DMA/initialData.RData")
# the initialData.RData loads in the variable "methylationData" and are M-values
colnames(methylationData) <- gsub(pattern = "X", replacement = "", colnames(methylationData))
row.names(methylationData) <- methylationData[,1] 
methylationData <- as.matrix(methylationData[,-1])

# M.values is used in script
M.values <- methylationData
rm(methylationData)



##############################################################################
# DMR.results.1 <- read.csv(file.path(input.dir, "08_DMR/DMR_identification_AsthmaVsHealthy.csv"),
#                           header = TRUE, row.names = 1)
# DMR.results.2 <- read.csv(file.path(input.dir, "08_DMR/DMR_identification_Th2HighvsHealthy.csv"),
#                           header = TRUE, row.names = 1)
# DMR.results.3 <- read.csv(file.path(input.dir, "08_DMR/DMR_identification_Th2LowvsHealthy.csv"),
#                           header = TRUE, row.names = 1)
DMR.results.3 <- read.csv(file.path(input.dir, "08_DMR/DMR_identification_ALL_Th2LowvsHealthy.csv"),
                          header = TRUE, row.names = 1)
# DMR.results.4 <- read.csv(file.path(input.dir, "08_DMR/DMR_identification_Th2HighvsTh2Low.csv"),
#                           header = TRUE, row.names = 1)
##############################################################################
CpG.annotation <- read.csv(file = file.path(input.dir, "07_anno/CpGLocation_854k.csv"),
                           header = TRUE, row.names = 1)

# Only used when quickly creating graphs to see the distribution between the groups
# eQTM.methylation.data <- read.csv(file = file.path(base.dir, "Data/005 eQTM analysis/methylationData_eQTM_samples_asthma_vs_control.csv"),
#                                   header = TRUE,
#                                   row.names = 1)
# colnames(eQTM.methylation.data) <- gsub(pattern = "X", replacement = "", colnames(eQTM.methylation.data))



# Change for preparation of different analyses
results <- DMR.results.3
# write.csv(rownames(results), file = file.path(output.dir2, "probenames_of_results.csv"))

# Manual inspection
# Some CpG sites cannot be found due to probe name not being the same
# gsub the names by removing one of the novel probes and their suffix
# Eg. "nv-GRCh38-chr17-7675160-7675160-G-A1" where the 1 after the A is not 
#   not part of the probe.
# There is a non-perl regex, but the perl version is 10% faster
#   Non-perl regex for comparison: "([A,T,C,G]-[A,T,C,G])\\d$"
# pattern <- "(?<=\\b[A,T,C,G]-[A,T,C,G])\\d"
# rownames(results) <- gsub(pattern = pattern, replacement = "", rownames(results), perl = TRUE)

DMR.mvalues <- data.frame()
DMR.annotation <- data.frame()

# ##############################################################################
# # Checking the output of the for-loop, debugging cause colmeans doesnt work
# ##############################################################################
# for (i in unique(results$DMR)) {
#   cat("\nProcessing DMR:", i, "\n")  # Debug: Track which DMR is being processed
#   
#   # Extract CpG sites within the current DMR
#   CpG.sites <- results[results$DMR == i, "CpG"]
#   
#   if (length(CpG.sites) == 0) {
#     warning("No CpG sites found for DMR:", i)
#     next  # Skip to the next iteration
#   }
#   
#   cat("CpG sites:", paste(CpG.sites, collapse = ", "), "\n")  # Debug: List CpG sites
#   
#   # Extract positional data
#   positions <- CpG.annotation[CpG.sites, c("chr", "pos", "strand"), drop = FALSE]
#   
#   if (nrow(positions) == 0) {
#     warning("No positional data found for CpG sites in DMR:", i)
#     next
#   }
#   
#   chr <- unique(positions$chr)
#   if (length(chr) > 1) {
#     warning("Multiple chromosomes detected for DMR:", i, "- Check data integrity.")
#   }
#   
#   start <- min(positions$pos, na.rm = TRUE)
#   end <- max(positions$pos, na.rm = TRUE)
#   middle <- start + round((end - start) / 2, digits = 0)
#   
#   cat("DMR range:", chr, "Start:", start, "End:", end, "Middle:", middle, "\n")  # Debug range
#   
#   positions <- data.frame(DMR = i, chr, start, middle, end)
#   DMR.annotation <- rbind(DMR.annotation, positions)
#   
#   # Extract M-values for the CpG sites
#   df <- M.values[CpG.sites, , drop = FALSE]
#   
#   if (nrow(df) == 0) {
#     warning("No M-values found for CpG sites in DMR:", i)
#     next
#   }
#   
#   cat("M-values matrix dimensions before transformation:", dim(df), "\n")  # Debug: Dimensions
#   
#   # Compute colMeans() only if there are multiple CpG sites
#   if (nrow(df) > 1) {
#     df <- as.data.frame(t(colMeans(df)))  # Ensure column means are calculated
#   } else {
#     df <- as.data.frame(t(df))  # Directly transpose single-row data
#   }
#   
#   cat("M-values matrix dimensions after transformation:", dim(df), "\n")  # Debug: New dimensions
#   # print(head(df))
#   
#   rownames(df) <- i
#   # print(head(df))
#   # df <- as.data.frame(t(df))
#   # print(head(df))
#   
#   DMR.mvalues <- rbind(DMR.mvalues, df)
#   # print(head(DMR.mvalues))
#   }
# ##############################################################################


####
# Check whether calculating the average M-value is done correctly
####
for (i in unique(results$DMR)) {
  # Withdraw CpG sites in DMR to get the M-values
  #   and printing CpG sites for troubleshooting and checking progress
  # CpG.sites <- rownames(results[results$DMR == i, ])
  CpG.sites <- results[results$DMR == i, "CpG"]

  DMR <- i
  cat(CpG.sites, "\n")

  # Calculating the important information per DMR
  # The middle position will be used for further analysis (eQTMR)
  positions <- CpG.annotation[CpG.sites, c("chr","pos", "strand")]
  chr <- unique(positions$chr)
  start <- min(positions$pos)
  end <- max(positions$pos)
  middle <- start + round((end - start) / 2, digits = 0)
  positions <- data.frame(DMR, chr, start, middle, end)
  DMR.annotation <- rbind(DMR.annotation, positions)

  # DMRs consisting of 1, will have the same value as the CpG site
  # DMRs consisting of >1 CpGs, the M-values will be averaged out
  # drop = FALSE, ensuring that DMRs consisting of single CpG site (1 entry)
  #   still keep it as a data.frame and not turn into a vector.
  df <- M.values[CpG.sites, , drop = FALSE]
  if (ncol(df) > 1) {
    df <- as.data.frame(t(colMeans(df)))
  } else {
    df <- as.data.frame(t(df))
  }
  rownames(df) <- i
  # df <- as.data.frame(t(df))
  DMR.mvalues <- rbind(DMR.mvalues, df)
}
rm(CpG.sites, DMR, positions, chr, start, end, middle, df)

# Extracting DMR data for all samples = 444 samples
write.csv(x = DMR.mvalues, file.path(output.dir2, "DMR_ALL_mval_Th2LowvsHealthy_allSamples.csv"))
# quit()

# Get the DMR name and position of eQTM analysis
DMR.annotation <- DMR.annotation[, c("DMR", "chr", "middle")]
colnames(DMR.annotation) <- c("DMR", "chr", "pos")


# Extracting DMR data for all eQTM samples = 365 samples
group_th <- read.csv(file.path(input.dir2, "covariates_eQTM_samples_IDs.csv"),
                     header = TRUE, row.names = 1)

# Asthma vs healthy 
# Don't change anything, all samples are taken
DMR.mvalues <- DMR.mvalues[, group_th$meth_file_id]

# Th2-high vs healthy
# group_th <- group_th[grep(pattern = "high|healthy", x = group_th$group_th), ]
# DMR.mvalues <- DMR.mvalues[, group_th$meth_file_id]

# Th2-low vs healthy
# group_th <- group_th[grep(pattern = "low|healthy", x = group_th$group_th), ]
# DMR.mvalues <- DMR.mvalues[, group_th$meth_file_id]

# Th2-high vs Th2-low
# group_th <- group_th[grep(pattern = "high|low", x = group_th$group_th), ]
# DMR.mvalues <- DMR.mvalues[, group_th$meth_file_id]

write.csv(x = DMR.annotation, file.path(output.dir, "DMR_ALL_pos_Th2LowvsHealthy.csv"))
write.csv(x = DMR.mvalues, file.path(output.dir, "DMR_eQTM_ALL_mval_Th2LowvsHealthy_allSamples.csv"))
          
quit()

df <- DMR.mvalues["chr15:76339951-76343297", colnames(DMR.mvalues) %in% group_th$meth_file_id]
df <- M2beta(df)
df <- eQTM.methylation.data["cg18341491",] 


df <- as.data.frame(t(df))
df$group_th <- "undeterm"
df$asthma <- group_th[group_th$meth_file_id %in% rownames(df), c("ASTHEA")]
df$group_th <- group_th[group_th$meth_file_id %in% rownames(df), c("group_th")]

df2 <- reshape2::melt(df)

COLOR_TEMP = c("#d5896f","#dab785","#70a288", "#6A9FB5")
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

max.cpg.limit <- ceiling_dec(max(df2$value), 1)
min.cpg.limit <- floor_dec(min(df2$value), 1)

p_CpG <- ggplot(data = df2) +
  aes(y = value, 
      x = variable,
      fill = group_th) + #
  geom_flat_violin(position = position_nudge(x = .20), 
                   alpha = .5) +
  geom_point(aes(color = group_th), 
             position = position_jitterdodge(jitter.width = 0.05, dodge.width = 0.3), 
               # position_jitter(width = .10),
             size = .3, 
             alpha = .5,
             show.legend = F) +
  
  geom_boxplot(width = .3, 
               outlier.shape = NA,
               alpha = .5) +
  labs(x = "CpG", y = "Beta-value") +
  scale_fill_manual(values = COLOR_TEMP) +
  scale_color_manual(values = COLOR_TEMP) +
  scale_y_continuous(limits = c(min.cpg.limit, max.cpg.limit), breaks = pretty_breaks(5)) +
  # expand_limits(x = 3, y = 10) + 
  theme_bw() 
  # theme(axis.ticks.x=element_blank(), 
  #       axis.title.x = element_blank(), 
  #       axis.title.y=element_blank(),
  #       axis.text.y=element_blank(), 
  #       axis.ticks.y=element_blank(),
  #       panel.border = element_blank(),
  #       panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank()
  # )