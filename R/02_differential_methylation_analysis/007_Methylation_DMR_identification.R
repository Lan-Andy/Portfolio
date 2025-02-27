# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
base.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation"


##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("GenomicRanges")
library("reshape2")
library("ggplot2")
library("PupillometryR")
library("scales")
library("plotly")

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
# cpg.location <- read.csv(file.path(base.dir, "Data/005 eQTM analysis/CpGLocation_860k.csv"),
#                          header = TRUE, row.names = 1)
# beta.values <- read.csv(file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/sig_ATLANTIS_betas_dasen.csv"),
#                         header = TRUE, row.names = 1)
# covariates <- read.csv(file = file.path(base.dir, "Data/005 eQTM analysis/covariates_eQTM_samples_original.csv"),
#                        header = TRUE, row.names = 1)
# 
# 
# DMR.result <- read.csv(file.path(base.dir, "Results/006 Differential Methylated Region Analysis DMR/DMR_results_750gap.csv"),
#                        header = TRUE, row.names = 1)
# DMA.result <- read.csv(file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3/Tt_CpG_asthma_vs_control_BH.csv"),
#                        header = TRUE, row.names = 1)
# 
# DMR.result.highvshealthy <- read.csv(file.path(base.dir, "Results/006 Differential Methylated Region Analysis DMR/DMR_results_Th2HighvsHealthy_750gap.csv"),
#                        header = TRUE, row.names = 1)
# DMA.result.highvshealthy <- read.csv(file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3 Th2 groups/DMA_Tt_CpG_Th2HighVsHealthy_BH.csv"),
#                                      header = TRUE, row.names = 1)
# 
# DMR.result.lowvshealthy <- read.csv(file.path(base.dir, "Results/006 Differential Methylated Region Analysis DMR/DMR_results_Th2LowvsHealthy_750gap.csv"),
#                        header = TRUE, row.names = 1)
# DMA.result.lowvshealthy <- read.csv(file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3 Th2 groups/DMA_Tt_CpG_Th2LowVsHealthy_BH.csv"),
#                                      header = TRUE, row.names = 1)
# 
# DMR.result.highvslow <- read.csv(file.path(base.dir, "Results/006 Differential Methylated Region Analysis DMR/DMR_results_Th2HighvsTh2Low_750gap.csv"),
#                        header = TRUE, row.names = 1)
# DMA.result.highvslow <- read.csv(file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3 Th2 groups/DMA_Tt_CpG_Th2HighvsTh2Low_BH.csv"),
#                                      header = TRUE, row.names = 1)

################################################################################
####                         !!!!      IMPORTANT    !!!!                    ####
####    Run these lines AFTER creating all the identification files !!!     ####
################################################################################
# base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/006 Differential Methylated Region Analysis DMR"
# # Create a subset of the raw beta values to create locally the boxplots for specific CpG sites
# DMR.ident.1 <- read.csv(file = file.path(paste0(base.dir, "/DMR_identification_", "AsthmaVsHealthy", ".csv")),
#                         header = TRUE, row.names = 1)
# DMR.ident.2 <- read.csv(file = file.path(paste0(base.dir, "/DMR_identification_", "Th2HighvsHealthy", ".csv")),
#                         header = TRUE, row.names = 1)
# DMR.ident.3 <- read.csv(file = file.path(paste0(base.dir, "/DMR_identification_", "Th2LowvsHealthy", ".csv")),
#                         header = TRUE, row.names = 1)
# DMR.ident.4 <- read.csv(file = file.path(paste0(base.dir, "/DMR_identification_", "Th2HighvsTh2Low", ".csv")),
#                         header = TRUE, row.names = 1)

# # Create a subset of the raw beta values to create locally the boxplots for specific CpG sites
# DMR.ident.1 <- read.csv(file = file.path(paste0(base.dir, "/08_DMR/DMR_identification_", "AsthmaVsHealthy", ".csv")),
#                         header = TRUE, row.names = 1)
# DMR.ident.2 <- read.csv(file = file.path(paste0(base.dir, "/08_DMR/DMR_identification_", "Th2HighvsHealthy", ".csv")),
#                         header = TRUE, row.names = 1)
# DMR.ident.3 <- read.csv(file = file.path(paste0(base.dir, "/08_DMR/DMR_identification_", "Th2LowvsHealthy", ".csv")),
#                         header = TRUE, row.names = 1)
# DMR.ident.4 <- read.csv(file = file.path(paste0(base.dir, "/08_DMR/DMR_identification_", "Th2HighvsTh2Low", ".csv")),
#                         header = TRUE, row.names = 1)
# 
# all.unique.CpGs <- unique(c(rownames(DMR.ident.1), 
#                             rownames(DMR.ident.2),
#                             rownames(DMR.ident.3), 
#                             rownames(DMR.ident.4))
# )
# loadRData <- function(fileName){
#   #loads an RData file, and returns it
#   load(fileName)
#   get(ls()[ls() != "fileName"])
# }
# 
# methylation <- loadRData(fileName = file.path(base.dir, "08_DMR/initialData_ATLANTIS_betas.RData"))
# 
# print(methylation[1:5, 1:5])
# 
# write.csv(x = methylation[all.unique.CpGs, ], file = file.path(base.dir, "08_DMR/DMR_all_CpGs_betas.csv"))
# 
# 
# quit()

################################################################################
cpg.location <- read.csv(file.path(base.dir, "07_anno/CpGLocation_854k.csv"),
                         header = TRUE, row.names = 1)
# beta.values <- read.csv(file.path(base.dir, "04_normalisation/ATLANTIS_betas_dasen.csv"),
#                         header = TRUE, row.names = 1)
# load(file = file.path(base.dir, "08_DMR/initialData_ATLANTIS_betas.RData"))
covariates <- read.csv(file = "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_eQTM/01_data_prep/covariates_eQTM_samples_original.csv",
                       header = TRUE, row.names = 1)

DMR.result <- read.csv(file.path(base.dir, "08_DMR/DMR_results_asthma_vs_control_750gap.csv"),
                       header = TRUE, row.names = 1)
DMA.result <- read.csv(file.path(base.dir, "06_DMA/Tt_CpG_asthma_vs_control_BH.csv"),
                       header = TRUE, row.names = 1)

DMR.result.highvshealthy <- read.csv(file.path(base.dir, "08_DMR/DMR_results_Th2HighVsHealthy_750gap.csv"),
                                     header = TRUE, row.names = 1)
DMA.result.highvshealthy <- read.csv(file.path(base.dir, "06_DMA/DMA_Tt_CpG_Th2HighVsHealthy_BH.csv"),
                                     header = TRUE, row.names = 1)

DMR.result.lowvshealthy <- read.csv(file.path(base.dir, "08_DMR/DMR_results_Th2LowVsHealthy_750gap.csv"),
                                    header = TRUE, row.names = 1)
DMA.result.lowvshealthy <- read.csv(file.path(base.dir, "06_DMA/DMA_Tt_CpG_Th2LowVsHealthy_BH.csv"),
                                    header = TRUE, row.names = 1)

DMR.result.highvslow <- read.csv(file.path(base.dir, "08_DMR/DMR_results_Th2HighvsTh2Low_750gap.csv"),
                                 header = TRUE, row.names = 1)
DMA.result.highvslow <- read.csv(file.path(base.dir, "06_DMA/DMA_Tt_CpG_Th2HighvsTh2Low_BH.csv"),
                                 header = TRUE, row.names = 1)

##############################################################################
#################### Transform data for analysis and visualisation ###########
##############################################################################
# beta.values <- methylation
# rm(methylation)
# colnames(beta.values) <- gsub(pattern = "X", replacement = "", x = colnames(beta.values))
# print(head(beta.values))

# # Filtering out the regions with a single CpG site
# DMR.result <- DMR.result[-which(DMR.result$n == 1), ]
# DMR.result.highvshealthy <- DMR.result.highvshealthy[-which(DMR.result.highvshealthy$n == 1), ]
# DMR.result.lowvshealthy <- DMR.result.lowvshealthy[-which(DMR.result.lowvshealthy$n == 1), ]
# DMR.result.highvslow <- DMR.result.highvslow[-which(DMR.result.highvslow$n == 1), ]

# head(DMR.result)
# reshape2::melt(DMR.result)


##############################################################################
#################### Identifying CpG sites within the DMRs ###################
##############################################################################
cpg.ranges <- GRanges(cpg.location$chr, IRanges(cpg.location$pos, cpg.location$pos))

# Create list of all the variables to iterate through them
listofDMRs <- list("AsthmaVsHealthy" = DMR.result,
                   "Th2HighvsHealthy" = DMR.result.highvshealthy,
                   "Th2LowvsHealthy" = DMR.result.lowvshealthy,
                   "Th2HighvsTh2Low" = DMR.result.highvslow
                   )

listofCpGs <- list("AsthmaVsHealthy" = DMA.result,
                   "Th2HighvsHealthy" = DMA.result.highvshealthy,
                   "Th2LowvsHealthy" = DMA.result.lowvshealthy,
                   "Th2HighvsTh2Low" = DMA.result.highvslow
                   )

comparisons <- c("Th2LowvsHealthy")

# c("AsthmaVsHealthy",
#   "Th2HighvsHealthy",
#   "Th2LowvsHealthy",
#   "Th2HighvsTh2Low"
# )

# save.image(file=file.path(base.dir, "08_DMR/DMR_identification_workshop.RData"))
# load(file.path(base.dir, "08_DMR/DMR_identification_workshop.RData"))

for (comparison in comparisons) {
  
  # To know which comparison we are in
  print(comparison)
  
  # Select the correct DMR results from the correct comparison
  DMRs <- listofDMRs[[comparison]]
  # Re-order the results based on the bonferroni adjusted p-value
  DMRs <- DMRs[order(DMRs$p.adjust, decreasing = FALSE), ]
  
  # Select the correct CpG topTable from the correct comparison
  CpGs <- listofCpGs[[comparison]]
  
  # Filter only the DMRs that are FDR (bonferroni) significant
  filteredDMRs <- DMRs
  # filteredDMRs <- DMRs[DMRs$p.adjust < 0.05, ]
  
  # Create dummy dataframe to append results to
  df <- data.frame()
  
  
  if (class(df) == "data.frame") {
    for (i in 1:nrow(filteredDMRs)) {
      print(i)
      # Create GRange object of the DMR, so we can find the CpG sites within that region
      DMR <- GRanges(filteredDMRs[i, "chr"], IRanges(filteredDMRs[i, "start"], filteredDMRs[i, "end"]))
      # Find the CpG sites within the DMR region
      findCpGs <- findOverlaps(cpg.ranges, DMR)
      # Extract the chromosome, position and strand of the CpG sites
      foundCpGs <- as.data.frame(cpg.location[findCpGs@from, c("chr", "pos", "strand")])
      # 2025 01 30 - When extracting all DMRs, some DMRs have overlapping CpGs (?)
      #   It was noticed that dmrff takes a window, but excludes the position itself
      #   Changes the position of the CpG only to 1bp instead of 2bp
      foundCpGs <- data.frame(
        "CpG" = rownames(foundCpGs), "chr" = foundCpGs$chr, "pos" = foundCpGs$pos, "strand" = foundCpGs$strand
      )
      # Extract the logFoldChange of the CpG sites
      foundCpGs$CpG_logFC <- CpGs[foundCpGs$CpG, "logFC"]
      # Extract the BH correct p-value of the CpG sites
      foundCpGs$CpG_pAdjust <- CpGs[foundCpGs$CpG, "adj.P.Val"]
      # Rename the DMR to following format "chr:start-end"
      foundCpGs$DMR <- paste0(filteredDMRs[i, "chr"], ":", filteredDMRs[i, "start"], "-", filteredDMRs[i, "end"])
      # Extract the leftover results of the DMRs from the DMR analysis
      foundCpGs$DMR_estimate <- filteredDMRs[i, "estimate"]
      foundCpGs$DMR_se <- filteredDMRs[i, "se"]
      foundCpGs$DMR_pAdjust <- filteredDMRs[i, "p.adjust"]
      # print(dim(foundCpGs))
      df <- rbind(df, foundCpGs)
      
    }
    print("DONE")
    write.csv(x = df, file = file.path(paste0(base.dir, "/08_DMR/DMR_identification_ALL_", comparison, ".csv")))
    # write.csv(x = df, file = file.path(paste0(base.dir, "/08_DMR/DMR_identification_", comparison, ".csv")))
    
    # # Filter the top 5 significant DMRs and create boxplots of the CpG sites within these DMRs
    # df2 <- df[which(df$DMR %in% unique(df$DMR)[1:5]), ]
    # 
    # n <- 0
    # 
    # for (j in unique(df2$DMR)) {
    #   n <- n + 1
    #   print(n)
    #   visual.data <- beta.values[rownames(df[df$DMR == j, ]), covariates$meth_file_id]
    #   visual.data <- as.data.frame(t(visual.data))
    #   visual.data$asthma <- covariates[which(covariates$meth_file_id %in% rownames(visual.data)), "ASTHEA"]
    #   visual.data$asthma <- gsub(pattern = "1", replacement = "Asthma", visual.data$asthma)
    #   visual.data$asthma <- gsub(pattern = "0", replacement = "Healthy", visual.data$asthma)
    # 
    #   head(visual.data)
    #   print(str(visual.data))
    #   print(colnames(visual.data))
    #   # visual.data$asthma <- as.factor(visual.data$asthma)
    #   # visual.data <- visual.data[, -grep(pattern = "NA*", x = colnames(visual.data))]
    # 
    #   # print("Printing the data.frame 'Visual.data'")
    #   # print(visual.data)
    #   # print("About to melt the data.frame 'visual.data'")
    #   # visual.data.m <- reshape2::melt(visual.data)
    #   # print("Printing the data.frame 'visual.data.m'")
    #   # print(visual.data.m)
    # 
    # 
    #   # print("About to melt the data.frame 'visual.data' with id.vars")
    #   # Transform the datafrane to create 1 figure with multiple CpG sites within
    #   visual.data.m <- reshape2::melt(visual.data, id.vars = "asthma")
    # 
    #   # print("Printing the data.frame 'visual.data.m'  with id.vars")
    #   # print(visual.data.m)
    # 
    #   boxplots <- ggplot(data = visual.data.m,
    #          aes(x = variable,
    #              y = value,
    #              fill = asthma)) +
    #     geom_flat_violin(position = position_nudge(x = .2),
    #                      alpha = .4) +
    #     geom_point(aes(colour = asthma),
    #                position = position_jitterdodge(jitter.width = 0.075, dodge.width = 0.25),
    #                size = 1,
    #                alpha = 0.4,
    #                show.legend = F) +
    #     geom_boxplot(width = .25,
    #                  outlier.shape = NA,
    #                  alpha = 0.5) +
    #     scale_fill_manual(values = c("#ff7f0e","#1f77b4")) +
    #     scale_color_manual(values = c("#ff7f0e","#1f77b4"))  +
    #     scale_y_continuous(#limits = c(0,1),
    #       breaks = pretty_breaks(n = 5)
    #     ) +
    #     labs(x = "",
    #          y = "",
    #          fill = "Status ",
    #          title = "") +
    #     guides(fill = guide_legend(nrow=1,
    #                                byrow=TRUE))+
    #     theme_bw() +
    #     theme(legend.position = "bottom",
    #           legend.margin = margin(-5, 0, 0, 0))
    # 
    #   # print(boxplots)
    #   # ggsave(plot = boxplots,
    #   #        path = file.path(base.dir, "Results/006 Differential Methylated Region Analysis DMR"),
    #   #        filename = paste0(comparison,"_DMR", n, ".tiff"))
    #   ggsave(plot = boxplots,
    #          path = file.path(base.dir, "08_DMR"),
    #          filename = paste0(comparison,"_DMR", n, ".tiff"))
    # 
    #   print(paste0("Saved plot DMR: ", j, " of [", n, "] of ", comparison))
    # }
    
    
  }
}


quit()


##############################################################################
#################### Visualisation of Th2-high/low vs healthy ################
##############################################################################
# Ran on computer NOT cluster

DMA.Th2HighHealthy <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/002 Differential Methylation Analysis DMA/Correct for PC1-3 Th2 groups/DMA_Tt_CpG_Th2HighVsHealthy_BH.csv",
                               header = TRUE,
                               row.names = 1)
DMR.Th2HighHealthy <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/006 Differential Methylated Region Analysis DMR/DMR_identification_Th2HighvsHealthy.csv",
                               header = TRUE,
                               row.names = 1)

DMA.Th2LowHealthy <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/002 Differential Methylation Analysis DMA/Correct for PC1-3 Th2 groups/DMA_Tt_CpG_Th2LowVsHealthy_BH.csv",
                              header = TRUE,
                              row.names = 1)
DMR.Th2lowHealthy <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/006 Differential Methylated Region Analysis DMR/DMR_identification_Th2LowvsHealthy.csv",
                              header = TRUE,
                              row.names = 1)

topDMR <- 100
DMR.Th2HighHealthy2 <- DMR.Th2HighHealthy[DMR.Th2HighHealthy$DMR %in% unique(DMR.Th2HighHealthy$DMR)[1:topDMR], ]
DMR.Th2lowHealthy2 <- DMR.Th2lowHealthy[DMR.Th2lowHealthy$DMR %in% unique(DMR.Th2lowHealthy$DMR)[1:topDMR], ]



Th2high.healthy.cpgs <- rownames(DMR.Th2HighHealthy2)
df <- data.frame()
df <- data.frame(DMA.Th2HighHealthy[Th2high.healthy.cpgs, "logFC"])
rownames(df) <- Th2high.healthy.cpgs
colnames(df) <- "HighHealthy_LogFC"
df$HighHealthy_adjPValue <- signif(DMA.Th2HighHealthy[Th2high.healthy.cpgs, "adj.P.Val"], 3)
df$LowHealthy_LogFC <- DMA.Th2LowHealthy[Th2high.healthy.cpgs, "logFC"]
df$LowHealthy_adjPValue <- signif(DMA.Th2LowHealthy[Th2high.healthy.cpgs, "adj.P.Val"], 3)
df$comparison <- "Th2High"

Th2low.healthy.cpgs <- rownames(DMR.Th2lowHealthy2)
df2 <- data.frame()
df2 <- data.frame(DMA.Th2HighHealthy[Th2low.healthy.cpgs, "logFC"])
rownames(df2) <- Th2low.healthy.cpgs
colnames(df2) <- "HighHealthy_LogFC"
df2$HighHealthy_adjPValue <- signif(DMA.Th2HighHealthy[Th2low.healthy.cpgs, "adj.P.Val"], 3)
df2$LowHealthy_LogFC <- DMA.Th2LowHealthy[Th2low.healthy.cpgs, "logFC"]
df2$LowHealthy_adjPValue <- signif(DMA.Th2LowHealthy[Th2low.healthy.cpgs, "adj.P.Val"], 3)
df2$comparison <- "Th2Low"

overlapping.CpGs <- intersect(rownames(DMR.Th2HighHealthy2), rownames(DMR.Th2lowHealthy2))

df[rownames(df) %in% overlapping.CpGs, "comparison"] <- "Both"
df2[rownames(df2) %in% overlapping.CpGs, "comparison"] <- "Both"

df3 <- rbind(df, df2)
df3$CpG <- rownames(df3)

scatterlogFC <- ggplot(data = df3,
       aes(x = LowHealthy_LogFC,
           y = HighHealthy_LogFC,
           colour = comparison,
           text = paste0("CpG: ", CpG, "\n",
                         "Th2-high vs H - pAdj ", HighHealthy_adjPValue, "\n",
                         "Th2-low vs H - pAdj ", LowHealthy_adjPValue))) + 
  # scale_x_continuous(limits = c(-1.5, 1.5)) + scale_y_continuous(limits = c(-1.5, 1.5)) +
  geom_hline(yintercept = 0, linetype = "solid") + geom_vline(xintercept = 0,  linetype = "solid") + 
  geom_abline(slope = 1, linetype = "dashed") + geom_abline(slope = -1, linetype = "dashed") +
  geom_point(size = 3, alpha = 0.4) + 
  scale_color_manual(values = c("#1f77b4", "#e6550d", "#ffb74d")) +
  guides(colour = guide_legend(title="Comparison")) + 
  labs(x = "Th2-Low vs Healthy LogFC", y = "Th2-High vs Healthy LogFC",
       caption = "Dashed lines are for visualisation purposes \n to clearly mark the quadrants and sections") + 
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        legend.position = "bottom")
# ggsave(plot = scatterlogFC, height = 10, width = 10, dpi = 300,
#        filename = "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/006 Differential Methylated Region Analysis DMR/HighvsLowvsHealthy_Comparison2.tiff")

plotly::ggplotly(scatterlogFC)
plotly.scatterlogFC <- ggplotly(scatterlogFC)
# save the widget at .html format
# library("htmlwidgets")
htmlwidgets::saveWidget(as_widget(plotly.scatterlogFC), file=file.path(base.dir, "Results/006 Differential Methylated Region Analysis DMR/scatterLogFC.html"))

# ##############################################################################
# #################### narrowPeak format for ENCODE  ################
# ##############################################################################
# # https://genome.ucsc.edu/cgi-bin/hgCustom?hgsid=2327924484_qwMSBhpnxs3ft4IU6BQh43afajdx
# 
# narrowPeakFormat.DMA <- DMA.result[DMA.result$adj.P.Val < 0.05, ]
# 
# narrowPeakFormat.DMA <- cpg.location[rownames(narrowPeakFormat.DMA), c("chr", "pos", "strand")]
# 
# 
# narrowPeakFormat.DMA2 <- DMA.result[rownames(narrowPeakFormat.DMA), c("P.Value", "adj.P.Val")]
# narrowPeakFormat.DMA2$P.Value <- -log10(narrowPeakFormat.DMA2$P.Value)
# narrowPeakFormat.DMA2$adj.P.Val <- -log10(narrowPeakFormat.DMA2$adj.P.Val)
# 
# narrowPeakFormat.DMA <- cbind(narrowPeakFormat.DMA, rownames(narrowPeakFormat.DMA), narrowPeakFormat.DMA2)
# colnames(narrowPeakFormat.DMA) <- c("chr", "start","strand", "name", "pValue", "qValue")
# narrowPeakFormat.DMA$end <- narrowPeakFormat.DMA$start + 1
# narrowPeakFormat.DMA$score <- 0
# narrowPeakFormat.DMA$signalValue <- 0
# narrowPeakFormat.DMA$peak <- 0
# 
# narrowPeakFormat.DMA <- narrowPeakFormat.DMA[, c("chr", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")]
# # write.table(x = narrowPeakFormat.DMA, file = file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3/DMA_narrowPeakFormat_ENCODE.csv"),
# #           sep = "\t", quote = FALSE)
# 
# eQTM.results <- read.csv(file.path(base.dir, "Data/005 eQTM analysis/EQTM_results_all.csv"),
#                          header = TRUE, row.names = 1)
# eQTM.sites <- unique(eQTM.results$snps)
# narrowPeakFormat.eQTM <- narrowPeakFormat.DMA[eQTM.sites, ]
# # write.table(x = narrowPeakFormat.eQTM, file = file.path(base.dir, "Data/005 eQTM analysis/eQTM_narrowPeakFormat_ENCODE.csv"),
# #           sep = "\t", quote = FALSE)
