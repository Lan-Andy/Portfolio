# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
input.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results"
base.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/"
output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/08_DMR"

##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("dplyr")
library("ggplot2")
library("ggrepel")
library("patchwork")
library("pheatmap")

library("PupillometryR")

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
# Asthma vs healthy   = DMR_results_asthma_vs_control_750gap.csv
# Th2-high vs healthy = DMR_results_Th2HighvsHealthy_750gap.csv
# Th2-low vs healthy  = DMR_results_Th2HighvsTh2Low_750gap.csv
# Th2-high vs Th2-low = DMR_results_Th2LowvsHealthy_750gap.csv

DMR.input <- file.path(input.dir, "Data/006 Differential Methylated Region Analysis DMR")
DMR.results <- read.csv(file = file.path(DMR.input, "DMR_results_asthma_vs_control_750gap.csv"),
                        header = TRUE,
                        row.names = 1)

# Filtering out the regions with a single CpG site
DMR.results <- DMR.results[-which(DMR.results$n==1),]

DMR.results$Legend <- ifelse(
  DMR.results$p.adjust < 0.05,
  ifelse(
    DMR.results$estimate < 0,
    "Hypomethylated",
    "Hypermethylated"
  ),
  "Not Significant"
)


annotation_data <- DMR.results[order(DMR.results$p.value), ]
annotation_data <- head(annotation_data, n = 10)
annotation_data$name <- paste0(annotation_data$chr, ":",annotation_data$start, "-", annotation_data$end)

################################################################################
##################### Volcano plot DMRs ########################################
################################################################################
y.bar.intercept <- DMR.results %>%
  dplyr::filter(
    as.character(Legend) == "Not Significant"
  ) %>%
  dplyr::pull(p.value) %>%
  min() %>%
  log10() * -1

DMR_volcano <- ggplot2::ggplot(
  data = DMR.results,
  mapping = ggplot2::aes(
    x = estimate,
    y = -log10(p.value)
  )
) +
  theme_bw() +
  ggplot2::geom_point(
    mapping = aes(
      color = Legend
    )
  ) +
  ggplot2::ylab("-log<sub>10</sub>(P-value)") +
  ggplot2::xlab("Estimate") +
  ggplot2::scale_color_manual(
    values = c(
      'Hypomethylated' = "blue",
      'Hypermethylated' = "red",
      'Not Significant' = "grey"
    )
  ) +
  ggplot2::geom_hline(
    yintercept = y.bar.intercept,
    colour = "black",
    linetype = "dashed"
  ) +
  ggrepel::geom_label_repel(
    data = annotation_data,
    mapping = ggplot2::aes(
      label = name
    ),
    size = 5,
    max.overlaps = 20,
    min.segment.length = 0
  ) +
  ggplot2::guides(
    color = guide_legend(
      override.aes = list(size = 3)
    )
  ) +
  # ggprism::theme_prism() +
  ggplot2::theme(
    legend.position = "bottom",
    axis.title.x = ggtext::element_markdown(size = 16),
    axis.title.y = ggtext::element_markdown(size = 16),
    axis.text.x = ggplot2::element_text(size = 14),
    axis.text.y = ggplot2::element_text(size = 14),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    plot.caption = ggtext::element_markdown(lineheight = 1.2)
  )
tiff(file.path("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS",
               "/Results/006 Differential Methylated Region Analysis DMR/volcanoplot_sig_CpGs_AsthmaHealthy_FINAL.tiff"), width = 3000, height = 3000, res = 300)
print(DMR_volcano)
dev.off()

ggsave(DMR_volcano, path = output.dir, file="006 Differential Methylated Region Analysis DMR/DMR_volcano_Th2LowvsHealthy.tiff",
       device = "png", width = 18, height = 18, units = "cm")

################################################################################
##################### Manhattan plot DMRs ######################################
################################################################################
# Manhattan for DMRs
# tT <- read.csv(file = file.path(input.dir, "Results/006 Differential Methylated Region Analysis DMR/DMR_results_750gap.csv"),
#                header = TRUE, row.names = 1)

# Asthma vs healthy   = 08_DMR/DMR_results_asthma_vs_control_750gap.csv
# Th2-high vs healthy = 08_DMR/DMR_results_Th2HighvsHealthy_750gap.csv
# Th2-low vs healthy  = 08_DMR/DMR_results_Th2HighvsTh2Low_750gap.csv
# Th2-high vs Th2-low = 08_DMR/DMR_results_Th2LowvsHealthy_750gap.csv

tT <- read.csv(file = file.path(DMR.input, "DMR_results_asthma_vs_control_750gap.csv"),
               header = TRUE, row.names = 1)

tT <- tT[-which(tT$n==1),]
selection <- which(tT$p.adjust<0.05)
tT2 <- tT[selection,]
rownames(tT2) <- paste0(tT2$chr, ":", tT2$start,"-", tT2$end)

# QQ plot and Manhattan plot
data <- tT
# qqPlot(data[, "adj.P.Val"])

SNP = paste0(tT$chr, ":", tT$start,"-", tT$end)
CHR = gsub(pattern = "chr", replacement = "", x = tT$chr)
BP = (tT$start + tT$end)/2
P = tT$p.value

Results_data <- as.data.frame(cbind(SNP, CHR, BP, P))
Results_data$CHR <- as.numeric(Results_data$CHR)
Results_data$BP  <- as.numeric(Results_data$BP)
Results_data$P   <- as.numeric(Results_data$P)

colnames(Results_data) <- c("SNP", "CHR", "BP", "P")


# Define custom function: https://stackoverflow.com/a/6468532
# Round the number up to the nearest table
roundUp <- function(x,to=10)
{
  to*(x%/%to + as.logical(x%%to))
}

########################################
### Change parameters if needed here ###
########################################
# Set annotation p-value threshold
# Amount of CpG you want to have labelled
topSnps <- 20
# Set significance threshold
signThreshold <- 0.05
######################################## 

# Selecting the threshold within the figure which CpGs should be labelled
AnnoThreshold <- -log(max(head(tT[order(tT$p.adjust), "p.value"], n=topSnps)), base = 10)
# Filtering the CpGs that should be highlighted in the figure
snpsOfInterest <- rownames(tT2)
# Nominal suggestive significance line
suggestiveLine <- -log(signThreshold, base = 10)
# Genomewide significance line
genomewideLine <- -log(max(tT[tT$p.adjust<signThreshold, "p.value"]), base = 10)

upperLimit <- roundUp(-log(min(tT[tT$p.adjust<signThreshold, "p.value"]), base = 10), to = 10)

# ggplot2 variant of the Manhattan plot is taken from 
# https://r-graph-gallery.com/101_Manhattan_plot.html Accessed 29/07/2024
# Uses the same data structure as qqman::manhattan
don <- Results_data %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(Results_data, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate(is_annotate=ifelse(-log10(P)>AnnoThreshold, "yes", "no")) 

axisdf = don %>%
  group_by(CHR) %>%
  summarize(center=(max(BPcum) + min(BPcum) ) / 2 )

manhattan_plot <- ggplot(don, aes(x = BPcum, y = -log10(P))) +
  
  # Show all points
  geom_point(data = subset(don, -log10(P) <= genomewideLine), aes(color = as.factor(CHR)), alpha = 0.8, size = 1.3) +
  scale_colour_manual(values = rep(c("grey60", "grey40"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
  scale_y_continuous(limits = c(0,upperLimit), expand = c(0, 0), breaks = seq(0, upperLimit, 5)) +     # remove space between plot area and x axis
  
  # Change axes titles
  labs(x = "Chromosome", y = expression("-Log"[10]*"(P-value)"),
       title = "") +
  
  # Add highlighted points
  geom_point(data = subset(don, -log10(P) > genomewideLine), 
             colour = "#0ca6f3", alpha = 0.8, size = 1.3) +
  # geom_point(data = subset(don, -log10(P) > genomewideLine), 
  #            color = "firebrick1", size = 1.3, pch = 16) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel(data = subset(don, is_annotate == "yes"), aes(label = SNP), size = 2) +
  
  # Add suggestive and genome-wide significance line
  # geom_hline(yintercept = suggestiveLine, colour = "blue") + 
  geom_hline(yintercept = genomewideLine, colour = "black", linetype = "dashed") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
# output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/08_DMR/"

ggsave(plot = manhattan_plot, 
       filename = file.path(output.dir, "006 Differential Methylated Region Analysis DMR/Manhattan_DMR_Th2LowvsHealthy.tiff"),
       dpi = 300)

################################################################################
##################### Filter DMRs to create boxplots for 370 samples ###########
################################################################################
input.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
meth.data <- read.csv(file = file.path(input.dir, "Data/006 Differential Methylated Region Analysis DMR/DMR_mval_Th2LowvsHealthy_allSamples.csv"),
                      header = TRUE,
                      row.names = 1)
DMR.results <- read.csv(file = file.path(input.dir, "Data/006 Differential Methylated Region Analysis DMR/DMR_results_Th2LowvsHealthy_750gap.csv"),
                        header = TRUE, row.names = 1)
Th2.groups <- read.csv(file = file.path(input.dir, "Data/007 DMR eQTM analysis/covariates_eQTM_samples_IDs.csv"),
                       header = TRUE, row.names = 1)
linkingtable <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/methID_rnaID_clinID.csv")
pheno.data <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/atlantis_patient_data.csv",
                       header = TRUE, row.names = 1, na.strings = "")

colnames(meth.data) <- gsub(pattern = "X", replacement = "", x = colnames(meth.data))

linkingtable$VISIT <- gsub(pattern = "Visit 1a", replacement = "VISIT 1", x = linkingtable$VISIT)
linkingtable$VISIT <- gsub(pattern = "Visit 2", replacement = "VISIT 2", x = linkingtable$VISIT)

linktable1 <- linkingtable[!is.na(linkingtable$PT) & linkingtable$VISIT == "VISIT 1",]

DMR.results <- DMR.results[-which(DMR.results$n==1),]
DMR.results <- DMR.results[DMR.results$p.adjust < 0.05, ]
DMR.results <- DMR.results[order(DMR.results$p.adjust), ]
DMR.results$name <- paste0(DMR.results$chr, ":",DMR.results$start, "-", DMR.results$end)

# Check the colnames of methylation samples is in the link table
methylation.patients <- base::intersect(colnames(meth.data), linktable1$meth_file_id)
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
    }
    else{
      # print("There is VISIT 1")
    }
  }
  rm(temp, patient.id)
}

linktable3 <- rbind(linktable2, linkingtable[linkingtable$PT %in% single.visit2.ids$PT, ])

pd.methyl <- pheno.data[pheno.data$PT %in% linktable2$PT, ]
pd.methyl.v1 <- pd.methyl[pd.methyl$VISIT %in% "VISIT 1", ]

pd.methyl <- pheno.data[pheno.data$PT %in% single.visit2.ids$PT, ]
pd.methyl <- pd.methyl %>% tidyr::fill(SMOKE)
pd.methyl.v2 <- pd.methyl[pd.methyl$VISIT %in% "VISIT 2", ]

pd.methyl.sub <- rbind(pd.methyl.v1, pd.methyl.v2)

# PC.all.cp <- PC.all.cp[PC.all.cp$meth_file_id %in% linktable3$meth_file_id, ]

pd.methyl.sub <- merge(pd.methyl.sub, linktable3, by = "PT")
# pd.methyl.sub <- merge(pd.methyl.sub, PC.all.cp, by = "meth_file_id")

pd.methyl.sub$SMOKE <- gsub(pattern = "-| ", replacement = "_", x = pd.methyl.sub$SMOKE)
pd.methyl.sub$SMOKE <- gsub(pattern = "Non", replacement = "Never", x = pd.methyl.sub$SMOKE)

pheno.data.Th2 <- pd.methyl.sub%>%
  dplyr::filter(ASTHEA == 'A')%>%
  # dplyr::mutate(include = if_else(((SYS_COR == 'Yes') | (BIO == 'Yes')), 'NO', 'YES'))%>%
  dplyr::mutate(include = if_else((SYS_COR == 'Yes'), 'NO', 'YES'))%>%
  dplyr::mutate(group_th = case_when(
    ((LABEOSV > 0.3) & (FENRES > 25)) ~ 'high',
    ((LABEOSV < 0.15) & (FENRES < 25)) ~ 'low',
    ((LABEOSV > 0.3) & (is.na(FENRES))) ~ 'high',
    ((LABEOSV < 0.15) & (is.na(FENRES))) ~ 'low')) %>%
  dplyr::mutate(group_th = dplyr::if_else(is.na(group_th), 'undeterm', group_th)) %>%
  dplyr::mutate(group_th = dplyr::if_else(include == 'NO', 'cor_bio', group_th))

# Adding the Th2 group into the main pheno table
#   to be able to compare the Th2-high/low to healthy participants
pd.methyl.sub$group_th <- "healthy" 
# Assign the Th2 group towards the right sample
pd.methyl.sub[which(pd.methyl.sub$meth_file_id %in% pheno.data.Th2$meth_file_id), "group_th"] <- pheno.data.Th2$group_th

# For creating better legends in figures? 
covariates <- pd.methyl.sub
covariates$ASTHEA <- ifelse(covariates$ASTHEA == "A", "Asthma", "Healthy")
covariates$SMOKE <- ifelse(covariates$SMOKE == "Never_smoker", 0,
                           ifelse(covariates$SMOKE == "Ex_smoker", 1,
                                  ifelse(covariates$SMOKE == "Current_Smoker", 2, NA)))
covariates$SEX <- ifelse(covariates$SEX == "F", 0, 1)
covariates$group_th <- ifelse(covariates$group_th %in% c("undeterm", "cor_bio"),
                              "undeterm",
                              covariates$group_th)

meth.data <- meth.data[, pd.methyl.sub$meth_file_id]

incl.endotypes = TRUE
plot.data3 <- as.data.frame(t(meth.data[, covariates$meth_file_id]))
if (incl.endotypes){
  plot.data3$status <- covariates$group_th
  plot.data3$status <- factor(plot.data3$status, levels = c("healthy", "low", "undeterm", "high"))
  
  plot.data3$status <- gsub(pattern = "healthy", replacement = "Healthy", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "high", replacement = "Th2-High", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "low", replacement = "Th2-Low", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "undeterm", replacement = "Th2-intermediate", x = plot.data3$status)
  plot.data3$status <- factor(plot.data3$status, levels = c("Healthy", "Th2-Low", "Th2-intermediate", "Th2-High"))
} else {
  plot.data3$status <- Th2.groups$ASTHEA
  plot.data3$status <- factor(plot.data3$status, levels = c("H", "A"))
  
  plot.data3$status <- gsub(pattern = "H", replacement = "Healthy", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "A", replacement = "Asthma", x = plot.data3$status)
  plot.data3$status <- factor(plot.data3$status, levels = c("Healthy", "Asthma"))
}  

dmr_columns <- colnames(plot.data3)[colnames(plot.data3) != "status"]

# Initialize a result data frame to store means for each group across all DMRs
result <- data.frame()

# Loop over each DMR column to calculate group means
for (dmr in dmr_columns) {
  # Calculate the mean methylation for each group
  group_means <- aggregate(
    plot.data3[[dmr]] ~ plot.data3$status,
    data = plot.data3,
    FUN = mean
  )
  # Rename columns for clarity
  colnames(group_means) <- c("Group", "Mean_Methylation")
  # Add the DMR name as a new column
  group_means$DMR <- dmr
  # Append to the result data frame
  result <- rbind(result, group_means)
}
# Reorganize the result for readability
result <- result[, c("DMR", "Group", "Mean_Methylation")]

# Create a new dataframe to store results
difference_results <- data.frame()

# Loop through each unique DMR in the result dataframe
for (dmr in unique(result$DMR)) {
  # Filter the data for the current DMR
  result2 <- result[result$DMR == dmr, ]
  
  # Extract the mean methylation for the relevant groups
  p_value <- DMR.results[DMR.results$name == dmr, "p.adjust"]
  mean_healthy <- result2[result2$Group == "Healthy", "Mean_Methylation"]
  mean_th2_low <- result2[result2$Group == "Th2-Low", "Mean_Methylation"]
  mean_th2_high <- result2[result2$Group == "Th2-High", "Mean_Methylation"]
  
  # Calculate the desired differences
  th2_low_minus_healthy <- mean_th2_low - mean_healthy
  th2_low_minus_th2_high <- mean_th2_low - mean_th2_high
  difference <- th2_low_minus_healthy - th2_low_minus_th2_high
  
  # Append the results to the new dataframe
  difference_results <- rbind(difference_results, data.frame(
    DMR = dmr,
    Th2_Low_minus_Healthy = th2_low_minus_healthy,
    Th2_Low_minus_Th2_High = th2_low_minus_th2_high,
    Difference = difference,
    adjust.pvalue = p_value
  ))
}
difference_results <- difference_results[order(abs(difference_results$Difference)), ]
filtered_results <- difference_results[difference_results$adjust.pvalue <= 10E-10, ]
filtered_results <- filtered_results[abs(filtered_results$Th2_Low_minus_Healthy) >= 0.02 
                                     & abs(filtered_results$Difference) <= 0.02
                                     , ]

plot.data <- plot.data3
top.cpgs <- "chr3:24494165-24495916" #DMR.results$name[1]

COLOR_TEMP_ALL <- c("#1f77b4", "#ffb74d", "black", "#e6550d")
names(COLOR_TEMP_ALL) <- c("Healthy", "Th2-Low", "Th2-intermediate", "Th2-High")
tiff(file.path("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS",
               "/Results/006 Differential Methylated Region Analysis DMR/topDMR_Th2_high.tiff"), width = 1000, height = 1000, res = 300)
for (cpg in top.cpgs) {
  # Wrap column name in backticks to handle special characters
  cpg_safe <- paste0("`", cpg, "`")
  
  
  # Create the boxplot
  p <- ggplot(plot.data, aes_string(x = "status", y = cpg_safe, fill = "status")) +
    geom_flat_violin(position = position_nudge(x = 0.2),
                     alpha = 0.4) +
    geom_point(aes(colour = status),
               position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.25),
               size = 1,
               alpha = 0.4,
               show.legend = F) +
    geom_boxplot(width = .25,
                 outlier.shape = NA,
                 alpha = 0.5) +
    scale_fill_manual(values = COLOR_TEMP_ALL) + # Apply custom colors
    scale_colour_manual(values = COLOR_TEMP_ALL) + # Apply custom colors
    labs(
      title = paste("Boxplot for", cpg),
      x = "Group",
      y = "Methylation (beta)"
    ) +
    theme_bw() +
    theme(legend.position = "none", # Optionally hide legend if redundant
          panel.grid = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14), 
          axis.text.x = element_text(size = 8),
          axis.title.x = element_blank()) 
  
  # Print the plot to display
  print(p)
  
  # Optionally save the plot to file
  # ggsave(filename = paste0("boxplot_", gsub("[:/]", "_", cpg), ".png"), plot = p)
}
dev.off()

################################################################################
##################### Boxplots of CpGs in DMR ##################################
################################################################################

# Asthma vs healthy   = /DMR_results_asthma_vs_control_750gap.csv
# Asthma vs healthy   = /DMR_identification_AsthmaVsHealthy.csv

# Th2-high vs healthy = /DMR_results_Th2HighvsHealthy_750gap.csv
# Th2-high vs healthy = /DMR_identification_Th2HighvsHealthy.csv

# Th2-low vs healthy  = /DMR_results_Th2LowvsHealthy_750gap.csv
# Th2-low vs healthy  = /DMR_identification_Th2LowvsHealthy.csv

# Th2-high vs Th2-low = /DMR_results_Th2HighvsTh2Low_750gap.csv
# Th2-high vs Th2-low = /DMR_identification_Th2HighvsTh2Low.csv
input.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"

DMR.input <- file.path(input.dir, "Data/006 Differential Methylated Region Analysis DMR")
DMR.results <- read.csv(file = file.path(DMR.input, "DMR_results_Th2LowvsHealthy_750gap.csv"),
               header = TRUE, row.names = 1)
DMR.ident <- read.csv(file = file.path(DMR.input, "DMR_identification_Th2LowvsHealthy.csv"),
                        header = TRUE, row.names = 1)
meth.data <- read.csv(file.path(DMR.input, "DMR_all_CpGs_betas.csv"),
                      header = TRUE, row.names = 1)
DMR.data <- read.csv(file.path(input.dir, "Data/006 Differential Methylated Region Analysis DMR/DMR_mval_Th2LowvsHealthy_allSamples.csv"),
                     header = TRUE, row.names = 1)
colnames(meth.data) <- gsub(pattern = "X", replacement = "", x = colnames(meth.data))
colnames(DMR.data) <- gsub(pattern = "X", replacement = "", x = colnames(DMR.data))

DMR.ident <- DMR.ident[order(DMR.ident$DMR, DMR.ident$pos), ]


source("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Documentation/Github/004 Bioinformatics/Scripts/M2beta.R")
DMR.data <- M2beta(DMR.data)
# Filtering out the regions with a single CpG site
DMR.results <- DMR.results[-which(DMR.results$n==1),]
DMR.results <- DMR.results[DMR.results$p.adjust < 0.05, ]
DMR.results <- DMR.results[order(DMR.results$p.adjust), ]
DMR.results$name <- paste0(DMR.results$chr, ":",DMR.results$start, "-", DMR.results$end)

# Get the name of the top DMR
# top.DMR <- unique(DMR.ident[DMR.ident$DMR_pAdjust == min(DMR.ident$DMR_pAdjust), "DMR"])

top.DMR <- unique(DMR.ident[DMR.ident$DMR_pAdjust == unique(head((DMR.ident$DMR_pAdjust))), "DMR"])
top.DMR <- unique(DMR.ident$DMR)[1:20]


top.DMR <- "chr17:28717304-28718284"
top.DMR <- "chr3:24494165-24495916"           
top.DMR <- DMR.results$name


# Getting the CpG sites from the top DMR
CpG.sites.top.DMR <- rownames(DMR.ident[DMR.ident$DMR %in% top.DMR, ])
# Ordering the CpG sites based on physical chr:pos (increasing)
CpG.sites.top.DMR <- CpG.sites.top.DMR[order(DMR.ident[CpG.sites.top.DMR, "pos"])]


# plot.data <- as.data.frame(t(meth.data[CpG.sites.top.DMR, Th2.groups$meth_file_id]))
# plot.data$status <- Th2.groups$ASTHEA
# plot.data$status <- gsub(pattern = "A", replacement = "Asthma", x = plot.data$status)
# plot.data$status <- gsub(pattern = "H", replacement = "Healthy", x = plot.data$status)
# plot.data$status <- as.factor(plot.data$status)


# Asthma vs healthy 
# Don't change anything, all samples are taken

# Th2-high vs healthy
# Th2.groups <- Th2.groups[grep(pattern = "high|healthy", x = Th2.groups$group_th), ]

# Th2-low vs healthy
Th2.groups <- Th2.groups[grep(pattern = "low|healthy", x = Th2.groups$group_th), ]

# Th2-high vs Th2-low
# Th2.groups <- Th2.groups[grep(pattern = "high|low", x = Th2.groups$group_th), ]

identification <- data.frame(
  "CpG" = rownames(DMR.ident),
  "DMR" = DMR.ident$DMR
)

## Am I looking at Asthma vs Healthy or T2-groups vs Healthy
incl.endotypes = TRUE

plot.data <- as.data.frame(t(meth.data[CpG.sites.top.DMR, Th2.groups$meth_file_id]))
## Missing some CpG sites?? 
## rownames(meth.data) == "cg23626733"
if (incl.endotypes){
  plot.data$status <- Th2.groups$group_th
  plot.data$status <- factor(plot.data$status, levels = c("healthy", "low", "undeterm", "high"))
  
  plot.data$status <- gsub(pattern = "healthy", replacement = "Healthy", x = plot.data$status)
  plot.data$status <- gsub(pattern = "high", replacement = "Th2-High", x = plot.data$status)
  plot.data$status <- gsub(pattern = "low", replacement = "Th2-Low", x = plot.data$status)
  plot.data$status <- gsub(pattern = "undeterm", replacement = "Th2-intermediate", x = plot.data$status)
  # plot.data <- plot.data[-(grep(pattern = "Th2-High|Th2-intermediate", x = plot.data$status)), ]
  plot.data$status <- factor(plot.data$status, levels = c("Healthy", "Th2-Low", "Th2-intermediate", "Th2-High"))
} else {
  plot.data$status <- Th2.groups$ASTHEA
  plot.data$status <- factor(plot.data$status, levels = c("H", "A"))
  
  plot.data$status <- gsub(pattern = "H", replacement = "Healthy", x = plot.data$status)
  plot.data$status <- gsub(pattern = "A", replacement = "Asthma", x = plot.data$status)
  plot.data$status <- factor(plot.data$status, levels = c("Healthy", "Asthma"))
}  
  



plot.data <- as.data.frame(t(plot.data))
plot.data$CpG <- rownames(plot.data)
plot.data <- merge(plot.data, identification, by.x = "CpG", by.y = "CpG", all.x = TRUE)
plot.data <- plot.data[match(identification$CpG, plot.data$CpG), ]

rownames(plot.data) <- plot.data$CpG

plot.data2 <- reshape2::melt(plot.data)

# Boxplots per asthmatic group as a whole and not divided in individual CpG sites
plot.data3 <- as.data.frame(t(DMR.data[top.DMR, Th2.groups$meth_file_id]))
if (incl.endotypes){
  plot.data3$status <- Th2.groups$group_th
  plot.data3$status <- factor(plot.data3$status, levels = c("healthy", "low", "undeterm", "high"))
  
  plot.data3$status <- gsub(pattern = "healthy", replacement = "Healthy", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "high", replacement = "Th2-High", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "low", replacement = "Th2-Low", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "undeterm", replacement = "Th2-intermediate", x = plot.data3$status)
  plot.data3$status <- factor(plot.data3$status, levels = c("Healthy", "Th2-Low", "Th2-intermediate", "Th2-High"))
} else {
  plot.data3$status <- Th2.groups$ASTHEA
  plot.data3$status <- factor(plot.data3$status, levels = c("H", "A"))
  
  plot.data3$status <- gsub(pattern = "H", replacement = "Healthy", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "A", replacement = "Asthma", x = plot.data3$status)
  plot.data3$status <- factor(plot.data3$status, levels = c("Healthy", "Asthma"))
}  

dmr_columns <- colnames(plot.data3)[colnames(plot.data3) != "status"]

# Initialize a result data frame to store means for each group across all DMRs
result <- data.frame()

# Loop over each DMR column to calculate group means
for (dmr in dmr_columns) {
  # Calculate the mean methylation for each group
  group_means <- aggregate(
    plot.data3[[dmr]] ~ plot.data3$status,
    data = plot.data3,
    FUN = mean
  )
  # Rename columns for clarity
  colnames(group_means) <- c("Group", "Mean_Methylation")
  # Add the DMR name as a new column
  group_means$DMR <- dmr
  # Append to the result data frame
  result <- rbind(result, group_means)
}
# Reorganize the result for readability
result <- result[, c("DMR", "Group", "Mean_Methylation")]

# Create a new dataframe to store results
difference_results <- data.frame()

# Loop through each unique DMR in the result dataframe
for (dmr in unique(result$DMR)) {
  # Filter the data for the current DMR
  result2 <- result[result$DMR == dmr, ]
  
  # Extract the mean methylation for the relevant groups
  p_value <- DMR.results[DMR.results$name == dmr, "p.adjust"]
  mean_healthy <- result2[result2$Group == "Healthy", "Mean_Methylation"]
  mean_th2_low <- result2[result2$Group == "Th2-Low", "Mean_Methylation"]
  mean_th2_high <- result2[result2$Group == "Th2-High", "Mean_Methylation"]
  
  # Calculate the desired differences
  th2_low_minus_healthy <- mean_th2_low - mean_healthy
  th2_low_minus_th2_high <- mean_th2_low - mean_th2_high
  difference <- th2_low_minus_healthy - th2_low_minus_th2_high
  
  # Append the results to the new dataframe
  difference_results <- rbind(difference_results, data.frame(
    DMR = dmr,
    Th2_Low_minus_Healthy = th2_low_minus_healthy,
    Th2_Low_minus_Th2_High = th2_low_minus_th2_high,
    Difference = difference,
    adjust.pvalue = p_value
  ))
}
difference_results <- difference_results[order(abs(difference_results$Difference)), ]
filtered_results <- difference_results[difference_results$adjust.pvalue <= 10E-10, ]
filtered_results <- filtered_results[abs(filtered_results$Th2_Low_minus_Healthy) >= 0.02 &
                                       abs(filtered_results$Difference) <= 0.02, ]

top.DMR <- filtered_results$DMR


COLOR_TEMP_ALL <- c("#1f77b4", "#ffb74d", "black", "#e6550d")
names(COLOR_TEMP_ALL) <- c("Healthy", "Th2-Low", "Th2-intermediate", "Th2-High")

output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/"
# pdf(file = file.path(output.dir, "/006 Differential Methylated Region Analysis DMR/ALL_DMR_Th2LowvsHealthy.pdf"),
#     width = 10, height = 8.5)
# Iterate through each DMR column and create a boxplot
for (dmr in filtered_results$DMR) {
  # Wrap column name in backticks to handle special characters
  dmr_safe <- paste0("`", dmr, "`")
  
  # Create the boxplot
  p <- ggplot(plot.data3, aes_string(x = "status", y = dmr_safe, fill = "status")) +
    geom_flat_violin(position = position_nudge(x = 0.2),
                     alpha = 0.4) +
    geom_point(aes(colour = status),
               position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.25),
               size = 1,
               alpha = 0.4,
               show.legend = F) +
    geom_boxplot(width = .25,
                 outlier.shape = NA,
                 alpha = 0.5) +
    scale_fill_manual(values = COLOR_TEMP_ALL) + # Apply custom colors
    scale_colour_manual(values = COLOR_TEMP_ALL) + # Apply custom colors
    labs(
      title = paste("Boxplot for", dmr),
      x = "Group",
      y = "Methylation (beta)"
    ) +
    theme_bw() +
    theme(legend.position = "none", # Optionally hide legend if redundant
          panel.grid = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 14),
          axis.title.x = element_blank()) 
  
  # Print the plot to display
  print(p)
  
  # Optionally save the plot to file
  # ggsave(filename = paste0("boxplot_", gsub("[:/]", "_", dmr), ".png"), plot = p)
}
dev.off()

plot.data4 <- reshape2::melt(plot.data3)

unique_variables <- unique(plot.data4$variable)[1:10]




global_min <- min(c(min(plot.data2$value, na.rm = TRUE), min(plot.data4$value, na.rm = TRUE)))
global_max <- max(c(max(plot.data2$value, na.rm = TRUE), max(plot.data4$value, na.rm = TRUE)))

# # Asthma  Healthy
COLOR_TEMP = c("#1f77b4","#ff7f0e")
# # Healthy   Th2-high
# COLOR_TEMP = c("#1f77b4", "#e6550d")
# # Healthy   Th2-Low
# COLOR_TEMP = c("#1f77b4", "#ffb74d")
# # Th2-high  Th2-low    
# COLOR_TEMP = c("#e6550d", "#ffb74d")
# Healthy   Th2-high  Th2-low   Undetermined
COLOR_TEMP_ALL = c("#1f77b4", "#ffb74d", "black", "#e6550d")

pdf(file = file.path(output.dir, "/006 Differential Methylated Region Analysis DMR/Top_DMR_CpGs_Th2LowvsHealthy.pdf"),
width = 10, height = 8.5)
cpgs.plot <- ggplot(data = plot.data2,
       aes(x = variable, y = value)) + 
  labs(title = top.DMR,
       y = "Beta value [normalized]",
       colour = "Status",
       fill = "Status") + 
  geom_boxplot(aes(fill = status), outliers = FALSE, alpha = 0.666) + 
  geom_point(aes(colour = status), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.7),
             alpha = 0.1) + 
  scale_fill_manual(values = COLOR_TEMP_ALL) +
  scale_colour_manual(values = COLOR_TEMP_ALL) +
  scale_y_continuous(limits = c(global_min, global_max)) +  # Set y-axis limits dynamically
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size =14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none")

overall.plot <- ggplot(data = plot.data4,
       aes(x = variable, y = value)) + 
  labs(title = top.DMR,
       y = "Beta value [normalized]",
       colour = "Status",
       fill = "Status") + 
  geom_boxplot(aes(fill = status), outliers = FALSE, alpha = 0.666) + 
  geom_point(aes(colour = status), 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.7),
             alpha = 0.1) + 
  scale_fill_manual(values = COLOR_TEMP_ALL) +
  scale_colour_manual(values = COLOR_TEMP_ALL) +
  scale_y_continuous(limits = c(global_min, global_max)) +  # Set y-axis limits dynamically
  theme_bw() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

layout <- "22222222211"
p_all <- overall.plot + cpgs.plot + 
  plot_layout(design = layout) 
tiff(file.path("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS",
               "/Results/006 Differential Methylated Region Analysis DMR/topDMR_sig_CpGs_Th2_high_healthy_FINAL.tiff"), width = 3000, height = 3000, res = 300)
print(p_all)
dev.off()


# Set colours for groups
COLOR_TEMP_ALL <- c("#1f77b4", "#ffb74d", "black", "#e6550d")

# Create a list to store plots
plot_list <- list()

# Get unique DMRs
unique_dmrs <- unique(plot.data4$variable)
for (dmr in unique_dmrs) {
  # Subset data for the current DMR
  dmr_data <- subset(plot.data4, variable == dmr)
  
  # Create a boxplot
  p <- ggplot(data = dmr_data, aes(x = status, y = value, fill = status)) +
    geom_flat_violin(position = position_nudge(x = .20), 
                     alpha = .8) +
    geom_point(aes(color = status), 
               position = position_jitter(width = .10),
               size = .3, 
               alpha = .5,
               show.legend = F) +
    geom_boxplot(width = .3, 
                 outlier.shape = NA,
                 alpha = .5) +
    # geom_boxplot(outlier.shape = NA, alpha = 0.7) + # Boxplot without outliers
    # geom_jitter(position = position_jitter(width = 0.2), size = 1.2, alpha = 0.5) + # Jittered points
    scale_fill_manual(values = COLOR_TEMP_ALL) + # Use custom colours
    scale_colour_manual(values = COLOR_TEMP_ALL) + # Use custom colours
    labs(title = dmr, x = "Status", y = "Methylation [beta]") + # Set titles
    theme_bw() + 
    theme(
      # axis.text.x = element_text(), # Rotate x-axis text
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.x = element_blank(), 
      legend.position = "none",
      panel.grid = element_blank(),
      
    )
  
  # Append the plot to the list
  plot_list[[dmr]] <- p
}

# Define plots per figure
plots_per_figure <- 4  # Maximum of 4 boxplots per figure
num_figures <- ceiling(length(plot_list) / plots_per_figure)

# Iterate through and create figures
for (fig_idx in seq_len(num_figures)) {
  # Extract the plots for the current figure
  start_idx <- (fig_idx - 1) * plots_per_figure + 1
  end_idx <- min(fig_idx * plots_per_figure, length(plot_list))
  current_plots <- plot_list[start_idx:end_idx]
  
  # Combine the current plots into a single figure
  combined_plot <- wrap_plots(current_plots, ncol = 2, nrow = 2) # Arrange in 2x2 grid
  
  # # Save the current figure
  # ggsave(
  #   filename = file.path(paste0("figure_", fig_idx, ".tiff")),
  #   plot = combined_plot,
  #   width = 12, height = 12
  # )
  
  tiff(paste0("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS",
                 "/Results/006 Differential Methylated Region Analysis DMR/Th2_low_DMR_unique_", fig_idx, ".tiff"), width = 3000, height = 3000, res = 300)
  
  # Print the current combined plot to the console (optional)
  print(combined_plot)
  dev.off()
}



################################################################################
##################### Calculate biggest gap between groups #####################
################################################################################

find_top_diff_DMRs <- function(DMR.data, Th2.groups, group1, group2, top_n = 1) {
  # Ensure Th2.groups only contains the relevant groups
  Th2.groups <- Th2.groups[Th2.groups$group_th %in% c(group1, group2), ]
  
  # Find common samples
  common_samples <- intersect(colnames(DMR.data), Th2.groups$meth_file_id)
  
  # Filter DMR.data and Th2.groups
  DMR.data_filtered <- DMR.data[, common_samples]
  Th2.groups_filtered <- Th2.groups[Th2.groups$meth_file_id %in% common_samples, ]
  
  # Ensure the order is the same
  Th2.groups_filtered <- Th2.groups_filtered[match(colnames(DMR.data_filtered), Th2.groups_filtered$meth_file_id), ]
  
  # Verify the matching
  if (!all(colnames(DMR.data_filtered) == Th2.groups_filtered$meth_file_id)) {
    stop("Sample mismatch between DMR.data and Th2.groups")
  }
  
  # Calculate means for each group
  DMR.means <- data.frame(
    group1 = rowMeans(DMR.data_filtered[, Th2.groups_filtered$group_th == group1], na.rm = TRUE),
    group2 = rowMeans(DMR.data_filtered[, Th2.groups_filtered$group_th == group2], na.rm = TRUE)
  )
  
  # Calculate the difference between groups
  DMR.diff <- DMR.means %>%
    mutate(diff = group2 - group1,
           abs_diff = abs(diff))
  
  # Sort by absolute difference and get top N
  top_DMRs <- DMR.diff %>%
    arrange(desc(abs_diff)) %>%
    head(top_n)
  
  # Return results as a list
  return(list(
    top_DMRs = top_DMRs,
    group1 = group1,
    group2 = group2
  ))
}

# For comparing Th2-low vs Th2-high, get top 5 DMRs
result_low_high <- find_top_diff_DMRs(DMR.data, Th2.groups, "low", "high", top_n = 20)

# Print results
cat("Top", nrow(result_low_high$top_DMRs), "DMRs with largest methylation differences between", 
    result_low_high$group1, "and", result_low_high$group2, ":\n")
print(result_low_high$top_DMRs)

# If you want to compare other groups and get only the top DMR:
result_healthy_high <- find_top_diff_DMRs(DMR.data, Th2.groups, "healthy", "high", top_n = 1)

# Print results
cat("Top DMR with largest methylation difference between", 
    result_healthy_high$group1, "and", result_healthy_high$group2, ":\n")
print(result_healthy_high$top_DMRs)


top.DMRs <- rownames(result_low_high$top_DMRs)


pdf(file = file.path(output.dir, "/006 Differential Methylated Region Analysis DMR/Top_DMR_CpGs_Th2LowvsHealthy.pdf"),
    width = 10, height = 8.5)
for (i in 1:length(top.DMRs)) {
  top.DMR <- top.DMRs[i]
  
  # Get CpG sites for the top DMR
  CpG.sites.top.DMR <- rownames(DMR.ident[DMR.ident$DMR == top.DMR, ])
  CpG.sites.top.DMR <- CpG.sites.top.DMR[order(DMR.ident[CpG.sites.top.DMR, "pos"])]
  
  # Prepare plot data for individual CpG sites
  plot.data <- as.data.frame(t(meth.data[CpG.sites.top.DMR, Th2.groups$meth_file_id]))
  plot.data$status <- Th2.groups$group_th
  plot.data$status <- factor(plot.data$status, levels = c("healthy", "low", "undeterm", "high"))
  
  plot.data$status <- gsub(pattern = "healthy", replacement = "Healthy", x = plot.data$status)
  plot.data$status <- gsub(pattern = "high", replacement = "Th2-High", x = plot.data$status)
  plot.data$status <- gsub(pattern = "low", replacement = "Th2-Low", x = plot.data$status)
  plot.data$status <- gsub(pattern = "undeterm", replacement = "Undetermined subtype", x = plot.data$status)
  plot.data <- plot.data[-(grep(pattern = "Th2-High|Undetermined subtype", x = plot.data$status)), ]
  plot.data$status <- factor(plot.data$status, levels = c("Healthy", "Th2-Low", "Undetermined subtype", "Th2-High"))
  
  plot.data2 <- reshape2::melt(plot.data)
  
  # Prepare plot data for overall DMR
  plot.data3 <- as.data.frame(t(DMR.data[top.DMR, Th2.groups$meth_file_id]))
  plot.data3$status <- Th2.groups$group_th
  plot.data3$status <- factor(plot.data3$status, levels = c("healthy", "low", "undeterm", "high"))
  
  plot.data3$status <- gsub(pattern = "healthy", replacement = "Healthy", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "high", replacement = "Th2-High", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "low", replacement = "Th2-Low", x = plot.data3$status)
  plot.data3$status <- gsub(pattern = "undeterm", replacement = "Undetermined subtype", x = plot.data3$status)
  plot.data3$status <- factor(plot.data3$status, levels = c("Healthy", "Th2-Low", "Undetermined subtype", "Th2-High"))
  
  plot.data4 <- reshape2::melt(plot.data3)
  
  global_min <- min(c(min(plot.data2$value, na.rm = TRUE), min(plot.data4$value, na.rm = TRUE)))
  global_max <- max(c(max(plot.data2$value, na.rm = TRUE), max(plot.data4$value, na.rm = TRUE)))
  
  COLOR_TEMP <- c("#1f77b4", "#ffb74d")
  COLOR_TEMP_ALL <- c("#1f77b4", "#ffb74d", "black", "#e6550d")
  
  cpgs.plot <- ggplot(data = plot.data2,
                      aes(x = variable, y = value)) + 
    labs(title = top.DMR,
         y = "Beta value [normalized]",
         colour = "Status",
         fill = "Status") + 
    geom_boxplot(aes(fill = status), outliers = FALSE, alpha = 0.666) + 
    geom_point(aes(colour = status), 
               position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.7),
               alpha = 0.1) + 
    scale_fill_manual(values = COLOR_TEMP) +
    scale_colour_manual(values = COLOR_TEMP) +
    scale_y_continuous(limits = c(global_min, global_max)) +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none")
  
  overall.plot <- ggplot(data = plot.data4,
                         aes(x = variable, y = value)) + 
    labs(title = top.DMR,
         y = "Beta value [normalized]",
         colour = "Status",
         fill = "Status") + 
    geom_boxplot(aes(fill = status), outliers = FALSE, alpha = 0.666) + 
    geom_point(aes(colour = status), 
               position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.7),
               alpha = 0.1) + 
    scale_fill_manual(values = COLOR_TEMP_ALL) +
    scale_colour_manual(values = COLOR_TEMP_ALL) +
    scale_y_continuous(limits = c(global_min, global_max)) +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  layout <- "22222222211"
  p_all <- overall.plot + cpgs.plot + 
    plot_layout(design = layout) 
  
  # Save the plot
  # plot_file <- file.path(output_dir, paste0("Top_DMR_CpGs_", group1, "vs", group2, "_", i, ".pdf"))
  # pdf(file = plot_file, width = 10, height = 8.5)
  print(p_all)
  # dev.off()
  # cat("Plot saved to:", plot_file, "\n")
}
dev.off()


################################################################################
##################### Heatmap with pheatmap #####################
################################################################################

input.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data"

DMR.results <- read.csv(file = file.path(input.dir, "006 Differential Methylated Region Analysis DMR/DMR_results_asthma_vs_control_750gap.csv"),
                        header = TRUE,
                        row.names = 1)
meth.data <- read.csv(file = file.path(input.dir, "006 Differential Methylated Region Analysis DMR/DMR_mval_AsthmaVsHealthy_allSamples.csv"),
                        header = TRUE,
                        row.names = 1)

linkingtable <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/methID_rnaID_clinID.csv")
pheno.data <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/atlantis_patient_data.csv",
                       header = TRUE, row.names = 1, na.strings = "")
PC.all.cp <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/ATLANTIS QC/Data/PCA technical probes/my_pca_data.csv")


colnames(meth.data) <- gsub(pattern = "X", replacement = "", colnames(meth.data))

DMR.results <- DMR.results[-which(DMR.results$n==1),]
DMR.results <- DMR.results[DMR.results$p.adjust < 0.05, ]
DMR.results$region <- paste0(DMR.results$chr,":", DMR.results$start, "-", DMR.results$end)
rownames(DMR.results) <- DMR.results$region

linkingtable$VISIT <- gsub(pattern = "Visit 1a", replacement = "VISIT 1", x = linkingtable$VISIT)
linkingtable$VISIT <- gsub(pattern = "Visit 2", replacement = "VISIT 2", x = linkingtable$VISIT)

linktable1 <- linkingtable[!is.na(linkingtable$PT) & linkingtable$VISIT == "VISIT 1",]

# Check the colnames of methylation samples is in the link table
methylation.patients <- intersect(colnames(meth.data), linktable1$meth_file_id)
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
    }
    else{
      # print("There is VISIT 1")
    }
  }
  rm(temp, patient.id)
}

linktable3 <- rbind(linktable2, linkingtable[linkingtable$PT %in% single.visit2.ids$PT, ])

pd.methyl <- pheno.data[pheno.data$PT %in% linktable2$PT, ]
pd.methyl.v1 <- pd.methyl[pd.methyl$VISIT %in% "VISIT 1", ]

pd.methyl <- pheno.data[pheno.data$PT %in% single.visit2.ids$PT, ]
pd.methyl <- pd.methyl %>% tidyr::fill(SMOKE)
pd.methyl.v2 <- pd.methyl[pd.methyl$VISIT %in% "VISIT 2", ]

pd.methyl.sub <- rbind(pd.methyl.v1, pd.methyl.v2)

PC.all.cp <- PC.all.cp[PC.all.cp$meth_file_id %in% linktable3$meth_file_id, ]

pd.methyl.sub <- merge(pd.methyl.sub, linktable3, by = "PT")
pd.methyl.sub <- merge(pd.methyl.sub, PC.all.cp, by = "meth_file_id")

pd.methyl.sub$SMOKE <- gsub(pattern = "-| ", replacement = "_", x = pd.methyl.sub$SMOKE)
pd.methyl.sub$SMOKE <- gsub(pattern = "Non", replacement = "Never", x = pd.methyl.sub$SMOKE)

# Some changes were made to create the distinction of just T2-high and T2-low
#   These changes were made for the comparison asthma vs healthy
#   Changes include cutt-off at 150. 
#   Not bothering about systemic corticosteroid use
#     Systemic corticosteroid users: 
#         3 individuals < 0.15 10^9 cells/L; 
#         10 individuals >= 0.15 10^9 cells/L
pheno.data.Th2 <- pd.methyl.sub%>%
  dplyr::filter(ASTHEA == 'A')%>%
  # dplyr::mutate(include = if_else(((SYS_COR == 'Yes') | (BIO == 'Yes')), 'NO', 'YES'))%>%
  dplyr::mutate(include = if_else((SYS_COR == 'Yes'), 'NO', 'YES'))%>%
  dplyr::mutate(group_th = case_when(
    # ((LABEOSV > 0.3) & (FENRES > 25)) ~ 'high',
    # ((LABEOSV < 0.15) & (FENRES < 25)) ~ 'low',
    # ((LABEOSV > 0.3) & (is.na(FENRES))) ~ 'high',
    # ((LABEOSV < 0.15) & (is.na(FENRES))) ~ 'low')) %>%
    ((LABEOSV >= 0.15) ) ~ 'high',
    ((LABEOSV < 0.15)) ~ 'low')) #%>%
# dplyr::mutate(group_th = dplyr::if_else(is.na(group_th), 'undeterm', group_th))%>%
# dplyr::mutate(group_th = dplyr::if_else(include == 'NO', 'cor_bio', group_th))

# Adding the Th2 group into the main pheno table
#   to be able to compare the Th2-high/low to healthy participants
pd.methyl.sub$group_th <- "healthy" 
# Assign the Th2 group towards the right sample
pd.methyl.sub[which(pd.methyl.sub$meth_file_id %in% pheno.data.Th2$meth_file_id), "group_th"] <- pheno.data.Th2$group_th

covariates <- pd.methyl.sub
covariates$ASTHEA <- ifelse(covariates$ASTHEA == "A", "Asthma", "Healthy")
covariates$SMOKE <- ifelse(covariates$SMOKE == "Never_smoker", 0,
                           ifelse(covariates$SMOKE == "Ex_smoker", 1, 
                                  ifelse(covariates$SMOKE == "Current_Smoker", 2, NA)))
covariates$SEX <- ifelse(covariates$SEX == "F", 0, 1)
covariates$group_th <- ifelse(covariates$group_th %in% c("undeterm", "cor_bio"), 
                              "undeterm/cor_bio", 
                              covariates$group_th)

meth.data <- meth.data[, pd.methyl.sub$meth_file_id]

heatmap.data <- as.matrix(meth.data[, colnames(meth.data)])
# heatmap.data <- heatmap.data[1:50, ]
# heatmap.data <- as.matrix(meth.data)


health.status <- as.character(covariates$ASTHEA)
health.status <- gsub(pattern = "1", replacement = "Asthma", health.status)
health.status <- gsub(pattern = "0", replacement = "Healthy", health.status)

smoking.status <- as.character(covariates$SMOKE)
smoking.status <- gsub(pattern = "2", replacement = "Current smoker", smoking.status)
smoking.status <- gsub(pattern = "1", replacement = "Ex-smoker", smoking.status)
smoking.status <- gsub(pattern = "0", replacement = "Never-smoker", smoking.status)

sex.status <- as.character(covariates$SEX)
sex.status <- gsub(pattern = "1", replacement = "Male", sex.status)
sex.status <- gsub(pattern = "0", replacement = "Female", sex.status)

Th2.group <- as.character(covariates$group_th)

# Filter the CpG sites whether hyper-/ or hypomethylated
DMR.results <- DMR.results[rownames(heatmap.data), ]
DMR.results$Direction <- "Hypermethylated"
DMR.results[DMR.results$estimate < 0, "Direction"] <- "Hypomethylated"
cpg.direction <- DMR.results$Direction

colColours <- cbind(health.status, smoking.status, sex.status, Th2.group) 

rowColours <- cbind(cpg.direction)
max_abs_value <- max(abs(4))
breaks <- seq(-max_abs_value, max_abs_value, length.out = 101)

# Define color palette for the heatmap
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(101)

# Ensure all annotation vectors match the dimensions of the data
n_cols <- ncol(heatmap.data)
n_rows <- nrow(heatmap.data)

# Ensure annotation vectors match the dimensions of heatmap.data
if (length(health.status) != n_cols) health.status <- rep(NA, n_cols)
if (length(smoking.status) != n_cols) smoking.status <- rep(NA, n_cols)
if (length(sex.status) != n_cols) sex.status <- rep(NA, n_cols)
if (length(Th2.group) != n_cols) Th2.group <- rep(NA, n_cols)
if (length(cpg.direction) != n_rows) cpg.direction <- rep(NA, n_rows)

# Define column and row annotations with names to match heatmap data
annotation_col <- data.frame(
  HealthStatus = health.status,
  SmokingStatus = smoking.status,
  Sex = sex.status,
  Th2Group = Th2.group
)
rownames(annotation_col) <- colnames(heatmap.data)

annotation_row <- data.frame(
  CpG_Direction = cpg.direction
)
rownames(annotation_row) <- rownames(heatmap.data)

# Define annotation colors
annotation_colors <- list(
  HealthStatus = c("Healthy" = "#1f77b4", "Asthma" = "#ff7f0e"),
  SmokingStatus = c("Current smoker" = "#238B45", "Ex-smoker" = "#74C476", "Never-smoker" = "#BAE4B3"),
  Sex = c("Male" = "#4B9CD3", "Female" = "#ff1493"),
  Th2Group = c("healthy" = "#1f77b4", "high" = "#e6550d", "low" = "#ffb74d", "undeterm/cor_bio" = "black"),
  CpG_Direction = c("Hypermethylated" = "red", "Hypomethylated" = "blue")
)

# Replace NA in annotations with a default color if desired
annotation_col[is.na(annotation_col)] <- "grey"
annotation_row[is.na(annotation_row)] <- "grey"

# Save the heatmap to a PDF file
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/"
pdf(file = file.path(base.dir, "Results/006 Differential Methylated Region Analysis DMR/heatmap_sig_DMRs_AsthmaHealthy_FINAL.pdf"),
    width = 5, height = 5)
tiff(file.path("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS",
               "/Results/006 Differential Methylated Region Analysis DMR/heatmap_sig_DMRs_AsthmaHealthy_FINAL.tiff"), width = 3000, height = 3000, res = 300)
# Generate the heatmap with annotations on both axes
pheatmap(
  heatmap.data,                           # Matrix of values to plot
  color = heatmap_colors,                 # Heatmap color palette
  breaks = breaks,                        # Define break points for colors
  cluster_rows = T,                    # Enable row clustering
  cluster_cols = T,                    # Enable column clustering
  annotation_col = annotation_col,        # Column annotations (top bar)
  annotation_row = annotation_row,        # Row annotations (side bar)
  annotation_colors = annotation_colors,  # Custom colors for annotations
  annotation_legend = F,               # Display annotation legend
  annotation_names_row = F,
  annotation_names_col = F,
  show_rownames = FALSE,                  # Hide row names
  show_colnames = FALSE,                  # Hide column names
  scale = "row"                           # Scale across rows
)
dev.off()
