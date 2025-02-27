# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
output.dir.plots <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results"


##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("tidyverse")
library("ggplot2")
library("ggrepel")
library("PupillometryR")
library("scales")
library("dplyr")

library("heatmap3")
library("scales")

library("pheatmap")
library("scales")
##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################

meth.data <- read.csv(file = file.path(base.dir, "Data/005 eQTM analysis/methylationData_eQTM_samples_asthma_vs_control.csv"),
                      header = TRUE,
                      row.names = 1)
DMA.result <- read.csv(file = file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3/Tt2_significant_CpG_asthma_vs_control_BH.csv"),
                       header = TRUE,
                       row.names = 1)

meth.data <- read.csv(file = file.path(base.dir, "Data/005 eQTM analysis/methylationData_eQTM_AllSamples_Th2HighvsTh2Low.csv"),
                      header = TRUE,
                      row.names = 1)
DMA.result <- read.csv(file = file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3 - Th2 groups/DMA_Tt2_significant_CpG_Th2HighvsTh2Low_BH.csv"),
                       header = TRUE,
                       row.names = 1)
covariates <- read.csv(file = file.path(base.dir, "Data/005 eQTM analysis/covariates_eQTM_samples_original.csv"),
                       header = TRUE,
                       row.names = 1)
################################################################################
# Need to get new Mvalues with the right probe ID's
meth.data <- read.csv(file = file.path(base.dir,"Data/002 Differential Methylation Analysis DMA/Correct for PC1-3/sig_AsthmaVsHealthy_ATLANTIS_Mvalues_dasen.csv"),
                      header = TRUE,
                      row.names = 1)

# pheno.data <- read.csv(file = file.path(base.dir, "Data/100 Characteristics tables/pheno_methylation_samples.csv"),
#                        header = TRUE,
#                        row.names = 1)

linkingtable <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/methID_rnaID_clinID.csv")

pheno.data <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/atlantis_patient_data.csv",
                       header = TRUE, row.names = 1, na.strings = "")

PC.all.cp <- read.csv("C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/ATLANTIS QC/Data/PCA technical probes/my_pca_data.csv")

DMA.result <- read.csv(file = file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3/Tt2_significant_CpG_asthma_vs_control_BH.csv"),
                       header = TRUE,
                       row.names = 1)
# DMA.result <- read.csv(file = file.path(base.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3 - Th2 groups/DMA_Tt2_significant_CpG_Th2HighvsTh2Low_BH.csv"),
#                        header = TRUE,
#                        row.names = 1)
##############################################################################
################ Change column names for easier interpretation ###############
##############################################################################
colnames(meth.data) <- gsub(pattern = "X", replacement = "", x = colnames(meth.data))

# colnames(m.values) <- gsub(pattern = "X", replacement = "", x = colnames(m.values))

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

# covariates$group_th <- pd.methyl.sub[which(pd.methyl.sub$PT %in% covariates$PT), "group_th"]
# covariates$group_th <- gsub(pattern = "cor_bio|undeterm", replacement = "undeterm/cor_bio" , covariates$group_th)
# covariates <- covariates[grep(pattern = "healthy|high|low", covariates$group_th), ]



# covariates <- covariates[order(covariates$group_th), ]
# meth.data <- meth.data[, covariates$meth_file_id]
# m.values <- m.values[, covariates$meth_file_id]


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


##############################################################################
################ Heatmap ###############
##############################################################################
## Create the heatmap
## Colv = dendrogram computed and reordered, we disabled here as we want to
##  to group the column with their respective experimental group (Asthma/Healthy)
## It is possible to cluster the CpG sites based on their methylation; Rowv
## labRow and labCol are the labels that belong respective to their row/column,
##  you should the give rownames()/colnames() of the table



heatmap.data <- as.matrix(meth.data[, colnames(meth.data)])
# heatmap.data <- heatmap.data[1:50, ]
# heatmap.data <- as.matrix(meth.data)


health.status <- as.character(covariates$ASTHEA)
health.status <- gsub(pattern = "1", replacement = "Asthma", health.status)
health.status <- gsub(pattern = "0", replacement = "Healthy", health.status)
# health.status[health.status=="Healthy"] <- "#1f77b4"
# health.status[health.status=="Asthma"] <- "#ff7f0e"

smoking.status <- as.character(covariates$SMOKE)
smoking.status <- gsub(pattern = "2", replacement = "Current smoker", smoking.status)
smoking.status <- gsub(pattern = "1", replacement = "Ex-smoker", smoking.status)
smoking.status <- gsub(pattern = "0", replacement = "Never-smoker", smoking.status)
# smoking.status[smoking.status=="Current smoker"] <- "#238B45"
# smoking.status[smoking.status=="Ex-smoker"] <- "#74C476"
# smoking.status[smoking.status=="Never-smoker"] <- "#BAE4B3"

sex.status <- as.character(covariates$SEX)
sex.status <- gsub(pattern = "1", replacement = "Male", sex.status)
sex.status <- gsub(pattern = "0", replacement = "Female", sex.status)
# sex.status[sex.status=="Male"] <- "#4B9CD3"
# sex.status[sex.status=="Female"] <- "#ff1493"

Th2.group <- as.character(covariates$group_th)
# Th2.group[Th2.group=="healthy"] <- "#1f77b4"
# Th2.group[Th2.group=="high"] <- "#e6550d"
# Th2.group[Th2.group=="low"] <- "#ffb74d"
# Th2.group[Th2.group=="undeterm/cor_bio"] <- "black"

# Filter the CpG sites whether hyper-/ or hypomethylated
DMA.result <- DMA.result[rownames(heatmap.data), ]
DMA.result$Direction <- "Hypermethylated"
DMA.result[DMA.result$logFC < 0, "Direction"] <- "Hypomethylated"
cpg.direction <- DMA.result$Direction
# cpg.direction[cpg.direction=="Hypermethylated"] <- "blue"
# cpg.direction[cpg.direction=="Hypomethylated"] <- "red"

## Incorporating kmeans clustering
# set.seed(22)
# num_clusters <- 10
# cluster_colors <- setNames(RColorBrewer::brewer.pal(n = num_clusters, name = "Set3"), as.character(1:num_clusters))
# obj <- pheatmap(t(heatmap.data), kmeans_k = num_clusters, cluster_cols = FALSE, display_numbers = TRUE)
# cluster.members <- as.data.frame(as.factor(obj$kmeans$cluster))


colColours <- cbind(health.status, smoking.status, sex.status, Th2.group) 
colColours <- cbind(health.status, smoking.status, sex.status) 

rowColours <- cbind(cpg.direction)
max_abs_value <- max(abs(4))
breaks <- seq(-max_abs_value, max_abs_value, length.out = 101)


# pdf(file = file.path(base.dir, "/Results/002 Differential Methylation Analysis DMA/heatmap_sig_CpGs_AsthmaHealthy.pdf"),
#     width = 5, height = 5)
# heatmap3(heatmap.data,
#          ColSideColors = colColours,
#          RowSideColors = rowColours,
#          scale = "row",
#          col = colorRampPalette(c("blue", "#FFFFFF", "#e00000"))(100),
#          breaks = breaks,
#          Rowv=TRUE, Colv=TRUE,
#          labRow = NA, labCol = NA)
# dev.off()

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
  Sex = sex.status
  # Th2Group = Th2.group
  # cluster = cluster.members
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
  # Cluster = cluster_colors
)

# Replace NA in annotations with a default color if desired
annotation_col[is.na(annotation_col)] <- "grey"
annotation_row[is.na(annotation_row)] <- "grey"

# Save the heatmap to a PDF file
# pdf(file = file.path(base.dir, "/Results/002 Differential Methylation Analysis DMA/heatmap_sig_CpGs_AsthmaHealthy_FINAL.pdf"),
# width = 5, height = 5)
# set.seed(22)
# Generate the heatmap with annotations on both axis

# Save as tiff file is faster than pdf
tiff(file.path(base.dir, "/Results/002 Differential Methylation Analysis DMA/heatmap_sig_CpGs_AsthmaHealthy_FINAL.tiff"), width = 3000, height = 3000, res = 300)

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


##############################################################################
################ Volcano plots | Asthma vs Control ###########################
##############################################################################
base.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/002 Differential Methylation Analysis DMA/"

tT <- read.csv(file = file.path(base.dir, "Correct for PC1-3/Tt_CpG_asthma_vs_control_BH.csv"),
               header = TRUE, row.names = 1)
tT2 <- read.csv(file = file.path(base.dir, "Correct for PC1-3/Tt_CpG_asthma_vs_control_bonferroni.csv"),
               header = TRUE, row.names = 1)


tT <- read.csv(file = file.path(base.dir, "Correct for PC1-3 - Th2 groups/DMA_Tt_CpG_Th2HighVsHealthy_BH.csv"),
               header = TRUE, row.names = 1)
tT2 <- read.csv(file = file.path(base.dir, "Correct for PC1-3 - Th2 groups/DMA_Tt_CpG_Th2HighVsHealthy_bonferroni.csv"),
                header = TRUE, row.names = 1)

tT <- read.csv(file = file.path(base.dir, "Correct for PC1-3 - Th2 groups/DMA_Tt_CpG_Th2LowVsHealthy_BH.csv"),
               header = TRUE, row.names = 1)
tT2 <- read.csv(file = file.path(base.dir, "Correct for PC1-3 - Th2 groups/DMA_Tt_CpG_Th2LowVsHealthy_bonferroni.csv"),
                header = TRUE, row.names = 1)

tT <- read.csv(file = file.path(base.dir, "Correct for PC1-3 - Th2 groups/DMA_Tt_CpG_Th2HighvsTh2Low_BH.csv"),
               header = TRUE, row.names = 1)
tT2 <- read.csv(file = file.path(base.dir, "Correct for PC1-3 - Th2 groups/DMA_Tt_CpG_Th2HighvsTh2Low_bonferroni.csv"),
                header = TRUE, row.names = 1)


tT$Legend <- ifelse(
  tT$adj.P.Val < 0.05,
  ifelse(
    tT$logFC < 0,
    "Hypomethylated",
    "Hypermethylated"
  ),
  "Not Significant"
)
# tT2$Legend <- ifelse(
#   tT2$adj.P.Val < 0.05,
#   ifelse(
#     tT2$logFC < 0,
#     "Hypomethylated",
#     "Hypermethylated"
#   ),
#   "Not Significant"
# )

annotation_data <- tT %>%
  dplyr::arrange(P.Value) %>%
  head(, n = 20) %>%
  tibble::rownames_to_column("cpg")

benjamini.intercept <- tT %>%
  dplyr::filter(
    as.character(Legend) == "Not Significant"
  ) %>%
  dplyr::pull(P.Value) %>%
  min() %>%
  log10() * -1

# bonferroni.intercept <- tT2 %>%
#   dplyr::filter(
#     as.character(Legend) == "Not Significant"
#   ) %>%
#   dplyr::pull(P.Value) %>%
#   min() %>%
#   log10() * -1

DMA_volcano <- ggplot2::ggplot(
  data = tT,
  mapping = ggplot2::aes(
    x = logFC,
    y = -log10(P.Value)
  )
) +
  theme_bw() +
  ggplot2::geom_point(
    mapping = aes(
      color = Legend
    )
  ) +
  ggplot2::ylab("-log<sub>10</sub>(P-value)") +
  ggplot2::xlab("log<sub>2</sub>(Fold Change)") +
  ggplot2::scale_color_manual(
    values = c(
      'Hypomethylated' = "blue",
      'Hypermethylated' = "red",
      'Not Significant' = "grey"
    )
  ) +
  ggplot2::geom_hline(
    yintercept = benjamini.intercept,
    colour = "black",
    linetype = "dashed"
  ) +
  # ggplot2::geom_hline(
  #   yintercept = bonferroni.intercept,
  #   colour = "forestgreen",
  #   linetype = "dashed"
  # ) +
  ggrepel::geom_label_repel(
    data = annotation_data,
    mapping = ggplot2::aes(
      label = cpg
    ),
    size = 4,
    min.segment.length = 0,
    max.overlaps = 20
  ) +
  ggplot2::guides(
    color = guide_legend(
      override.aes = list(size = 3)
    )
  ) +
  # ggprism::theme_prism() +
  # ggplot2::xlim(x_limits) + # Set x-axis limits
  # ggplot2::ylim(y_limits) + # Set y-axis limits

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

ggsave(DMA_volcano, #path = base.dir, 
       file="C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/002 Differential Methylation Analysis DMA/DMA_volcano_Th2HighVsHealthy.tiff",
       device = "tiff", width = 18, height = 18, units = "cm")

# If you want to have the same axes through all volcano plot,
#   decide which limits you want to use 
plot_build <- ggplot2::ggplot_build(DMA_volcano)
x_limits <- plot_build$layout$panel_params[[1]]$x.range
y_limits <- plot_build$layout$panel_params[[1]]$y.range


##############################################################################
################ If boxplots are needed from certain CpG sites ###############
##############################################################################
# The code below has been incorpated into the eQTM figure, as the marginal plot
#   (raincloud plot) are on the side and on the top of the eQTM scatter plot
# Might want to use this part when just want to look at certain CpG sites.

# ##############################################################################
# ################ Filter the top 20 CpG sites from DMA analysis ###############
# ##############################################################################
# m.values.20 <- as.data.frame(t(m.values[1:20, ]))
# m.values.20 <- m.values.20[covariates$meth_file_id, ]
# m.values.20$asthma <- covariates[which(covariates$meth_file_id %in% rownames(m.values.20)), "ASTHEA"]
# 
# m.values.20$asthma <- gsub(pattern = "0", replacement = "Healthy", x = m.values.20$asthma)
# m.values.20$asthma <- gsub(pattern = "1", replacement = "Asthma", x = m.values.20$asthma)
# 
# 
# m.values.20$smoke <- covariates[which(covariates$meth_file_id %in% rownames(m.values.20)), "SMOKE"]
# m.values.20$smoke <- gsub(pattern = "0", replacement = "Never-smoker", x = m.values.20$smoke)
# m.values.20$smoke <- gsub(pattern = "1", replacement = "Ex-smoker", x = m.values.20$smoke)
# m.values.20$smoke <- gsub(pattern = "2", replacement = "Current smoker", x = m.values.20$smoke)
# 
# m.values.20$sex <- covariates[which(covariates$meth_file_id %in% rownames(m.values.20)), "SEX"]
# m.values.20$sex <- gsub(pattern = "0", replacement = "Female", x = m.values.20$sex)
# m.values.20$sex <- gsub(pattern = "1", replacement = "Male", x = m.values.20$sex)
# 
# m.values.20$age <- covariates[which(covariates$meth_file_id %in% rownames(m.values.20)), "AGE"]
# 
# head(rownames(DMA.result), n=20)
# meth.data.20 <- meth.data[head(rownames(DMA.result), n=20), ]
# 
# meth.data.20 <- as.data.frame(t(meth.data.20))
# meth.data.20 <- meth.data.20[covariates$meth_file_id, ]
# meth.data.20$asthma <- covariates[which(covariates$meth_file_id %in% rownames(meth.data.20)), "ASTHEA"]
# 
# meth.data.20$asthma <- gsub(pattern = "0", replacement = "Healthy", x = meth.data.20$asthma)
# meth.data.20$asthma <- gsub(pattern = "1", replacement = "Asthma", x = meth.data.20$asthma)
# 
# 
# meth.data.20$smoke <- covariates[which(covariates$meth_file_id %in% rownames(meth.data.20)), "SMOKE"]
# meth.data.20$smoke <- gsub(pattern = "0", replacement = "Never-smoker", x = meth.data.20$smoke)
# meth.data.20$smoke <- gsub(pattern = "1", replacement = "Ex-smoker", x = meth.data.20$smoke)
# meth.data.20$smoke <- gsub(pattern = "2", replacement = "Current smoker", x = meth.data.20$smoke)
# 
# meth.data.20$sex <- covariates[which(covariates$meth_file_id %in% rownames(meth.data.20)), "SEX"]
# meth.data.20$sex <- gsub(pattern = "0", replacement = "Female", x = meth.data.20$sex)
# meth.data.20$sex <- gsub(pattern = "1", replacement = "Male", x = meth.data.20$sex)
# 
# meth.data.20$age <- covariates[which(covariates$meth_file_id %in% rownames(meth.data.20)), "AGE"]
# 
# ##############################################################################
# ################ Create the raincloud plot from the top 20 CpG sites #########
# ##############################################################################
# 
# # https://jtr13.github.io/cc21fall2/raincloud-plot-101-density-plot-or-boxplotwhy-not-do-both.html
# # Raincloud plot
# # Boxplot and density plot next to each other
# pdf(file = file.path(base.dir, "/Results/002 Differential Methylation Analysis DMA/Top20CpG.pdf"),
#     width = 5, height = 5)
# for (i in 1:(ncol(meth.data.20)-4)) {
# plot2 <- ggplot(meth.data.20) +
#   aes(x = asthma,
#       y = meth.data.20[, i],
#       fill = asthma) +
#   geom_flat_violin(position = position_nudge(x = .2),
#                    alpha = .4) +
#   geom_point(aes(color = asthma),
#              position = position_jitter(w = .075),
#              size = 1,
#              alpha = 0.4,
#              show.legend = F) +
#   geom_boxplot(width = .25,
#                outlier.shape = NA,
#                alpha = 0.5) +
#   scale_fill_manual(values = c("#ff7f0e","#1f77b4")) +
#   scale_color_manual(values = c("#ff7f0e","#1f77b4"))  +
#   # expand_limits(x = 3, y = 10) +
#   scale_y_continuous(#limits = c(0,1),
#                      breaks = pretty_breaks(n = 5)
#                      ) +
#   labs(x = "",
#        y = colnames(meth.data.20[i]),
#        fill = "Status ",
#        title = "") +
#   guides(fill = guide_legend(nrow=1,
#                              byrow=TRUE))+
#   theme_bw() +
#   theme(legend.position = "bottom",
#         legend.margin = margin(-5, 0, 0, 0))
# 
#   print(plot2)
# }
# dev.off()
