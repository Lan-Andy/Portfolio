# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
# output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/"
# output.dir.plots <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/"

output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/05_pca/"
##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("dplyr")
library("ggplot2")
library("minfi")
library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
dataDirectory <- "/groups/umcg-griac/tmp02/projects/ATLANTIS_2023-396/idats"
# head(list.files(dataDirectory, recursive = TRUE))

## Read samples from sample sheet
targets <- read.metharray.sheet(
  dataDirectory,
  pattern = "EPIC2023-396_ATLANTIS.csv",
  recursive = T,
  verbose = TRUE
)
# targets[1:5,]

##############################################################################
#################### Substitute the missing values into 'NULL' ###############
##############################################################################
targets <- apply(targets, 2, function(x) {
  x[is.na(x)] <- 'NULL'
  x
})
targets <- data.frame(targets)
# targets[1:5, ]

filenames <- targets
# filenames

##############################################################################
#################### Alternative RAMWAS PCA control probes ###################
##############################################################################
# Paragraph 1.6.1.
# https://www.bioconductor.org/packages/release/bioc/vignettes/ramwas/inst/doc/RW5a_matrix.html#:~:text=i.e.%2C%20power).-,1.6.1,Principal%20components%20analysis%20(PCA)%20on%20control%20probes,-We%20extract%20red
# However, they are using all control types, even the background and negative control probes.

for (i in 1:nrow(filenames)) {
  print(i)

  ## Creating a new rgSet, where no samples are deleted.
  ## Dataset contains good and bad quality samples, should not make a big difference
  ## Calculating the principal components per sample, bad quality samples will be
  ##  deleted from the principal component values
  RGset <- read.metharray(filenames$Basename[i], verbose = TRUE)
  RGset@annotation <- c(
    array = 'IlluminaHumanMethylationEPIC',
    annotation = '20a1.hg38'
  )
  # Due to correcting for background, we are not pulling out the probes
  #   that are expected to see background intensity
  #   - Andy
  RGset <- bgcorrect.illumina(RGset)

  ## Checking whether changing the array name will pull out different control probes
  ##    This is not the case, still got all 633 original control probes
  # print("Printing controlType:")
  # controlType = unique(getManifest(RGset)@data$TypeControl$Type)
  # Different controlTypes like STAINING, EXTENSION, BISULFITE CONVERSION I, etc.

  # print("Printing controlSet:")
  # controlSet = getControlAddress(RGset, controlType = controlType)
  # length(controlSet) #633 CpG sites

  # print("Printing probeinfo of controls")
  # control.info <- getProbeInfo(RGset, type = "Control")
  # print(control.info)
  # write.csv(x = control.info, file = "control_info_Lehne.csv")

  #BSC1 control probes
  BSCI.Green.Name = getProbeInfo(RGset, type = "Control")[
    c(14:16),
  ]$ExtendedType
  # print(BSCI.Green.Name)
  BSCI.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(BSCI.Green.Name),
    dimnames = list(BSCI.Green.Name, sampleNames(RGset))
  )
  # print(head(BSCI.Green))
  BSCI.Green[BSCI.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[c(14:16, 19:21), ]$Address,
  ]
  # print(head(BSCI.Green))
  BSCI.Red.Name = getProbeInfo(RGset, type = "Control")[c(17:18), ]$ExtendedType
  # print(BSCI.Red.Name)
  BSCI.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(BSCI.Red.Name),
    dimnames = list(BSCI.Red.Name, sampleNames(RGset))
  )
  BSCI.Red[BSCI.Red.Name, ] <- getRed(RGset)[
    getProbeInfo(RGset, type = "Control")[c(17:18, 22:23), ]$Address,
  ]

  # #BSC2 control probes
  BSCII.Red.Name = getProbeInfo(RGset, type = "Control")[24:27, ]$ExtendedType
  # print(BSCII.Red.Name)
  BSCII.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(BSCII.Red.Name),
    dimnames = list(BSCII.Red.Name, sampleNames(RGset))
  )
  BSCII.Red[BSCII.Red.Name, ] <- getRed(RGset)[
    getProbeInfo(RGset, type = "Control")[24:27, ]$Address,
  ]

  # #STAINING
  stain.Red.Name = getProbeInfo(RGset, type = "Control")[3, ]$ExtendedType
  # print(stain.Red.Name)
  stain.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(stain.Red.Name),
    dimnames = list(stain.Red.Name, sampleNames(RGset))
  )
  stain.Red[stain.Red.Name, ] <- getRed(RGset)[
    getProbeInfo(RGset, type = "Control")[3, ]$Address,
  ]
  stain.Green.Name = getProbeInfo(RGset, type = "Control")[1, ]$ExtendedType
  # print(stain.Green.Name)
  stain.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(stain.Green.Name),
    dimnames = list(stain.Green.Name, sampleNames(RGset))
  )
  stain.Green[stain.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[1, ]$Address,
  ]

  # #EXTENSION
  extensionA.Red.Name = getProbeInfo(RGset, type = "Control")[7, ]$ExtendedType
  # print(extensionA.Red.Name)
  extensionA.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(extensionA.Red.Name),
    dimnames = list(extensionA.Red.Name, sampleNames(RGset))
  )
  extensionA.Red[extensionA.Red.Name, ] <- getRed(RGset)[
    getProbeInfo(RGset, type = "Control")[7, ]$Address,
  ]
  extensionT.Red.Name = getProbeInfo(RGset, type = "Control")[6, ]$ExtendedType
  # print(extensionT.Red.Name)
  extensionT.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(extensionT.Red.Name),
    dimnames = list(extensionT.Red.Name, sampleNames(RGset))
  )
  extensionT.Red[extensionT.Red.Name, ] <- getRed(RGset)[
    getProbeInfo(RGset, type = "Control")[6, ]$Address,
  ]
  extensionC.Green.Name = getProbeInfo(RGset, type = "Control")[
    8,
  ]$ExtendedType
  # print(extensionC.Green.Name)
  extensionC.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(extensionC.Green.Name),
    dimnames = list(extensionC.Green.Name, sampleNames(RGset))
  )
  extensionC.Green[extensionC.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[8, ]$Address,
  ]
  extensionG.Green.Name = getProbeInfo(RGset, type = "Control")[
    5,
  ]$ExtendedType
  # print(extensionG.Green.Name)
  extensionG.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(extensionG.Green.Name),
    dimnames = list(extensionG.Green.Name, sampleNames(RGset))
  )
  extensionG.Green[extensionG.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[5, ]$Address,
  ]

  # #HYBRIDISATION
  hybridH.Green.Name = getProbeInfo(RGset, type = "Control")[11, ]$ExtendedType
  # print(hybridH.Green.Name)
  hybridH.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(hybridH.Green.Name),
    dimnames = list(hybridH.Green.Name, sampleNames(RGset))
  )
  hybridH.Green[hybridH.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[11, ]$Address,
  ]
  hybridM.Green.Name = getProbeInfo(RGset, type = "Control")[10, ]$ExtendedType
  # print(hybridM.Green.Name)
  hybridM.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(hybridM.Green.Name),
    dimnames = list(hybridM.Green.Name, sampleNames(RGset))
  )
  hybridM.Green[hybridM.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[10, ]$Address,
  ]
  hybridL.Green.Name = getProbeInfo(RGset, type = "Control")[9, ]$ExtendedType
  # print(hybridL.Green.Name)
  hybridL.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(hybridL.Green.Name),
    dimnames = list(hybridL.Green.Name, sampleNames(RGset))
  )
  hybridL.Green[hybridL.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[9, ]$Address,
  ]

  # #TARGET REMOVAL
  target.Green.Name = getProbeInfo(RGset, type = "Control")[
    12:13,
  ]$ExtendedType
  # print(target.Green.Name)
  target.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(target.Green.Name),
    dimnames = list(target.Green.Name, sampleNames(RGset))
  )
  target.Green[target.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[12:13, ]$Address,
  ]

  # #Specificity I (PM)
  specI.Green.Name = getProbeInfo(RGset, type = "Control")[28:30, ]$ExtendedType
  # print(specI.Green.Name)
  specI.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(specI.Green.Name),
    dimnames = list(specI.Green.Name, sampleNames(RGset))
  )
  specI.Green[specI.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[28:30, ]$Address,
  ]
  specI.Red.Name = getProbeInfo(RGset, type = "Control")[34:36, ]$ExtendedType
  # print(specI.Red.Name)
  specI.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(specI.Red.Name),
    dimnames = list(specI.Red.Name, sampleNames(RGset))
  )
  specI.Red[specI.Red.Name, ] <- getRed(RGset)[
    getProbeInfo(RGset, type = "Control")[34:36, ]$Address,
  ]

  # #Specificity II
  specII.Red.Name = getProbeInfo(RGset, type = "Control")[40:42, ]$ExtendedType
  # print(specII.Red.Name)
  specII.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(specII.Red.Name),
    dimnames = list(specII.Red.Name, sampleNames(RGset))
  )
  specII.Red[specII.Red.Name, ] <- getRed(RGset)[
    getProbeInfo(RGset, type = "Control")[40:42, ]$Address,
  ]

  # #NON POLYMORPHIC
  np.Red.Name = getProbeInfo(RGset, type = "Control")[43:44, ]$ExtendedType
  # print(np.Red.Name)
  np.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(np.Red.Name),
    dimnames = list(np.Red.Name, sampleNames(RGset))
  )
  np.Red[np.Red.Name, ] <- getRed(RGset)[
    getProbeInfo(RGset, type = "Control")[43:44, ]$Address,
  ]
  np.Green.Name = getProbeInfo(RGset, type = "Control")[45:46, ]$ExtendedType
  # print(np.Green.Name)
  np.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(np.Green.Name),
    dimnames = list(np.Green.Name, sampleNames(RGset))
  )
  np.Green[np.Green.Name, ] <- getGreen(RGset)[
    getProbeInfo(RGset, type = "Control")[45:46, ]$Address,
  ]

  # #Normalisation
  control = getProbeInfo(RGset, type = "Control")
  normC.Green.Name = control[control[, 2] == 'NORM_C', 4]
  # print(normC.Green.Name)
  normC.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(normC.Green.Name),
    dimnames = list(normC.Green.Name, sampleNames(RGset))
  )
  normC.Green[normC.Green.Name, ] <- getGreen(RGset)[
    control[control[, 2] == 'NORM_C', 1],
  ]
  normG.Green.Name = control[control[, 2] == 'NORM_G', 4]
  # print(normG.Green.Name)
  normG.Green <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(normG.Green.Name),
    dimnames = list(normG.Green.Name, sampleNames(RGset))
  )
  normG.Green[normG.Green.Name, ] <- getGreen(RGset)[
    control[control[, 2] == 'NORM_G', 1],
  ]
  normA.Red.Name = control[control[, 2] == 'NORM_A', 4]
  # print(normA.Red.Name)
  normA.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(normA.Red.Name),
    dimnames = list(normA.Red.Name, sampleNames(RGset))
  )
  normA.Red[normA.Red.Name, ] <- getRed(RGset)[
    control[control[, 2] == 'NORM_A', 1],
  ]
  normT.Red.Name = control[control[, 2] == 'NORM_T', 4]
  # print(normT.Red.Name)
  normT.Red <- matrix(
    NA_real_,
    ncol = ncol(RGset),
    nrow = length(normT.Red.Name),
    dimnames = list(normT.Red.Name, sampleNames(RGset))
  )
  normT.Red[normT.Red.Name, ] <- getRed(RGset)[
    control[control[, 2] == 'NORM_T', 1],
  ]

  #combine ctrl probe intensities
  ctrl = rbind(
    as.matrix(BSCI.Green),
    as.matrix(BSCI.Red),
    as.matrix(BSCII.Red),
    as.matrix(stain.Red),
    as.matrix(stain.Green),
    as.matrix(extensionA.Red),
    as.matrix(extensionT.Red),
    as.matrix(extensionC.Green),
    as.matrix(extensionG.Green),
    as.matrix(hybridH.Green),
    as.matrix(hybridM.Green),
    as.matrix(hybridL.Green),
    as.matrix(target.Green),
    as.matrix(specI.Green),
    as.matrix(specI.Red),
    as.matrix(specII.Red),
    as.matrix(np.Red[1, ]),
    as.matrix(np.Red[2, ]),
    as.matrix(np.Green[1, ]),
    as.matrix(np.Green[2, ]),
    as.matrix(normC.Green),
    as.matrix(normG.Green),
    as.matrix(normA.Red),
    as.matrix(normT.Red)
  )

  #add data for the new samples
  if (exists("ctrl.all")) {
    ctrl.all <- rbind(ctrl.all, t(ctrl))
  } else {
    ctrl.all <- t(ctrl)
  }
}

print("Printing all control probes:")
ctrl.all

write.csv(
  x = ctrl.all,
  file = file.path(output.dir, "ctrlProbes/all_control_info_Lehne.csv")
)
ctrl.all <- read.csv(
  file = file.path(output.dir, "ctrlProbes/all_control_info_Lehne.csv"),
  header = TRUE,
  row.names = 1
)
# head(ctrl.all)
# str(ctrl.all)
# dim(ctrl.all)
# class(ctrl.all)
################################################################################
# Lehne's way of calculating the principal components,
#   within GRIAC we are using a different method, see below
# # PCA of control-probe intensities
# pca <- prcomp(na.omit(ctrl.all))
# ctrlprobes.scores = pca$x
# colnames(ctrlprobes.scores) = paste(colnames(ctrlprobes.scores), '_cp', sep='')
# write.csv(x = ctrlprobes.scores, file = "control_probes_PCs_Lehne.csv")
################################################################################

# ##############################################################################
# ##################### Principal component analysis ###########################
# ##############################################################################
"
PCA is a standard function in R (including many options)
We will perform the PCA on log transformed data.
But note that the log function in R, by default, will take natural log of the value
In proteomics and transcriptomics data is always log2 transformed

Important is to know that the pca function in R needs the genes in columns and experiments in row
We can simply transform the table using the t function in R
"

do.center = TRUE
do.scale = TRUE

# pca <- prcomp(pca.data, center=TRUE, scale=TRUE)
methylation.pca <- ctrl.all %>%
  stats::prcomp(
    center = do.center,
    scale. = do.scale
  )

save(methylation.pca, file = file.path(output.dir, "ctrlProbes/PCA.RData")) # Variable name==methylation.pca
# Write summary of PCAs to files
# "The function 'summary' will show the contribution of the Principle Components to the variance"
pcs.summary <- summary(methylation.pca)

# "The content of the data structure methylation.pca"
head(methylation.pca$rotation)
head(methylation.pca$center)
head(methylation.pca$x)

pcs.summary$importance %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("PC.name") %>%
  readr::write_csv(
    file.path(output.dir, "ctrlProbes/importance.csv")
  )

pcs.summary$x %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("CpG.sites") %>%
  readr::write_csv(
    file.path(output.dir, "ctrlProbes/values.csv")
  )

pcs.summary$rotation %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column("sample.id") %>%
  readr::write_csv(
    file.path(output.dir, "ctrlProbes/rotation.csv")
  )

data.frame(
  rownames = names(pcs.summary$center),
  center = pcs.summary$center,
  scale = pcs.summary$scale
) %>%
  readr::write_csv(
    file.path(output.dir, "ctrlProbes/rest.csv")
  )


##############################################################################
############### Calculating principal components for visualisation ###########
##############################################################################
# Calculating the eigen values of these Principal Components
#   together with the cumulative variance of principal component
eigen_values <- methylation.pca[["sdev"]]^2
eigen_values <- as.data.frame(cbind(
  "Eigen" = eigen_values,
  "Variance_percent" = eigen_values * 100 / sum(eigen_values),
  "Cumulative_variance_percent" = cumsum(eigen_values * 100 / sum(eigen_values))
))
## Rename eigen values towads PC1 PC2 etc.
rownames(eigen_values) <- paste0(
  "PC",
  seq(1, dim(eigen_values)[1], 1),
  sep = ""
)
eigen_values <- cbind(eigen_values, "PC" = seq(1, dim(eigen_values)[1], 1))


# Screeplot
screeplot <-
  ggplot(
    data = eigen_values[eigen_values$PC %in% c(1:15), ],
    aes(x = PC, y = Variance_percent)
  ) +
  scale_x_continuous(breaks = 1:15) +
  scale_y_continuous(breaks = c(seq(from = 0, to = 100, by = 20))) +
  xlab("Principal Component") +
  ylab("Percentage of explained variances") +
  geom_bar(stat = "identity", fill = "lightblue") +
  geom_hline(yintercept = 86.39, linetype = "dashed", colour = "red") +
  annotate(
    "text",
    x = 6,
    y = 73,
    label = "Explaining 86.39% of total variance"
  ) +
  geom_point(colour = "grey50") +
  geom_line(colour = "grey50") +
  geom_point(
    data = eigen_values[eigen_values$PC %in% c(1:15), ],
    aes(x = PC, y = Cumulative_variance_percent)
  ) +
  geom_line(
    data = eigen_values[eigen_values$PC %in% c(1:15), ],
    aes(x = PC, y = Cumulative_variance_percent)
  ) +
  ggtitle("Screeplot") +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

png(
  file = file.path(output.dir, "ctrlProbes/screePlot.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300,
  pointsize = 12,
  bg = "white"
)
print(screeplot)
dev.off()

# "If we can plot the contribution of the Princple Components to the variance
#   we can observe that the first 2 PC's contribute the most to the Variance"
# graph1 <- plot(methylation.pca, type = "l")
#
# ggsave(plot = graph1, filename = file.path(output.dir, "ctrlProbes/pca_lines.tiff"),
#        dpi = 300)

# "We can plot the first 2 Principle Components for all experiments using the plot function of R"
graph2 <- plot(methylation.pca$x, type = "p")

ggsave(
  plot = graph2,
  filename = "pca_points.tiff",
  path = output.dir,
  dpi = 300
)

# "Grab the principal component data for further visualisation"
my_pca_data <- as.data.frame(methylation.pca$x)

# "Removing the X from columnnames if they have been imported differently
rownames(methylation.pca$x) <- gsub(
  pattern = "X",
  replacement = "",
  x = rownames(methylation.pca$x)
)

# "Load in the linking table to find the samples we need
linkingtable <- read.csv(
  "/groups/umcg-griac/tmp02/projects/TatianaKarp/ATL_methyaltion/methID_rnaID_clinID.csv",
  header = TRUE,
  row.names = 1
)
my_pca_data$meth_file_id <- rownames(methylation.pca$x)[
  rownames(methylation.pca$x) %in% linkingtable$meth_file_id
]

# Linking table containing the principal components of control probes
# This file will be used for the differential methylation analysis
readr::write_csv(
  x = my_pca_data,
  file = file.path(output.dir, "ctrlProbes/my_pca_data.csv")
)


# ##############################################################################
# ############### Visualisation principal components ###########################
# ##############################################################################

# Merge the principal component together with pheno data in linkingtable
my_pca_data <- merge(my_pca_data, linkingtable, by = "meth_file_id")

PCAplot <- ggplot(
  my_pca_data,
  aes(x = PC1, y = PC2, label = factor(Array), color = factor(Array))
) +
  geom_point(pch = 18, cex = 4)
PCAplot <- PCAplot +
  geom_text(
    hjust = -0.4,
    size = 3,
    check_overlap = TRUE,
    aes(colour = factor(Array))
  )
PCAplot <- PCAplot + ggtitle("PCA plot") + theme_bw()

png(
  file = file.path(output.dir, "ctrlProbes/PCA_plot_Array.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300,
  pointsize = 12,
  bg = "white"
)
print(PCAplot)
dev.off()

PCAplot <- ggplot(
  my_pca_data,
  aes(x = PC1, y = PC2, label = factor(Slide), color = factor(Slide))
) +
  geom_point(pch = 18, cex = 4)
PCAplot <- PCAplot +
  geom_text(
    hjust = -0.4,
    size = 3,
    check_overlap = TRUE,
    aes(colour = factor(Slide))
  )
PCAplot <- PCAplot + ggtitle("PCA plot") + theme_bw()

png(
  file = file.path(output.dir, "ctrlProbes/PCA_plot_Slide.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300,
  pointsize = 12,
  bg = "white"
)
print(PCAplot)
dev.off()

PCAplot <- ggplot(
  my_pca_data,
  aes(x = PC1, y = PC2, label = factor(GENDER), color = factor(GENDER))
) +
  geom_point(pch = 18, cex = 4)
PCAplot <- PCAplot +
  geom_text(
    hjust = -0.4,
    size = 3,
    check_overlap = TRUE,
    aes(colour = factor(GENDER))
  )
PCAplot <- PCAplot + ggtitle("PCA plot") + theme_bw()

png(
  file = file.path(output.dir, "ctrlProbes/PCA_plot_GENDER.png"),
  width = 10,
  height = 10,
  units = "in",
  res = 300,
  pointsize = 12,
  bg = "white"
)
print(PCAplot)
dev.off()
