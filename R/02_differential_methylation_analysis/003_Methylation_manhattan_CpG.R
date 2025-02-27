# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
# output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
# output.dir.plots <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/"


base.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/04_normalisation/"
output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/06_DMA/"

# tT2 <- read.csv(file = file.path(output.dir, "DMA_Tt2_significant_CpG_Th2HighVsHealthy_BH.csv"),
#                 header = TRUE, row.names = 1)
# 
# M.values <- read.csv(file = file.path(base.dir, "ATLANTIS_Mvalues_dasen.csv"),
#                      header = TRUE, row.names = 1)
# M.values <- M.values[rownames(tT2), ]
# write.csv(x = M.values, file = file.path(output.dir, "sig_Th2HighVsHealthy_ATLANTIS_Mvalues_dasen.csv"))
# 
# # B.values <- read.csv(file = file.path(base.dir, "ATLANTIS_betas_dasen.csv"),
# #                      header = TRUE, row.names = 1)
# # B.values <- B.values[rownames(tT2), ]
# # write.csv(x = B.values, file = file.path(output.dir, "sig_AsthmaVsHealthy_ATLANTIS_betas_dasen.csv"))
# 
# quit()
#############################################################################
################## Load used libraries within script #########################
##############################################################################
library("dplyr")
library("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
# library("biomaRt")
# library("org.Hs.eg.db")
library("ggplot2")
library("ggrepel")

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
input.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/06_DMA"
output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/07_anno"

# This gives extensive annotation information about the CpG sites
# annotation <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# tT2 <- read.csv(file = file.path(input.dir, "Tt2_significant_CpG_asthma_vs_control.csv"),
#                header = TRUE, row.names = 1)

##############################################################################
#################### Extract annotation data from CpG site ###################
##############################################################################
# annotation.sub <- annotation[annotation$Name %in% rownames(tT2), c("chr", "pos", "strand",
#                                                "UCSC_RefGene_Group", "UCSC_RefGene_Name",
#                                                "GencodeV41_Group", "GencodeV41_Name")]
# write.csv(annotation.sub, file = file.path(output.dir, "CpG_annotation.csv"))
# rm(annotation.sub)

# probe.anno <- Locations[rownames(tT), ]
# head(probe.anno)

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
# Check whether the CpG site is within a gene 
# This would direct influence the proteinformation
# Create dataframe of ensembl information from Homo sapiens

# mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")
# 
# # Pull out the gene symbol, chromosome, start and end position of the gene.
# gene.anno <- getBM(attributes=c("hgnc_symbol","chromosome_name","start_position","end_position",
#                                 "strand", "band"),
#      filters=c("chromosome_name","biotype"),
#      values=list(chromosome_name=c(1:22,"X","Y"),biotype="protein_coding"),
#      mart=mart) # ensembl human genes
# 
# library(GenomicRanges)
# probesRanges = GenomicRanges::GRanges(seqnames = probe.anno$chr,
#                        ranges = IRanges(start = probe.anno$pos,
#                                       end = probe.anno$pos+1)
#                        )
# probesRanges
# genesRanges = GenomicRanges::GRanges(seqnames = gene.anno$chromosome_name,
#                       ranges = IRanges(start = gene.anno$start_position,
#                                      end = gene.anno$end_position)
#                       )
# genesRanges
# IRanges::findOverlaps(genesRanges, probesRanges)
# 
# rm(probesRanges, genesRanges, gene.anno)

##############################################################################
#################### Manhattan plot asthma vs control ########################
##############################################################################
input.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/002 Differential Methylation Analysis DMA/Correct for PC1-3"
output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Results/003 Annotation"
# Changing this and the standard error file, depending on the hypothesis/input files
# Asthma vs Healthy   = 06_DMA/Tt_CpG_asthma_vs_control_BH.csv
# Th2-high vs Healthy = 06_DMA/DMA_Tt_CpG_Th2HighVsHealthy_BH.csv
# Th2-low vs Healthy  = 06_DMA/DMA_Tt_CpG_Th2LowVsHealthy_BH.csv
# Th2-high vs Th2-low = 06_DMA/DMA_Tt_CpG_Th2HighvsTh2Low_BH.csv

# tT <- read.csv(file = file.path(output.dir, "Data/002 Differential Methylation Analysis DMA/Correct for PC1-3 - Th2 groups/DMA_Tt_CpG_Th2LowVsHealthy_BH.csv"),
#                header = TRUE, row.names = 1)



# tT <- read.csv(file = file.path(input.dir, "Tt_CpG_asthma_vs_control_BH.csv"),
#                header = TRUE, row.names = 1)
# 
# selection <- which(tT$adj.P.Val<0.05)
# tT2 <- tT[selection,]

output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/08_DMR"
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}


data <- loadRData(fileName = file.path(output.dir, "initialData_ATLANTIS_betas.RData"))

# QQ plot and Manhattan plot
# data <- tT
# car::qqPlot(data[, "adj.P.Val"])
# 
annotation <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

annotation$Name <- gsub(pattern = "_.*$", replacement = "", annotation$Name)
annotation <- annotation[!duplicated(annotation$Name),]
rownames(annotation) <- annotation$Name

sel <- rownames(data)
index_anno <- annotation@rownames %in% sel
anno_cor <- subset(x = annotation, index_anno==TRUE)
dim(anno_cor)
# Saving the location of the DMA analysis
location.covariates <- c("chr","pos","strand","UCSC_RefGene_Group","UCSC_RefGene_Name","GencodeV41_Group","GencodeV41_Name")
# write.csv(anno_cor[, location.covariates], file = file.path("/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation/07_anno/CpGLocation_854k.csv"))

# write.csv(anno_cor[, location.covariates], file = file.path(output.dir, "Data/005 eQTM analysis/CpGLocation_854k.csv"))
quit()
# Saving the location of only the significant CpG sites
location.covariates <- c("Name", "chr","pos")
eQTM.location <- anno_cor[rownames(tT2), location.covariates]
colnames(eQTM.location) <- c("cpg", "chr","pos")
# rownames(eQTM.location) <- 1:nrow(eQTM.location)
# write.csv(eQTM.location, file.path(output.dir, "Data/005 eQTM analysis/CpGLocation_asthma_vs_control.csv"))
# write.csv(eQTM.location, file.path(output.dir, "CpGLocation_Th2LowVsHealthy.csv"))

print("Printing data")
data <- data[order(rownames(data)),]
print("Printing anna_cor")
anno_cor <- anno_cor[order(anno_cor@rownames), location.covariates]

print("Printing Man_data")
Man_data <- cbind(data[,1:5], anno_cor)


# Man_data$check <- ifelse(Man_data$x==Man_data$Name, 1,2)
#
SNP = as.data.frame(Man_data$Name)
CHR = as.data.frame(Man_data$chr)
BP= as.numeric(Man_data$pos)
P = Man_data$P.Value
#B = Man_data$Meta_summary

print("Creating Results_data")
Results_data <- as.data.frame(cbind(SNP, CHR, BP, P))
colnames(Results_data) <- c("SNP", "CHR", "BP", "P")
Results_data$CHR <- as.numeric(gsub(pattern = "chr", replacement = "", x = Results_data$CHR))


###########################################################################################################################
# Adapted from qqman R package, adjusted by Andy Lan
# MyManhattan <- function (x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", col = c("gray10", 
#                                                                                   "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
#                          genomewideline = -log10(5e-08), highlight = NULL, logp = TRUE, 
#                          annotatePval = NULL, annotateTop = TRUE, ...) 
# {
#   CHR = BP = P = index = NULL
#   if (!(chr %in% names(x))) 
#     stop(paste("Column", chr, "not found!"))
#   if (!(bp %in% names(x))) 
#     stop(paste("Column", bp, "not found!"))
#   if (!(p %in% names(x))) 
#     stop(paste("Column", p, "not found!"))
#   if (!(snp %in% names(x))) 
#     warning(paste("No SNP column found. OK unless you're trying to highlight."))
#   if (!is.numeric(x[[chr]])) 
#     stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
#   if (!is.numeric(x[[bp]])) 
#     stop(paste(bp, "column should be numeric."))
#   if (!is.numeric(x[[p]])) 
#     stop(paste(p, "column should be numeric."))
#   if (!is.null(x[[snp]])) 
#     d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
#                    pos = NA, index = NA, SNP = x[[snp]], stringsAsFactors = FALSE)
#   else d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]], 
#                       pos = NA, index = NA)
#   d <- d[order(d$CHR, d$BP), ]
#   if (logp) {
#     d$logp <- -log10(d$P)
#   }
#   else {
#     d$logp <- d$P
#   }
#   d$index = rep.int(seq_along(unique(d$CHR)), times = tapply(d$SNP, 
#                                                              d$CHR, length))
#   nchr = length(unique(d$CHR))
#   if (nchr == 1) {
#     d$pos = d$BP
#     xlabel = paste("Chromosome", unique(d$CHR), "position")
#   }
#   else {
#     lastbase = 0
#     ticks = NULL
#     for (i in unique(d$index)) {
#       if (i == 1) {
#         d[d$index == i, ]$pos = d[d$index == i, ]$BP
#       }
#       else {
#         lastbase = lastbase + max(d[d$index == (i - 
#                                                   1), "BP"])
#         d[d$index == i, "BP"] = d[d$index == i, "BP"] - 
#           min(d[d$index == i, "BP"]) + 1
#         d[d$index == i, "pos"] = d[d$index == i, "BP"] + 
#           lastbase
#       }
#     }
#     ticks <- tapply(d$pos, d$index, quantile, probs = 0.5)
#     xlabel = "Chromosome"
#     labs <- unique(d$CHR)
#   }
#   xmax = ceiling(max(d$pos) * 1.03)
#   xmin = floor(max(d$pos) * -0.03)
#   def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", 
#                    las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(0, 
#                                                                      ceiling(max(d$logp))), xlab = xlabel, ylab = expression(-log[10](italic(p))))
#   dotargs <- list(...)
#   do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% 
#                                             names(dotargs)]))
#   if (!is.null(chrlabs)) {
#     if (is.character(chrlabs)) {
#       if (length(chrlabs) == length(labs)) {
#         labs <- chrlabs
#       }
#       else {
#         warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
#       }
#     }
#     else {
#       warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
#     }
#   }
#   if (nchr == 1) {
#     axis(1, ...)
#   }
#   else {
#     axis(1, at = ticks, labels = labs, ...)
#   }
#   col = rep_len(col, max(d$index))
#   if (nchr == 1) {
#     with(d, points(pos, logp, pch = 20, col = col[1], ...))
#   }
#   else {
#     icol = 1
#     for (i in unique(d$index)) {
#       points(d[d$index == i, "pos"], d[d$index == i, "logp"], 
#              col = col[icol], pch = 20, ...)
#       icol = icol + 1
#     }
#   }
#   if (suggestiveline) 
#     abline(h = suggestiveline, col = "blue")
#   if (genomewideline) 
#     abline(h = genomewideline, col = "red")
#   if (!is.null(highlight)) {
#     if (any(!(highlight %in% d$SNP))) 
#       warning("You're trying to highlight SNPs that don't exist in your results.")
#     d.highlight = d[which(d$SNP %in% highlight), ]
#     with(d.highlight, points(pos, logp, col = "#ff7f0e", 
#                              pch = 20, ...))
#   }
#   if (!is.null(annotatePval)) {
#     if (logp) {
#       topHits = subset(d, P <= annotatePval)
#     }
#     else topHits = subset(d, P >= annotatePval)
#     par(xpd = TRUE)
#     if (annotateTop == FALSE) {
#       if (logp) {
#         with(subset(d, P <= annotatePval), textxy(pos, 
#                                                   -log10(P), offset = 0.625, labs = topHits$SNP, 
#                                                   cex = 0.45), ...)
#       }
#       else with(subset(d, P >= annotatePval), textxy(pos, 
#                                                      P, offset = 0.625, labs = topHits$SNP, cex = 0.45), 
#                 ...)
#     }
#     else {
#       topHits <- topHits[order(topHits$P), ]
#       topSNPs <- NULL
#       for (i in unique(topHits$CHR)) {
#         chrSNPs <- topHits[topHits$CHR == i, ]
#         topSNPs <- rbind(topSNPs, chrSNPs[1, ])
#       }
#       if (logp) {
#         textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, 
#                labs = topSNPs$SNP, cex = 0.5, ...)
#       }
#       else textxy(topSNPs$pos, topSNPs$P, offset = 0.625, 
#                   labs = topSNPs$SNP, cex = 0.5, ...)
#     }
#   }
#   par(xpd = FALSE)
# }

##########################################################################################################################

# pdf(file = file.path(output.dir, "Results/002 Differential Methylation Analysis DMA/Manhattan plot.pdf"))
# MyManhattan(Results_data, ylim=c(0,10),
#             col = c("gray60", "gray80"),
#             highlight = rownames(tT2),
#             genomewideline = -log(1.37E-05, base = 10),
#             suggestiveline = -log(0.05, base = 10),
#             annotatePval = 1.569234e-07
# )#, main = "Meta_analysis 4 cohorts")) 
# dev.off()


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
AnnoThreshold <- -log(max(head(tT[order(tT$P.Value), "P.Value"], n=topSnps)), base = 10)
# Filtering the CpGs that should be highlighted in the figure
snpsOfInterest <- rownames(tT2)
# Nominal suggestive significance line
suggestiveLine <- -log(signThreshold, base = 10)
# Genomewide significance line
genomewideLine <- -log(max(tT[tT$adj.P.Val<signThreshold, "P.Value"]), base = 10)

upperLimit <- roundUp(-log(min(tT[tT$adj.P.Val<signThreshold, "P.Value"]), base = 10), to = 10)

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
  scale_y_continuous(limits = c(0,upperLimit), expand = c(0, 0), breaks = seq(0, upperLimit, 2)) +     # remove space between plot area and x axis
  
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
  # geom_hline(yintercept = suggestiveLine, colour = "forestgreen") + 
  geom_hline(yintercept = genomewideLine, colour = "black", linetype = "dashed") +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


ggsave(plot = manhattan_plot, 
       filename = file.path(output.dir, "Manhattan_AsthmaVsHealthy.tiff"),
       dpi = 300)
