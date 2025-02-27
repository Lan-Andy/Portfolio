# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
# output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/"

##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("dplyr")
library("limma")
# library("edgeR")
# library("DESeq2")
# library("biomaRt")

##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
input.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_methylation"
output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_eQTM/01_data_prep"


linkingtable <- read.csv("/groups/umcg-griac/tmp02/projects/TatianaKarp/ATL_methyaltion/methID_rnaID_clinID.csv")
pheno.data <- read.csv("/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_gene_expression/data/atlantis_patient_data.csv",
                       header = TRUE,
                       row.names = 1,
                       na.strings = ""
                       )
GE.data <- read.csv("/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_gene_expression/data/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv",
                    header = TRUE,
                    row.names = 1
                    )
# meth.data <- as.matrix(
#     read.csv(file = file.path(input.dir, "04_normalisation/ATLANTIS_betas_dasen.csv"),
#                header = TRUE,
#                row.names = 1)
#     )
load(file = file.path(input.dir, "06_DMA/initialData.RData")) #methylationData
meth.data <- methylationData
rm(methylationData)
rownames(meth.data) <- meth.data[, 1]
meth.data <- meth.data[, -1]
meth.data <- as.matrix(meth.data)

PC.all.cp <- read.csv(file = file.path(input.dir, "05_pca/ctrlProbes/my_pca_data.csv"),
                      header = TRUE)
CpG.sites <- read.csv(file = file.path(input.dir, "06_DMA/DMA_Tt2_significant_CpG_Th2HighvsTh2Low_BH.csv"),
                      header = TRUE,
                      row.names = 1)
# CpG.annotation <- read.csv(file = file.path(input.dir, "07_anno/CpG_annotation.csv"),
#                         header = TRUE,
#                         row.names = 1)
Th2_groups <- read.csv(file = file.path(output.dir, "ATLANTIS_linkingtable_with_group_Th.csv"),
                       header = TRUE,
                       row.names = 1)



# linkingtable <- read.csv(file.path(output.dir, "methID_rnaID_clinID.csv"))
# 
# pheno.data <- read.csv(file.path(output.dir, "atlantis_patient_data.csv"),
#                        header = TRUE, row.names = 1, na.strings = "")
# 
# GE.data <- read.csv(file.path(output.dir, "Data/004 Differential Gene Expression Analysis DGA/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"),
#                     header = TRUE,
#                     row.names = 1)
# meth.data <- as.matrix(
#   read.csv(file = file.path(output.dir, "Data/sliced_methylationData_Mvalues.csv"),
#              header = TRUE,
#              row.names = 1
#              )
# )
# 
# PC.all.cp <- read.csv(file = file.path(output.dir, "ATLANTIS QC/Data/PCA technical probes/my_pca_data.csv"),
#                       header = TRUE)
# 
# CpG.sites <- read.csv(file = "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/Data/002 Differential Methylation Analysis DMA/Correct for PC1-3/Tt2_significant_CpG_asthma_vs_control_BH.csv",
#                       header = TRUE,
#                       row.names = 1)
# Th2_groups <-read.csv(file.path(output.dir, "Data/100 Characteristics tables/ATLANTIS_linkingtable_with_group_Th.csv"),
#                       header = TRUE,
#                       row.names = 1)

##############################################################################
#################### Get a list of locations for the genes and CpG ###########
##############################################################################
# Get the gene locations, by extracting HGNC symbol, ENSEMBL id, start and end position
# # Get the ensembl information of Homo sapiens
# ensembl <- useEnsembl(biomart='ensembl',
#                       dataset="hsapiens_gene_ensembl")
# 
# # Check what version of ensembl we using, what build it has build upon
# searchDatasets(mart = ensembl, pattern = "hsapiens") # version GRCh38.p14
# 
# # Extract the gene symbol, ensembl id, chromosome, start and end position fo the genes
# # Deleted the hgnc symbol, as not all ensembl id genes have a gene symbol 13/07/2024
# results <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
#                  filters = c("chromosome_name"),
#                  values = c(1:22,"X","Y"),
#                  mart = ensembl)
# # Adding "chr" prefix to chromosome numbers/X/Y
# results$chromosome_name <- paste0("chr", results$chromosome_name)
# # Save the gene locations for eQTM analysis
# write.csv(x = results, file = file.path(output.dir, "geneLocation.csv"))
# rm(ensembl, results)

# # Change the CpG annotation file to CpG location file ONLY
# CpG.annotation$cpg <- rownames(CpG.annotation)
# select.cols <- c("cpg", "chr", "pos")
# CpG.annotation <- CpG.annotation[, select.cols]
# 
# write.csv(x = CpG.annotation, file = file.path(output.dir, "cpgLocation.csv"))
# rm(CpG.annotation, select.cols)


##############################################################################
#################### Select and transform data.frames ########################
##############################################################################
# Change sample names (located in columns) to match the linking table
colnames(GE.data) <- gsub(pattern = "X", replacement = "", colnames(GE.data))
colnames(meth.data) <- gsub(pattern = "X", replacement = "", colnames(meth.data))
linkingtable$GenomeScan_ID <- gsub(pattern = "X", replacement = "", linkingtable$GenomeScan_ID)

# Change visit values to be more clear and easier access
linkingtable$VISIT <- gsub(pattern = "Visit 1a", replacement = "VISIT 1", x = linkingtable$VISIT)
linkingtable$VISIT <- gsub(pattern = "Visit 2", replacement = "VISIT 2", x = linkingtable$VISIT)

# Selecting samples that have patient id & RNAseq & methylation data
sum(is.na(linkingtable$PT)) # 53 methylation samples that do not have patient ID

# Select samples with patient id
linktable1 <- linkingtable[!is.na(linkingtable$PT) & linkingtable$VISIT == "VISIT 1",]

# Filter out samples that have bad quality RNAseq
linktable2 <- linktable1[linktable1$RNA_qc_check == "YES",]

# Identifying samples with only visit 2
unique.patient.ids <- linkingtable[!is.na(linkingtable$PT), ]
single.visit2.ids <- data.frame()
for (patient.id in unique.patient.ids$PT) {
  temp <- linkingtable[linkingtable$PT %in% patient.id, c("PT", "VISIT", "RNA_qc_check")]
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

# Filter out samples with only visit 2 and have good RNAseq check
single.visit2.ids <- single.visit2.ids[single.visit2.ids$RNA_qc_check == "YES", ]

# linktable3 consist of:
#   1. Samples that either have visit 1 or visit 1 and 2, and have a good RNAseq check
#   2. Samples that only have visit 2 methylation data and a good RNAseq check
linktable3 <- rbind(linktable2, 
                    linkingtable[linkingtable$PT %in% single.visit2.ids$PT &
                                   linkingtable$RNA_qc_check == "YES", ])

pd.methyl <- pheno.data[pheno.data$PT %in% linktable2$PT, ]
pd.methyl.v1 <- pd.methyl[pd.methyl$VISIT %in% "VISIT 1", ]

pd.methyl <- pheno.data[pheno.data$PT %in% single.visit2.ids$PT, ]
pd.methyl <- pd.methyl %>% tidyr::fill(SMOKE)
pd.methyl.v2 <- pd.methyl[pd.methyl$VISIT %in% "VISIT 2", ]

# pheno.sub consist of clinical data of samples that have
#   gene expression and methylation data
#   They have the corresponding pheno data with their visit
pheno.sub <- rbind(pd.methyl.v1, pd.methyl.v2)

# Selecting the PCs only from the samples we are interested in (see pd.methyl.sub patient.ids)
PC.all.cp <- PC.all.cp[PC.all.cp$meth_file_id %in% linktable3$meth_file_id, ]

# Merge all the dataframes to 1 big dataframe
pheno.sub <- merge(pheno.sub, linktable3, by = "PT")
pheno.sub <- merge(pheno.sub, PC.all.cp, by = "meth_file_id")

## Adding the group_th to the big pheno_data
pheno.sub <- merge(pheno.sub, Th2_groups, by = "meth_file_id")

# write.csv(x = pheno.sub, file = file.path(output.dir, "pheno_eQTM_samples.csv"))
# print("Printing colnames of pheno.sub")
# colnames(pheno.sub)

# Extract the variables only needed for the eQTM analysis and save
covariates <- c("PT.x", "GenomeScan_ID.x", "meth_file_id", 
                "ASTHEA.x", "AGE", "SEX", "SMOKE", "group_th",
                "PC1", "PC2", "PC3")

pheno.sub <- pheno.sub[, covariates]

covariates <- c("PT", "GenomeScan_ID", "meth_file_id", 
                "ASTHEA", "AGE", "SEX", "SMOKE", "group_th",
                "PC1", "PC2", "PC3")
colnames(pheno.sub) <- covariates
# write.csv(x = pheno.sub, file = file.path(output.dir, "covariates_eQTM_samples_IDs.csv"))


# Replace the values into integers
pheno.sub$SEX <- gsub(pattern = "M", replacement = 1, pheno.sub$SEX) # male
pheno.sub$SEX <- gsub(pattern = "F", replacement = 0, pheno.sub$SEX) # female
pheno.sub$SEX <- as.integer(pheno.sub$SEX)

pheno.sub$AGE <- as.integer(pheno.sub$AGE)

pheno.sub$ASTHEA <- gsub(pattern = "A", replacement = 1, pheno.sub$ASTHEA) # asthma
pheno.sub$ASTHEA <- gsub(pattern = "H", replacement = 0, pheno.sub$ASTHEA) # healthy
pheno.sub$ASTHEA <- as.integer(pheno.sub$ASTHEA)

pheno.sub$SMOKE <- gsub(pattern = "Non-smoker", replacement = 0, pheno.sub$SMOKE) # never-smoker
pheno.sub$SMOKE <- gsub(pattern = "Ex-smoker", replacement = 1, pheno.sub$SMOKE) # ex-smoker
pheno.sub$SMOKE <- gsub(pattern = "Current Smoker", replacement = 2, pheno.sub$SMOKE) # current smoker
pheno.sub$SMOKE <- as.integer(pheno.sub$SMOKE)


pheno.sub$group_th <- gsub(pattern = "healthy", replacement = 0, pheno.sub$group_th)
pheno.sub$group_th <- gsub(pattern = "low", replacement = 1, pheno.sub$group_th)
pheno.sub$group_th <- gsub(pattern = "high", replacement = 2, pheno.sub$group_th)
pheno.sub$group_th <- gsub(pattern = "undeterm", replacement = 3, pheno.sub$group_th)
pheno.sub$group_th <- as.integer(pheno.sub$group_th)


# Filter methylation data with eQTM samples & signifant CpG sites
meth.data.sub <- meth.data[rownames(CpG.sites), pheno.sub$meth_file_id]

# source(file.path(input.dir, "M2beta.R"))
# meth.data.sub <- beta2M(meth.data.sub)


# # Filter gene expression data with eQTM samples 
GE.data.sub <- GE.data[, pheno.sub$GenomeScan_ID]

# Normalize (log2cpm) the gene expression with edgeR 
#   Tatiana has done this similarly
# 2024-09-30: Alen disapproves this method of TMM
#   use limma::voom() as whole pipeline is CpG site related
GE.data.sub <- as.matrix(GE.data.sub)

# # Transform matrix into usable DGEList object for edgeR
# DGEL <- edgeR::DGEList(GE.data.sub)
# # Calculating the factors which should be normalized for
# # Chose for TMM : Trimmed Mean of the M-values
# DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
# # Perform log2normalization based on TMM method
# logcpm.ge.sub <- edgeR::cpm(DGEL, log=TRUE)


### Old way of filtering genes, only do this when performing 
###   differential gene expression.
# Calculate the total library sizes per sample
# colsums.GE.data.sub <- colSums(GE.data.sub)
# 
# # Calculate the FPM
# FPM.GE.data.sub <- t((t(GE.data.sub)/colsums.GE.data.sub) * 1E6)
# 
# # Check how the table looks like
# FPM.GE.data.sub[1:5, 1:5]
# range(colSums(FPM.GE.data.sub)) # 1e+06 1e+06
# 
# ## Select genes with mean >= 10FPM in complete group
# # Get the rowmeans of every row (=gene)
# rowmeans.FPM.GE.data.sub <- rowMeans(FPM.GE.data.sub)
# 
# # Filter the genes that have mean FPM >= 10, result = indices of the genes
# eQTM.selected.GE.data.sub <- GE.data.sub[which(rowmeans.FPM.GE.data.sub >= 10),]



# # Normalize the gene expression RNAseq data with limma::voom()
# normalized.GE.data.sub <- limma::voom(eQTM.selected.GE.data.sub)
# # Extract the normalized gene expression data
# normalized.GE.data.sub <- normalized.GE.data.sub$E

############################################################################
normalized.GE.data.sub <- limma::voom(GE.data.sub)
normalized.GE.data.sub <- normalized.GE.data.sub$E

# norm.ge.allgenes <- normalized.GE.data.sub
# norm.ge.subgenes <- normalized.GE.data.sub
# 
# plot(x = norm.ge.allgenes["ENSG00000074657", colnames(norm.ge.subgenes)],
#      y = norm.ge.subgenes["ENSG00000074657", ])
# cor.test(x = norm.ge.allsamples["ENSG00000074657", colnames(norm.ge.subsamples)],
#                                 y = norm.ge.subsamples["ENSG00000074657", ])
############################################################################

# Order the samples the same as covariates table
#   https://cran.r-project.org/web/packages/MatrixEQTL/MatrixEQTL.pdf
# See Matrix_eQTL_main() -> cvrt argument definition
print("Ordering gene expression and methylation data")
print(pheno.sub)
dim(pheno.sub)
class(pheno.sub)

meth.data.sub <- meth.data.sub[, pheno.sub$meth_file_id]
normalized.GE.data.sub <- normalized.GE.data.sub[, pheno.sub$GenomeScan_ID]
covariates <- c("PT",
                # "ASTHEA", # Use when asthma vs healthy
                "AGE", "SEX", "SMOKE", 
                "group_th", # Use when other comparison
                "PC1", "PC2", "PC3")
pheno.sub <- pheno.sub[, covariates]

# Transposing covariates as samples should be in columns for Matrix_eQTL_main()
pheno.sub <- as.data.frame(t(pheno.sub))

colnames(pheno.sub) <- paste0("ATLANTIS_", pheno.sub["PT", ])
pheno.sub <- pheno.sub[-which(rownames(pheno.sub) == "PT"), ]

print("Saving covariates, genexpression and methylationdata")
# Save covariates, methylation data, gene expression data for eQTM analysis
write.csv(x = pheno.sub, file = file.path(output.dir, "covariates_eQTM_samples_Th2HighLowUndeterm.csv"))

write.csv(x = meth.data.sub, file = file.path(output.dir, "methylationData_eQTM_AllSamples_Th2HighvsTh2Low.csv"))

write.csv(x = normalized.GE.data.sub, file = file.path(output.dir, "normalized_geneExpressionData_eQTM_Th2HighLowUndeterm_samples.csv"))





