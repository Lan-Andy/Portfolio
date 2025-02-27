.libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
setwd("C:/Users/24976197/OneDrive - UTS/Desktop/R/")
root.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS"
expression.dir <- file.path(root.dir, "Data/004 Differential Gene Expression Analysis DGA")
## This script will perform DGE analysis for ATLANTIS dataset comparing asthma vs healthy subjects 
library(dplyr)
library(edgeR)
library(ggplot2)
library(tibble)
library(biomaRt)
library(conflicted)

# metafile 
# New meta.data 
pheno.data <- read.csv(file.path(root.dir, "atlantis_patient_data.csv"), header = TRUE, row.names = 1)
masterTable <- read.csv(file.path(root.dir, "methID_rnaID_clinID.csv"), header = TRUE) %>%
  dplyr::filter(RNA_qc_check == "YES")

masterTable$VISIT <- gsub(pattern = "Visit 1a", replacement = "VISIT 1", masterTable$VISIT)
masterTable$VISIT <- gsub(pattern = "Visit 2", replacement = "VISIT 2", masterTable$VISIT)

pheno.data.filt <- pheno.data %>%
  semi_join(masterTable, by = c("PT", "VISIT"))

#### Lancet Respiratory Medicine paper, 2019
#### "Shared and Distinct Genetic Risk Factors for Childhood Onset and 
####    Adult Onset Asthma: Genome- and Transcriptome-wide Studies"
####  https://doi.org/10.1016/S2213-2600(19)30055-4
#### They have defined:
####    Childhood onset < 12 years 
####    Adulthood onset > 25 and < 66
#### JACI paper, 2024
#### "Blood eosinophils and fractional exhaled nitric oxide are prognostic and 
####    predictive biomarkers in childhood asthma"
#### https://doi.org/10.1016/j.jaci.2023.09.044
#### Definition of childhood onset asthma beween 6-11 years
#### European respiratory review (ERS), 2013
#### "Adult-onset asthma: is it really different?"
#### https://doi.org/10.1183/09059180.00007112
#### Mentioning late onset asthma from 12-65 years
#### Thorax, 2016
#### "Clinical and functional differences between early-onset and late-onset adult asthma: 
####  a population-based Tasmanian Longitudinal Health Study"
#### https://doi.org/10.1136/thoraxjnl-2015-208183
#### Mentioning that age of 13 was their cut-off  
#### NEJM, 2003
#### "A longitudinal, population-based, cohort study of childhood asthma 
####  followed to adulthood"
#### https://doi.org/10.1056/nejmoa022363 
#### Mentioning that their childhood to adulthood asthma is from 9 - 26
#### JACI, 2018 Ian and Kian Fan are on this paper
#### "Pathway discovery using transcriptomic profiles in adult-onset severe asthma"
#### https://doi.org/10.1016/j.jaci.2017.06.037 
#### Mentioning cut-off 18 year for EO and AO

# Early onset asthma = < 13 yr diagnosed  | EO
# Adult onset asthma = >= 25 - 65         | AO

pheno.data.H <- pheno.data.filt %>%
  dplyr::filter(ASTHEA == "H")
pheno.data.H$onset <- "H"

pheno.data.EO <- pheno.data.filt %>%
  dplyr::filter(AGE_DIAG < 13)
pheno.data.EO$onset <- "EO"

pheno.data.AO <- pheno.data.filt %>%
  dplyr::filter(AGE_DIAG >= 25 & AGE_DIAG < 66)
pheno.data.AO$onset <- "AO"

# pheno.data.LO <- pheno.data.filt %>%
#   dplyr::filter(AGE_DIAG >= 40)
# pheno.data.LO$onset <- "LO"


# Lost 1 patient, due to no previous asthma diagnosis

common_PT_VISIT <- rbind(pheno.data.H[, c("PT", "VISIT", "onset", "AGE_DIAG")],
                         pheno.data.EO[, c("PT", "VISIT", "onset", "AGE_DIAG")],
                         pheno.data.AO[, c("PT", "VISIT", "onset", "AGE_DIAG")]
                         )

masterTable.filt <- masterTable %>%
  inner_join(common_PT_VISIT, by = c("PT", "VISIT"))

masterTable.filt <- masterTable.filt %>%
  group_by(GenomeScan_ID) %>%
  mutate(preferred_visit = ifelse(any(VISIT == "VISIT 1"), "VISIT 1", "VISIT 2")) %>%
  ungroup() %>%
  dplyr::filter(VISIT == preferred_visit)

rownames(masterTable.filt) <- masterTable.filt$GenomeScan_ID

# UMI deduplicated read count table
expression.data <- read.csv(file.path(expression.dir, "/20201107_ATLANTIS_raw_readcount_dedup_FINAL.csv"), header =TRUE) %>%
  tibble::column_to_rownames("Gene")
expression.data <- as.matrix(expression.data[, unique(masterTable.filt$GenomeScan_ID)])

# create design martix, covariates include sex, age, smoking status
design <- model.matrix(~0 + onset + age + GENDER + smoking.status, data = masterTable.filt)

# perform DE analysis
DGEL <- edgeR::DGEList(expression.data)
keep <- edgeR::filterByExpr(DGEL,design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
DGEL <- edgeR::estimateDisp(DGEL, design)
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) 
# legacy = TRUE: for edgeR 4.0 and above- to reproduce the results of edgeR 3.0)

# define contrast for the comparison
contrasts <- limma::makeContrasts(
  onset.status = onsetAO - onsetEO,
  levels = design
)
qlf   <- edgeR::glmQLFTest(fit, contrast = contrasts[,"onset.status"])

################################################################################

# create design martix, covariates include sex, age, smoking status
design <- model.matrix(~0 + AGE_DIAG + age + GENDER + smoking.status, data = masterTable.filt)

expression.data.filt <- expression.data[, rownames(design)]

colnames(expression.data.filt) %in% rownames(design)
# perform DE analysis
DGEL <- edgeR::DGEList(expression.data.filt)
keep <- edgeR::filterByExpr(DGEL,design) 
DGEL <- DGEL[keep, , keep.lib.sizes=FALSE]
DGEL <- edgeR::calcNormFactors(DGEL, method = "TMM")
DGEL <- edgeR::estimateDisp(DGEL, design)
fit <- edgeR::glmQLFit(DGEL, design, legacy = TRUE) 
# legacy = TRUE: for edgeR 4.0 and above- to reproduce the results of edgeR 3.0)

# define contrast for the comparison
contrasts <- limma::makeContrasts(age_effect = AGE_DIAG, 
                                  levels = design)
qlf2 <- edgeR::glmQLFTest(fit, 
                          contrast = contrasts[,"age_effect"])
################################################################################


# Add gene names
## add gene names 
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_ids <- rownames(DGEL$counts)  # here are the ensembl gene IDs
all_new_gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)

#create the result table
de.results <- edgeR::topTags(
  qlf,
  n=nrow(DGEL))$table %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::left_join(
    y = all_new_gene,
    by = c("Gene" = "ensembl_gene_id")) %>%
  mutate (logFC = round (logFC, 2),
          PValue = signif(PValue, digits = 3),
          FDR = signif(FDR, digits = 3))# %>%
#readr::write_csv("./Umi_dedup/Dif_expr/DE.genes.ATLANTIS.csv")

#summarize:
summary(decideTests(qlf)) 

# 1*onsetAO -1*onsetEO
# Down      7
# NotSig    17985
# Up        1



