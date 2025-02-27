# .libPaths("C:/Users/24976197/OneDrive - UTS/Desktop/R/Packages")
# output.dir <- "C:/Users/24976197/OneDrive - UTS/Desktop/PhD UTS Data folder/004 Bioinformatics/ATLANTIS/"

##############################################################################
################## Load used libraries within script #########################
##############################################################################
library("MatrixEQTL")
library("limma")
library("dplyr")
##############################################################################
#################### Load in raw data for analysis ###########################
##############################################################################
#---NOTES---
#The  code performs eQTL analysis of  data set consisting of three files: genotype, expression, and covariates. Also have gene location and SNP location file
#For every gene-SNP pair it runs linear regression analysis accounting for the set of covariates.
#GE, SNP and Covariates data files must have columns corresponding to samples and with one gene/SNP/covariate in each row. 
#The columns of all three files must have matching order. 
#All measurements must be numeric and the values in the genotype data set do not have to be descrete.
#-----------

## Location of the folder with the data files.
base.dir = "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_eQTM/01_data_prep/"
output.dir <- "/groups/umcg-griac/tmp02/projects/AndyLan/ATLANTIS_eQTM/02_eQTM_results/"

useModel = modelLINEAR;


################################################################################
# #### Asthma vs healthy #### 
# # Genotype file name
# SNP_file_name = paste(base.dir, "methylationData_eQTM_samples_asthma_vs_control.csv", sep="");
# snps_location_file_name = paste(base.dir, "CpGLocation_asthma_vs_control.csv", sep="");
# 
# # Gene expression file name
# expression_file_name = paste(base.dir, "normalized_geneExpressionData_eQTM_asthma_vs_control_samples.csv", sep="");
# gene_location_file_name = paste(base.dir, "geneLocation.csv", sep="");
# 
# # Covariates file name
# # Set to character() for no covariates
# covariates_file_name = paste(base.dir, "covariates_eQTM_samples_asthma_vs_control.csv", sep="");
################################################################################
# #### Th2-high vs healthy #### 
# Genotype file name
# SNP_file_name = paste(base.dir, "methylationData_eQTM_AllSamples_Th2HighVsHealthy.csv", sep="");
# snps_location_file_name = paste(base.dir, "CpGLocation_Th2HighVsHealthy.csv", sep="");
# 
# # Gene expression file name
# expression_file_name = paste(base.dir, "normalized_geneExpressionData_eQTM_Th2HighLowUndeterm_samples.csv", sep="");
# gene_location_file_name = paste(base.dir, "geneLocation.csv", sep="");
# 
# # Covariates file name
# # Set to character() for no covariates
# covariates_file_name = paste(base.dir, "covariates_eQTM_samples_Th2HighLowUndeterm.csv", sep="");
################################################################################
#### Th2-high vs Th2-low ####
# Genotype file name
SNP_file_name = paste(base.dir, "methylationData_eQTM_AllSamples_Th2HighvsTh2Low.csv", sep="");
snps_location_file_name = paste(base.dir, "CpGLocation_Th2HighvsTh2Low.csv", sep="");

# Gene expression file name
expression_file_name = paste(base.dir, "normalized_geneExpressionData_eQTM_Th2HighLowUndeterm_samples.csv", sep="");
gene_location_file_name = paste(base.dir, "geneLocation.csv", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "covariates_eQTM_samples_Th2HighLowUndeterm.csv", sep="");
################################################################################

# snps <- read.csv(SNP_file_name, header= TRUE, row.names = 1)
# head(snps)
# snpspos <- read.csv(snps_location_file_name, header= TRUE, row.names = 1)
# head(snpspos)
# gene <- read.csv(expression_file_name, header= TRUE, row.names = 1)
# head(gene)
# genepos <- read.csv(gene_location_file_name, header= TRUE, row.names = 1)
# head(genepos)
# quit()

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0;#1e-2;

# Error covariance matrix
errorCovariance = numeric();

# Distance for local gene-SNP pairs
cisDist = 50000;

## Load genotype data
print("Genotype")
snps = SlicedData$new();
snps$fileDelimiter = ",";      # the TAB character, instead of /t it was just spaces (could have saved with /t as delimeter by adding ,sep ="/t" when using write.table earlier)
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data
print("Gene expression")
gene = SlicedData$new();
gene$fileDelimiter = ",";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates
print("Covariates")
cvrt = SlicedData$new();
cvrt$fileDelimiter = ",";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 2;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

#Run the analysis
snpspos = read.csv(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1);
genepos = read.csv(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 1);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot",
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);


unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)
me$cis$eqtls %>%
  readr::write_csv(
    file.path(output.dir, "all-after-matrixeqtl_Th2HighvsTh2Low_AllSamples.csv")
  )


####EQTM MULTIPLE TESTING CORRECTION ####

# CISeqtls_sig <- read.csv(file.path(output.dir, "all-after-matrixeqtl.csv"))
CISeqtls_sig <- me$cis$eqtls

list=as.matrix(unique(CISeqtls_sig$snps))

finaltable=matrix(ncol=7,nrow=0)
colnames(finaltable)=c(colnames(CISeqtls_sig),"finalFDR")

for (i in 1:nrow(list)){
  cpg=list[i,]
  sub=CISeqtls_sig[which(CISeqtls_sig$snps==cpg),]
  finalFDR=p.adjust(sub$pvalue, method="BH")
  finaltable=rbind(finaltable,cbind(sub,finalFDR))
  print(i/nrow(list)*100)
}

#plot(me)

#SAVE ALL EQTMs
write.csv(finaltable, file = file.path(output.dir, "all_genomewideEQTMS_Th2HighvsTh2Low_AllSamples.csv"))

##SAVE ALL SIGNIFICANT EQTMS WITH MULTIPLE TESTING CORRECTION
finaltable_sig <- finaltable[which(finaltable$finalFDR<0.05),]
finaltable_sig <- finaltable_sig[order(finaltable_sig$finalFDR),]
write.csv(finaltable_sig, file = file.path(output.dir, "final_genomewideEQTMS_Th2HighvsTh2Low_AllSamples.csv"))

save.image(file = file.path(output.dir, "toEQTM_Th2HighvsTh2Low_AllSamples.rdata")) #ALL DATA UP TO HERE

# finaltable_sig_bak <- finaltable_sig
#This is to only look at the sites that were upregulated
#intersect(finaltable_sig$gene,row.names(tT2_1_DEA[which(tT2_1_DEA$logFC>0),]))  
# ditto, downregulated
#intersect(finaltable_sig$gene,row.names(tT2_1_DEA[which(tT2_1_DEA$logFC<0),]))

# save.image("Workspaces/toEQTM.rdata") #ALL DATA UP TO HERE
# load("Workspaces/toEQTM.rdata")
# beepr::beep(sound = 2)