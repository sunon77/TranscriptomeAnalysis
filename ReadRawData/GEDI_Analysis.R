# TCGA Breast cancer RNASeq Gene expression data from Theo
# Breast cancer cell clustering together
# totally 1040 Breast Tissue Samples
# - 7 metastasis gene expression
# - 111 normal samples
# 
# Author: Joseph X. Zhou
# Date:   1/23/2015
# Modify History:
# 1) Date: 06/04/2014,

## Linux server: osiris
cd /proj/huanglab/orthomcl/CancerExpEnrichAnalysis

ptm <- proc.time();
#setwd("C:/Users/jx.zhou/OneDrive/Research/Projects/CancerIntercellNet")
setwd("/proj/huanglab/orthomcl/CancerExpEnrichAnalysis")
inFileName = "BCDiff_15400x1165_TCGA_BC_RNASeq.txt"
geneIn= read.csv(inFileName, sep="\t",header=T)
nCol = dim(geneIn)[2];
head(geneIn[,1:12])
print(dim(geneIn))
rowNum = dim(geneIn)[1]
colNum = dim(geneIn)[2]

## Rename colums
for(i in 2:colNum)
{
	colnames(geneIn)[i] = paste(as.character(geneIn[1,i]),as.character(colnames(geneIn)[i]),"");
}
head(geneIn)



##  Quantile normalizing after filtering background
library(preprocessCore)
geneExpMat = as.matrix(geneIn[3:15402,2:1166])
geneExpMat_QN = quantile(geneExpMat, na.rm = TRUE)



## outFileName = "BCDiff_1404x47_14CTL_SAM_MultiClass_cutoff30.txt"
outFileName = "geneExpMatrix_GEDI.txt"
write.table(geneExpMat_QN, file=outFileName, sep=",", 
              col.name=T, row.names=F, quote=FALSE)
