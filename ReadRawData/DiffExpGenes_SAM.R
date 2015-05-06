# From Gene Exp Profile after removing background noise using at least one column P < 0.05
# Use SAM, LIMMA and RankPod to find differentially expressed genes in Yan Cui's data Set
# http://compbio.uthsc.edu/microarray/lecture2.htm
# Author: Joseph X. Zhou
# Date:   03/02/2013
# Modify History:
# 1) Date:  03/25/2013, really test 3 methods for diff
# 2) Date:  03/27/2013, Find commonly diffExpGene for BCD 
# 3) Date:  11/14/2013, LIMMA DEGs for 1DayP and 2DayP  ###

ptm <- proc.time()
setwd("C:/Users/joseph/Dropbox/Research/Projects/BreastCancerDiff/Data")

# Quantile normalize first, then 

inFileName = "BCDiff_13956x60_ColNameOrder_bk1DayP_QuntileNormalized.txt";
fileStr = strsplit(inFileName, "\\.")
geneIn= read.csv(inFileName, sep="\t",header=TRUE,row.names=NULL)


#######################################################
## Define classes
dayGroup        = rep(2,1,47); # 67: remove DMSO
dayGroup[1:15]  = 0;
dayGroup[(16:(15+17))]   = 1;
print(dayGroup)

#######################################################
## SAM Analysis: Multi-classes: Untreated, Day 1, 3, 5
library(siggenes);
geneclass = dayGroup;
genenames = geneIn[,1] 
geneExpMat = geneIn[,c(3:49)]


ptm <- proc.time();
sam.out <-sam(geneExpMat,geneclass, B=1000, rand = 923, 
  gene.names = genenames, na.replace=TRUE)
cutoff    = 80;
sum.sam.out = summary(sam.out, cutoff);

#print(sum.sam.out)
sigGenes = sum.sam.out@mat.sig;
# examine significant gene list
#plot(sam.out,cutoff) #multi class, plot has no meaning
print(cutoff)
print(dim(sigGenes))
print(proc.time()-ptm);

#print(sigGenes[,c(5,6)])

sigBool     <- geneInDayOrdered[,1] %in% rownames(sigGenes)
genesChangedSig <- subset(geneInDayOrdered, sigBool)
dim(genesChangedSig)

## outFileName = "BCDiff_1404x47_14CTL_SAM_MultiClass_cutoff30.csv"
outFileName = "BCDiff_138x47_14CTL_SAM_MultiClass_cutoff90.csv"
write.table(genesChangedSig, file=outFileName, sep=",", 
              col.name=T, row.names=F, quote=FALSE)
