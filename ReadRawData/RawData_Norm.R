# Read Raw microarray data and process them to gene expression Matrix
# Author: Joseph X. Zhou
# Date:   05/06/2015
# Modify History:
# 1) Date: 11/05/2013, 

read.chunk <- function(file, readLines, skipLine){
	if(skipLine > 1){
		p = read.csv(file, skip = skipLine,
	   #p = read.csv(file, skip = (lines*(clump-1))+1 if not a textConnection
			nrows = readLines, sep="\t",header=TRUE)
	} else {
		p = read.csv(file, skip = 0, nrows = lines)
	}
	# prepare column number
	# GeneID and ProbID - column 1,2
	selectCol = c(1,2);
	#First signal colume, first p-value,
	iFirstSignalCol = 3
	iFirstpValueCol = 8
	iSpaceCol       = 6
	iTotalSampleNum = 72
	for (i in 1:iTotalSampleNum)
	{
	  iCol = c(iFirstSignalCol+(i-1)*iSpaceCol, iFirstpValueCol+(i-1)*iSpaceCol);
	  selectCol =  c(selectCol, iCol);
	}
	# Last three GO Term
	iLastCol = length(p)
	return(p[selectCol])
}

setwd("C:/Users/jx.zhou/Dropbox/Programming/R/MCF_Reprg_GeneAnalysis/data")
computerName = Sys.getenv('computername') 
inFileName  = "WeiWu72samples_RawData_Feb2012.txt"
# Title lines to be skipped
skipLine  = 8;
# How many lines to read   ( stop line = 29383  - 8 = 29375
readLines = 29375;
geneIn = read.chunk(inFileName, readLines, skipLine)


################################################################
## Direct read
inFileName  = "WeiWu72samples_RawData_Feb2012.txt"
geneIn= read.csv(inFileName, sep="\t",header=T)
print(head(geneIn))
print(dim(geneIn))


# Look at gene expression distribution and p values distribution
require(graphics)
# before normalization
boxplot(geneIn[,3:72], what=c(1,1,1,0), log="",main="Raw Data", 
	   xlab = "Different samples", ylab = "gene expression") #,ylim=c(0,6e3))

#####################################################
# Quantile normalizing after filtering background
library(preprocessCore)
geneMat = data.matrix(geneIn[,3:62])
geneMatNorm = normalize.quantiles(geneMat)
# before normalization
boxplot(geneMatNorm, what=c(1,1,1,0), log="",main="Raw Data", 
	   xlab = "Different samples", ylab = "gene expression") 

geneOutGediNorm = geneIn
geneOutGediNorm[,3:62] = data.frame(geneMatNorm)

#####################################################
# median normalizing after 
ptm <- proc.time()
library(limma)
geneMat = data.matrix(geneIn[,3:62])
#geneMat = replicate(11,rnorm(105,0,100))
colMedian= apply(geneMat,2,median)
matMedian = median(colMedian)
geneMatNormMedian = apply(geneMat,2,function(col) col/median(col)*matMedian);


# geneMatNormMedian = maNorm(geneMat,norm="median", echo=TRUE)
geneOutGediNorm = geneIn;
geneOutGediNorm[,3:62] = data.frame(geneMatNormMedian)

outFileName = "BCDiff_13955x60_ColOrdered_medianNormalized.txt"
write.table(geneOutGediNorm, file=outFileName, sep="\t", row.names=F, quote=FALSE)
print(proc.time()-ptm);


