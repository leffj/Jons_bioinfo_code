#!/usr/bin/Rscript

library(vegan)

# ./raup-crick.r <input_file> <output_file>

otuTable.f=commandArgs()[6]
outFP=commandArgs()[7]

otuTable.wTaxa = read.table(otuTable.f,header=TRUE,row.names=1,sep='\t',check.names=FALSE,skip=1,comment.char="")
# strip taxonomy
taxonomy = otuTable.wTaxa[,ncol(otuTable.wTaxa)]
otuTable = otuTable.wTaxa[,1:(ncol(otuTable.wTaxa)-1)]
# transpose otu table
otuTable.xform = t(otuTable)
# create dissimilarity matrix from otu table
otuTable.dist = vegdist(otuTable.xform,method='raup')
# write dissimilarity matrix to file
write.table(as.matrix(otuTable.dist),file=outFP,sep='\t',row.names=TRUE,col.names=NA)
