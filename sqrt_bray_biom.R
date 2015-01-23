#!/usr/bin/Rscript

library(vegan)
require(biom)

# ./sqrt_bray.r <input_file> <output_file>

otuTable.f=commandArgs()[6]
outFP=commandArgs()[7]

otuTable = read_biom(otuTable.f)
otuTable = as.data.frame(as.matrix(biom_data(otuTable)))
# transform otu table (square root transformation)
otuTable.xform = t(sqrt(otuTable))
# create dissimilarity matrix from otu table
otuTable.dist = vegdist(otuTable.xform,method='bray')
# write dissimilarity matrix to file
write.table(as.matrix(otuTable.dist),file=outFP,sep='\t',row.names=TRUE,col.names=NA)
