#!/usr/bin/Rscript

# ./sqrt_bray.r <input_file> <output_file>

inFP=commandArgs()[6]
outFP=commandArgs()[7]

# inFP = '/Users/leffj/software/Jons_bioinfo_code/test_gc.txt'

classes = c('character','integer','numeric')

input = read.table(inFP,header=FALSE,sep='\t',colClasses=classes,comment.char="")
means = tapply(input[,3],INDEX=input[,1],FUN=mean)
stds = tapply(input[,3],INDEX=input[,1],FUN=sd)
smry = data.frame(means,stds)

write.table(smry,file=outFP,sep='\t',row.names=TRUE,col.names=NA,quote=F)
