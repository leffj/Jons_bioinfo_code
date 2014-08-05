# R script to take a distance matrix and average distances by a given factor provided in a mapping file

library(vegan)

# Take mean or median dissimilarities by indicated factor.

# input values
dmat.f = ''
map.f = ''
cat = ""
fun = mean
out.f = ''

# import
dmat = read.table(dmat.f,header=TRUE,row.names=1,sep='\t',check.names=FALSE,comment.char="")
map = read.table(map.f,header=TRUE,row.names=1,sep='\t',check.names=FALSE,comment.char="")

# convert dmat to 3 column format
dmat.clmns = data.frame(t(combn(unlist(labels(dmat)[1]),2)),as.numeric(as.dist(dmat)))
names(dmat.clmns) = c('x1','x2','dist')

# list x1 and x2 categories in new clmns
cat1 = map[match(dmat.clmns$x1,row.names(map)),cat]
cat2 = map[match(dmat.clmns$x2,row.names(map)),cat]
dmat.clmns = cbind(dmat.clmns,cat1,cat2)

# only take samples in mapping file
dmat.clmns = dmat.clmns[!is.na(dmat.clmns$cat1) & !is.na(dmat.clmns$cat2),]

# retain rows where distances are comparing samples from the same cat
dmat.clmns.reduced = dmat.clmns[dmat.clmns$cat1 != dmat.clmns$cat2,]

# combine cat columns to form comparison categories
# however, need to realize 'cat1__cat2' is equal to 'cat2__cat1'
# make df with every unique combination of cat levels
uniqueLevels = levels(factor(c(as.character(dmat.clmns.reduced$cat1),as.character(dmat.clmns.reduced$cat2))))
categories = data.frame(t(combn(uniqueLevels,2)))
# make a vector with each combination (first,second and second,first) for each combination
categorieslookup = paste(categories[,1],categories[,2])
categorieslookup = c(categorieslookup,paste(categories[,2],categories[,1]))
# make a vector with an index number for each combination of cat levels that associates with lookup vec
catIndex = c(seq.int(1,nrow(categories)),seq.int(1,nrow(categories)))
# now match the category combination (cat1 + cat2) in the dataset to the combination in 
#    the lookup vec and pull the category index
joinedCats = paste(dmat.clmns.reduced$cat1,dmat.clmns.reduced$cat2)
distCats = catIndex[match(joinedCats,categorieslookup)]
dmat.clmns.reduced = cbind(dmat.clmns.reduced,distCats)

# get measure of center for each distance category
byCatDists = tapply(dmat.clmns.reduced$dist,INDEX=dmat.clmns.reduced$distCats,FUN=fun)
byCatDists.clms = cbind(as.character(categories[,1]),as.character(categories[,2]),as.vector(byCatDists))
row.names(byCatDists.clms) = NULL

# convert 3 column distances back to matrix format
uNames = as.character(uniqueLevels)
byCatDists.dmat = data.frame(matrix(ncol=length(uNames),nrow=length(uNames)))
names(byCatDists.dmat) = uNames
row.names(byCatDists.dmat) = uNames
for(i in 1:nrow(byCatDists.clms)){
  byCatDists.dmat[byCatDists.clms[i,1],byCatDists.clms[i,2]] = as.character(byCatDists.clms[i,3])
}
for(i in 1:nrow(byCatDists.dmat)){
  for(k in 1:ncol(byCatDists.dmat)){
    if(is.na(byCatDists.dmat[i,k])){
      byCatDists.dmat[i,k] = byCatDists.dmat[k,i]
    }
  }
}
byCatDists.dmat[is.na(byCatDists.dmat)] = 0

# write output matrix
write.table(byCatDists.dmat,file=out.f,sep='\t',row.names=TRUE,col.names=NA)
