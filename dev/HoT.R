library("ape")
library("phylobase")
library("adephylo")

# generate random tree
plot(tree <- rtree(15))

partitions <- function(tree){
# get all partitions of the tree
partition <- treePart(tree, "dummy")
partition

# remove duplicated (first bifurcation from node) since its duplicated
# if the first split includes a single sequences, then do not take it as first
# comparision
rs <- rowSums(partition)
if(any(rs==0)) {
  part2 <- partition[!(rs==0),]
}

# find dupliate and remove from partitions
dup <- list()
res <- list()
for(j in 1:ncol(part2)){
  part <- part2[,-j]
  for(i in 1:ncol(part)){
    dup[[i]] <- (part2[,j]>0) == !(part[,i]>0)
  }
  res[[j]] <- dup
}
dup <- lapply(res, function(x) do.call(cbind, x))
del <- unlist(lapply(dup, function(x) which(apply(x, 2, all))))
partition <- partition[-del[2]]
return(partition)
}


parts_for_MSA <- partitions(tree)


