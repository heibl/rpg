#'  Helper funcitons for HoT
#'  @export


# create N-3 partitions from tree
partitions <- function(tree){

  # get all partitions of the tree
  partition <- treePart(tree, "dummy")

  # remove duplicated partition (first bifurcation from node)
  # if the first split includes a single sequences, then do not take it as
  # first comparision

  partition <- partition[!(colSums(partition) %in% c(nrow(partition)-1, 1))]

  # find dupliate and remove from partitions
  ## [FK - 2017-03-27]
  ## this can be done faster if only the combinations are checked.
  dup <- list()
  res <- list()
  for(j in 1:ncol(partition)){
    part <- partition[,-j]
    for(i in 1:ncol(part)){
      dup[[i]] <- (partition[,j]>0) == !(part[,i]>0)
    }
    res[[j]] <- dup
  }
  dup <- lapply(res, function(x) do.call(cbind, x))
  del <- unlist(lapply(dup, function(x) which(apply(x, 2, all))))
  if(length(del)>0)
    partition <- partition[-del[2]]
  return(partition)
}


# Reverse DNA sequences
rev_DNA <- function(x){
  if(is.matrix(x)){
    x_mat <- as.character(x)
    rev <- t(apply(x_mat, 1, rev))
  }
  if(is.list(x)){
    x_l <- as.character(x)
    rev <- lapply(x_l, rev)
  }
  if(inherits(x, "DNAbin"))
    rev <- as.DNAbin(rev)
  if(inherits(x, "AAbin"))
    rev <- as.AAbin(rev)
  return(rev)
}


# Produce alternative MSAs from starting tree
align_part_set <- function(x, partition_set, mafft_exec, method){
  seq_left <- x[partition_set>0]
  seq_right <- x[partition_set==0]

  ## Aligning left HEADS and TAILS
  if (length(seq_left) == 1){
    headsA <-  seq_left
    tailsA <- rev_DNA(seq_left)
  }else{
    headsA <- mafft(seq_left, exec = mafft_exec, method = method)
    tailsA <- mafft(rev_DNA(seq_left), exec = mafft_exec, method = method)
  }

  ## Aligning right HEADS and TAILS
  if (length(seq_right) ==1){
    headsB <-  seq_right
    tailsB <- rev_DNA(seq_right)
  }else{
    headsB <- mafft(seq_right, exec = mafft_exec, method = method)
    tailsB <- mafft(rev_DNA(seq_right), exec = mafft_exec, method = method)
  }

  # aling 4 combinations of the basic MSAs (heads) and also
  # align in revers direction (tails)

  # headsA - headsB - HEADS
  msa1 <- mafft(x = headsA, y = headsB, add = "add",
    method= method, exec = mafft_exec)
  # headsA - headsB - TAILS
  msa2 <- mafft(x = rev_DNA(headsA), y = rev_DNA(headsB), add = "add",
    method= method, exec = mafft_exec)
  msa2 <- rev_DNA(msa2)

  # headsA - tailsB - HEADS
  msa3 <- mafft(headsA, rev_DNA(tailsB), add = "add",
    method= method, exec = mafft_exec)
  # headsA - tailsB - TAILS
  msa4 <- mafft(rev_DNA(headsA), tailsB, add = "add",  # rev(rev())=> normal
    method= method, exec = mafft_exec)
  msa4 <- rev_DNA(msa4)

  # tailsA - headsB - HEADS
  msa5 <- mafft(rev_DNA(tailsA), headsB, add = "add",
    method= method, exec = mafft_exec)
  # tailsA - headsB - TAILS
  msa6 <- mafft(tailsA, rev_DNA(headsB), add = "add",
    method= method, exec = mafft_exec)
  msa6 <- rev_DNA(msa6)

  # tailsA - tailsB - HEADS
  msa7 <- mafft(rev_DNA(tailsA), rev_DNA(tailsB), add = "add",
    method= method, exec = mafft_exec)
  msa8 <- mafft(tailsA, tailsB, add = "add",
    method= method, exec = mafft_exec)
  msa8 <- rev_DNA(msa8)

  comb_msas <- list(msa1, msa2, msa3, msa4,
    msa5, msa6, msa7, msa8)

  return(comb_msas)
}
