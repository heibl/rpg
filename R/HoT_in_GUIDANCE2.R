#' @export

Hot_GUIDANCE2 <- function(msa, n.coopt, type, td,
  files, raw_seq, msa.program, method,
  exec){

  # create start_tree
  li <- msa
  msa <- read.fas(msa, type = type)
  msa.nj <- ape::nj(dist.ml(as.phyDat(msa)))
  msa.nj <- root(msa.nj, outgroup = msa.nj$tip.label[1])
  msa.nj <- multi2di(msa.nj)
  msa.nj <- compute.brlen(msa.nj)

  ## produce MSA partitions
  align_parts <- partitions(msa.nj)

  ## sample 4 or n co-optimal
  nt <- Ntip(msa.nj)
  n.co <- sample((nt-3)*8, n.coopt)
  n.co_mat <- data.frame(n = 1:((nt-3)*8),
    part = rep(1:(nt-3), each = 8),
    n.in.part = rep(1:8, times = (nt-3)))
  n.co.sub <- n.co_mat[n.co_mat$n %in% n.co,]

  # reduce partitions to the randomly choosen co-optimal MSA number
  align_parts <- align_parts[,n.co.sub$part]
  # number of random MSA within partition (remember, each partition has 8 alignments)
  n.co.sub <- n.co.sub$n.in.part

  # make the 4 or n alignments
  alt_msas <- foreach(z = 1:ncol(align_parts)) %do% {
    # setTxtProgressBar(pb, i)
    align_part_set(x = raw_seq, partition_set = align_parts[,z],
      method = method, exec = exec, msa.program = msa.program,
      coopt.sub = n.co.sub[z])
  }

  ## unlist
  alt_msas <- foreach(k = 1:length(alt_msas), .combine = c) %do% {
    alt_msas[[k]]
  }

  ## write to files
  for(j in 1:n.coopt)
    write.fas(alt_msas[[j]], files[j])
}


