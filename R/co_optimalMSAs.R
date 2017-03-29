#' Co-optimal alignments from tree partitions
# Produce co-optimal MSAs from the N-3 bipartitions of the starting tree
# This function performs one iteration i
# 1: Take i partition and align both sides Heads (H) and Tails (T)
# and we get H1,H2 and T1, T2
# 2: align combinations again Heads and Tails: H1H2, H1T1, H1T2, T1T2 => 8 combinations
#' @export


align_part_set <- function(x, partition_set, exec, method, msa.program, type =type){

  seq_left <- x[partition_set>0]
  seq_right <- x[partition_set==0]

  ## MAFFT
  ###########
  if(msa.program == "mafft"){
  ## Aligning left HEADS and TAILS
  if (length(seq_left) == 1){
    headsA <-  seq_left
    tailsA <- rev_DNA(seq_left)
  }else{
    headsA <- mafft(seq_left, exec = exec, method = method)
    tailsA <- mafft(rev_DNA(seq_left), exec = exec, method = method)
  }
  ## Aligning right HEADS and TAILS
  if (length(seq_right) ==1){
    headsB <-  seq_right
    tailsB <- rev_DNA(seq_right)
  }else{
    headsB <- mafft(seq_right, exec = exec, method = method)
    tailsB <- mafft(rev_DNA(seq_right), exec = exec, method = method)
  }
  # aling 4 combinations of the basic MSAs (heads) and also
  # align in revers direction (tails)
  # headsA - headsB - HEADS
  msa1 <- mafft(x = headsA, y = headsB, add = "add",
    method= method, exec = exec)
  # headsA - headsB - TAILS
  msa2 <- mafft(x = rev_DNA(headsA), y = rev_DNA(headsB), add = "add",
    method= method, exec = exec)
  msa2 <- rev_DNA(msa2)
  # headsA - tailsB - HEADS
  msa3 <- mafft(headsA, rev_DNA(tailsB), add = "add",
    method= method, exec = exec)
  # headsA - tailsB - TAILS
  msa4 <- mafft(rev_DNA(headsA), tailsB, add = "add",  # rev(rev())=> normal
    method= method, exec = exec)
  msa4 <- rev_DNA(msa4)
  # tailsA - headsB - HEADS
  msa5 <- mafft(rev_DNA(tailsA), headsB, add = "add",
    method= method, exec = exec)
  # tailsA - headsB - TAILS
  msa6 <- mafft(tailsA, rev_DNA(headsB), add = "add",
    method= method, exec = exec)
  msa6 <- rev_DNA(msa6)
  # tailsA - tailsB - HEADS
  msa7 <- mafft(rev_DNA(tailsA), rev_DNA(tailsB), add = "add",
    method= method, exec = exec)
  msa8 <- mafft(tailsA, tailsB, add = "add",
    method= method, exec = exec)
  msa8 <- rev_DNA(msa8)
  ## list of 8 combs
  comb_msas <- list(msa1, msa2, msa3, msa4,
    msa5, msa6, msa7, msa8)
  }

  ## MUSCLE
  ###########
  if(msa.program == "muscle"){
    if (length(seq_left) == 1){
      headsA <-  seq_left
      tailsA <- rev_DNA(seq_left)
    }else{
      headsA <- muscle2(seq_left, exec = exec, type = type)
      tailsA <- muscle2(rev_DNA(seq_left), exec = exec, type = type)
    }
    if (length(seq_right) ==1){
      headsB <-  seq_right
      tailsB <- rev_DNA(seq_right)
    }else{
      headsB <- muscle2(seq_right, exec = exec, type = type)
      tailsB <- muscle2(rev_DNA(seq_right), exec = exec, type = type)
    }
    msa1 <- muscle2(x = headsA, y = headsB, exec = exec, type =type)
    msa2 <- muscle2(x = rev_DNA(headsA), y = rev_DNA(headsB),
      exec = exec, type =type)
    msa2 <- rev_DNA(msa2)
    msa3 <- muscle2(headsA, rev_DNA(tailsB),  exec = exec, type =type)
    msa4 <- muscle2(rev_DNA(headsA), tailsB,  exec = exec, type =type)
    msa4 <- rev_DNA(msa4)
    msa5 <- muscle2(rev_DNA(tailsA), headsB, exec = exec, type =type)
    msa6 <- muscle2(tailsA, rev_DNA(headsB), exec = exec, type =type)
    msa6 <- rev_DNA(msa6)
    msa7 <- muscle2(rev_DNA(tailsA), rev_DNA(tailsB), exec = exec, type =type)
    msa8 <- muscle2(tailsA, tailsB, exec = exec, type =type)
    msa8 <- rev_DNA(msa8)
    comb_msas <- list(msa1, msa2, msa3, msa4,
      msa5, msa6, msa7, msa8)
  }


  ## CLUSTAL
  ###########
  if(msa.program == "clustalo"){
    if (length(seq_left) == 1){
      headsA <-  seq_left
      tailsA <- rev_DNA(seq_left)
    }else{
      headsA <- clustalo(seq_left, exec = exec, type =type)
      tailsA <- clustalo(rev_DNA(seq_left), exec = exec, type =type)
    }
    if (length(seq_right) ==1){
      headsB <-  seq_right
      tailsB <- rev_DNA(seq_right)
    }else{
      headsB <- clustalo(seq_right,  exec = exec, type =type)
      tailsB <- clustalo(rev_DNA(seq_right),exec = exec, type =type)
    }
    msa1 <- clustalo(x = headsA, y = headsB,exec = exec, type =type)
    msa2 <- clustalo(x = rev_DNA(headsA), y = rev_DNA(headsB),
      exec = exec, type =type)
    msa2 <- rev_DNA(msa2)
    msa3 <- clustalo(headsA, rev_DNA(tailsB), exec = exec, type =type)
    msa4 <- clustalo(rev_DNA(headsA), tailsB, exec = exec, type =type)
    msa4 <- rev_DNA(msa4)
    msa5 <- clustalo(rev_DNA(tailsA), headsB, exec = exec, type =type)
    msa6 <- clustalo(tailsA, rev_DNA(headsB), exec = exec, type =type)
    msa6 <- rev_DNA(msa6)
    msa7 <- clustalo(rev_DNA(tailsA), rev_DNA(tailsB), exec = exec, type =type)
    msa8 <- clustalo(tailsA, tailsB,  exec = exec, type =type)
    msa8 <- rev_DNA(msa8)
    comb_msas <- list(msa1, msa2, msa3, msa4,
      msa5, msa6, msa7, msa8)
  }


  ## CLUSTALW2
  #############
  if(msa.program == "clustalw2"){
    if (length(seq_left) == 1){
      headsA <-  seq_left
      tailsA <- rev_DNA(seq_left)
    }else{
      headsA <- clustalw2(seq_left, exec = exec, type =type)
      tailsA <- clustalw2(rev_DNA(seq_left), exec = exec, type =type)
    }
    if (length(seq_right) ==1){
      headsB <-  seq_right
      tailsB <- rev_DNA(seq_right)
    }else{
      headsB <- clustalw2(seq_right,  exec = exec, type =type)
      tailsB <- clustalw2(rev_DNA(seq_right),exec = exec, type =type)
    }
    msa1 <- clustalw2(x = headsA, y = headsB,exec = exec, type =type)
    msa2 <- clustalw2(x = rev_DNA(headsA), y = rev_DNA(headsB),
      exec = exec, type =type)
    msa2 <- rev_DNA(msa2)
    msa3 <- clustalw2(headsA, rev_DNA(tailsB), exec = exec, type =type)
    msa4 <- clustalw2(rev_DNA(headsA), tailsB, exec = exec, type =type)
    msa4 <- rev_DNA(msa4)
    msa5 <- clustalw2(rev_DNA(tailsA), headsB, exec = exec, type =type)
    msa6 <- clustalw2(tailsA, rev_DNA(headsB), exec = exec, type =type)
    msa6 <- rev_DNA(msa6)
    msa7 <- clustalw2(rev_DNA(tailsA), rev_DNA(tailsB), exec = exec, type =type)
    msa8 <- clustalw2(tailsA, tailsB,  exec = exec, type =type)
    msa8 <- rev_DNA(msa8)
    comb_msas <- list(msa1, msa2, msa3, msa4,
      msa5, msa6, msa7, msa8)
  }
  return(comb_msas)
}
