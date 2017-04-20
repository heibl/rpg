## This code is part of the rpg package
## © C. Heibl 2016 (last update 2017-04-12)

#' @title Ultra-Large Multiple Sequence Alignment with PASTA
#' @description Provides a complete reimplementation of the PASTA algorithm (Mirarab, Nguyen, and Warnow 2014) in R.
#' @param seqs An object of class \code{\link{DNAbin}} or \code{\link{AAbin}}
#'   containing unaligned sequences of DNA or amino acids.
#' @param gt \emph{Currently unused.}
#' @param k An integer giving the size of cluster in which the dataset is split.
#' @importFrom igraph as_edgelist
#' @importFrom ips mafft mafft.merge
#' @export

pasta <- function(seqs, gt, k = 200, cutoff = 0.93, parallel = FALSE,
                  bootstrap = 100, msa.program = "mafft", method = "auto",
                  exec, ncore){

  ## remove gaps from aligned sequences
  ## ----------------------------------
  if (is.matrix(seqs)) {
    seqs <- del.gaps(seqs)
  }

  ## less than k species will be aligned with MAFFT-LINSI
  ## ----------------------------------------------------
  if (length(seqs) <= k){
    cat(length(seqs), "species will be aligned with MAFFT L-INS-i\n")

    seqs <- guidance(seqs, parallel = parallel, ncore = ncore,
                     bootstrap = bootstrap, msa.program = msa.program,
                     method = method, exec = exec)

    ## more than k species will be aligned with PASTA
    ## ----------------------------------------------
  } else {

    ## This is a quick hack to get an inital guide tree
    ## Should be replaced by the method used by Mirarab and Warnow
    ## or perhaps a hybrid with taxonomy.
    if (missing(gt)){
      gt <- mafft(seqs, method = "auto")
      if (inherits(gt, "DNAbin")){
        gt <- dist.dna(gt, model = "F81")
      } else {
        gt <- dist.aa(gt)
      }
      gt <- nj(gt)
    }

    ## split dataset in subsets of size <= k
    ## -------------------------------------
    subtrees <- centroidDecomposition(gt, k = k)
    subtrees <- lapply(subtrees, function(z) z$tip.label)
    names(subtrees) <- paste0("S", seq_along(subtrees))

    ## alignment of subtrees
    ## ---------------------
    foo <- function(seqs, taxa){
      guidance(seqs[taxa], parallel = parallel, ncore = ncore,
               bootstrap = bootstrap, msa.program = msa.program,
               method = method, exec = exec)
    }
    s <- lapply(subtrees, foo, seqs = seqs)
    names(s) <- names(subtrees)
    seqs <- lapply(s, function(z) z$base_msa)

    ## compute spanning tree of subsets
    ## --------------------------------
    st <- spanningTree(gt, subtrees)

    ## do profile-alignment
    ## --------------------
    e <- as_edgelist(st)
    merger <- function(seqlist, index, exec){
      mafft.merge(seqlist[index], exec = exec)
    }
    seqs <- apply(e, 1, merger, seqlist = seqs, exec = exec)
    names(seqs) <- paste(e[, 1], e[, 2], sep = "-")

    save.image("devworkspace.rda")

    ## do transitivity merging
    ## -----------------------
    load("devworkspace.rda")


  }
  seqs
}
