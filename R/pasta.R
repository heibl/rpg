## This code is part of the rpg package
## © C. Heibl 2016 (last update 2017-04-06)

#' @title Ultra-Large Multiple Sequence Alignment with PASTA
#' @description Provides a complete reimplementation of the PASTA algorithm (Mirarab, Nguyen, and Warnow 2014) in R.
#' @param seqs An object of class \code{\link{DNAbin}} or \code{\link{AAbin}}
#'   containing unaligned sequences of DNA or amino acids.
#' @param gt \emph{Currently unused.}
#' @param k An integer giving the size of cluster in which the dataset is split.
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

    seqs <- guidance(seqs, cutoff = cutoff, parallel = parallel, ncore = ncore,
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
    subtree <- centroidDecomposition(gt, k = k)
    subtree <- lapply(subtree, function(z) z$tip.label)
    names(subtree) <- paste0("S", seq_along(subtree))

    ## compute spanning tree of subsets
    ## --------------------------------
    st <- spanningTreeh(gt, subtree)

    ## do profile-alignment
    ## --------------------


  }
  seqs
}
