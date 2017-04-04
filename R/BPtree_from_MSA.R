#' Generate NJ BP tree from MSA
#'
#' @param msa multiple sequence alignment of class AAbin or DNAbin
#' @return tree of class phylo
#' @import phangorn
#' @import ape
#' @export


msaBP_nj_tree <- function(msa, outgroup){

  base.msa.bp <- msaBP(msa)
  # convert to class phyDAT
  base.msa.ml <- as.phyDat(as.character(base.msa.bp))
  # find ML distance as input to nj tree search
  ml.dist.msa <- dist.ml(base.msa.ml)
  # NJ
  tr <- ape::nj(ml.dist.msa)
  # root
  if(outgroup == "auto")
    outgroup <- tr$tip.label[1]
  tr <- root(tr, outgroup = outgroup)
  # resolve polytomies
  tr <- multi2di(tr)
  ## Rescale branch lengths
  tr <- compute.brlen(tr)

  return(tr)
}
