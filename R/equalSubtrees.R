## This code is part of the rpg package
## © C. Heibl 2016 (last update 2017-03-15)

#' @title Split Phylogeny Equally
#' @description Split a phylogenetic tree into two equally sized subtrees.
#' @param phy A phylogenetic tree of class \code{\link{phylo}}.
#' @param k An integer, only trees with > k tips will be split.
#' @return A list containing two subtrees of class \code{\link{"phylo}}.
#' @details The phylogenetic tree is treated as unrooted, no matter whether is
#'   rooted or not.
#' @seealso \code{\link{decomposePhylo}} to decompose large phylogenetic trees
#'   into subsets of a certain size.
#' @importFrom ips descendants
#' @export

equalSubtrees <- function(phy, k){

  if (!missing(k)){
    if (Ntip(phy) <= k) return(phy)
  }

  ## extract edge matrix (it will be modified)
  e <- phy$edge # edge matrix

  ## delete terminal edges
  e <- e[!e[, 2] %in% 1:Ntip(phy), ]

  ## calculate number of tips of second-column nodes
  ## (tip number of first-column nodes is simply Ntip - nn)
  n <- lapply(e[, 2], descendants, phy = phy)
  nn <- sapply(n, length)

  ## identigy most balances split
  id <- which.min(abs(nn/Ntip(phy) - 0.5))
  # cat("most balanced split at edge", id, ":", nn[id], "versus", Ntip(phy) - nn[id], "tips")

  ## split into subclades
  list(extract.clade(phy, e[id, 2]),
       drop.tip(phy, n[[id]]))
}
