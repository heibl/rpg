## This code is part of the rpg package
## © C. Heibl 2017 (last update 2017-04-06)

#' @title Minimum-Spanning-Tree for PASTA
#' @description Creates miminum spanning tree connecting subsets obtained by
#'   centroid decomposition.
#' @param phy An object of class \code{\link{phylo}} containing a complete phylogenetic tree.
#' @param subtrees A list containing object of class \code{\link{phylo}} as
#'   obtained with \code{\link{centroidDecomposition}}.
#' @importFrom ape drop.tip
#' @importFrom ips terminalSisters
#' @import igraph
#' @export

spanningTree <- function(phy, subtrees){

  # subtrees <- lapply(subtrees, function(z) z$tip.label)
  # names(subtrees) <- paste0("S", seq_along(subtrees))

  ## replace tiplabels by their subtree name
  ## ---------------------------------------
  for (i in names(subtrees)){
    phy$tip.label[phy$tip.label %in% subtrees[[i]]] <- i
  }

  ## shrink monophyletic subtrees
  ## ----------------------------
  repeat {
    ts <- terminalSisters(phy, labels = FALSE)
    id <- sapply(ts, function(z, phy) phy$tip.label[z[1]] == phy$tip.label[z[2]], phy = phy)
    if (!any(id)) break
    id <- sapply(ts[id], head, n = 1)
    phy <- drop.tip(phy, id)
  }
  id <- lapply(names(subtrees), function(z, phy) which(phy$tip.label == z), phy = phy)
  id <- lapply(id, tail, n = -1)
  phy <- drop.tip(phy, unlist(id))

  ## convert to igraph and do propagating of tip labels
  ## to internal nodes
  ## -----------------
  phy <- as.igraph.phylo(phy, directed = FALSE)
  repeat {
    vertex_names <- vertex_attr(phy)$name
    int <- sample(grep("Node", vertex_names), 1)
    n <- vertex_names[neighbors(phy, int)]
    id <- grep("Node", n)
    if (length(id)) n <- n[-id]
    if (!length(n)) next
    if (length(n)) n <- sample(n, 1)
    deg <- degree(phy, n)
    ## handle nodes of degree 1
    if (deg == 1){
      phy <- delete.vertices(phy, n)
      vertex_attr(phy)$name[match(vertex_names[int], vertex_attr(phy)$name)] <- n
    } else {
        nn <- vertex_attr(phy)$name[neighbors(phy, n)]
        phy <- delete.vertices(phy, n)
        phy <- add_edges(phy, nn)
        vertex_attr(phy)$name[match(vertex_names[int], vertex_attr(phy)$name)] <- n
    }
    ## stopping conditition
    vertex_names <- vertex_attr(phy)$name
    int <- grep("Node", vertex_names)
    if (!length(grep("Node", vertex_attr(phy)$name))) break
  }
  phy
}

