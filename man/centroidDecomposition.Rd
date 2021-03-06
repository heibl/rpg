% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/centroidDecomposition.R
\name{centroidDecomposition}
\alias{centroidDecomposition}
\title{Centroid Decomposition of Phylogenetic Trees}
\usage{
centroidDecomposition(phy, k = 200)
}
\arguments{
\item{phy}{An object of class \code{\link{phylo}}.}

\item{k}{An integer giving the upper size limit of subtrees.}
}
\description{
Splits a phylogenetic tree into a number of subtree
with size at most k.
}
