#' Reverse DNA sequences from class DNAbin or AAbin
#' @export
#'
rev_DNA <- function(x){

  if (!inherits(x, c("DNAbin","AAbin")))
    stop("sequences not of class DNAbin or AAbin (ape)")

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
