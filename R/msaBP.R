#' Bootstrap a multiple sequence alignnment
#' @import parallel
#' @import doSNOW
#' @export

msaBP <- function(msa){

  if (!inherits(msa, c("DNAbin", "AAbin")))
    stop("msa not of class DNAbin or AAbin (ape)")

  msa <- msa[,sample(ncol(msa), replace = TRUE)]

  return(msa)
}

