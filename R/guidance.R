#################################################
##### GUIDANCE
##### Franz-S. Krah
##### 02/26/2017
##### after Penn et al. 2010


# require("stringr")
# require("doSNOW")
# # require("bigmemory")
# # require("biganalytics")
# # require("bigtabulate")
# require("Rcpp")
# require("pbmcapply")
# require("plyr")
# require("ips")
# require("AlignStat")
#
#
#
# source("R/functions.R")
#
#

#   require("stringr")
#   # require("bigmemory")
#   # require("biganalytics")
#   # require("bigtabulate")
#   require("pbmcapply")
#   require("AlignStat")

#' @import doSNOW
#' @import foreach
#' @import parallel
#' @import pbmcapply

guidance <- function(dna, G.CS.th = 0.93, parallel = FALSE, ncore,
  bootstrap, msa.programm = "mafft", method = "auto"){

  if (!is.object(dna)){
    dna <- read.FASTA(dna)
  }

  ##############################################
  ## SOME CHECKS
  ##############################################
  if (!inherits(dna, "DNAbin"))
    stop("sequences not of class DNAbin (ape)")

  ##############################################
  ## PART I
  ##############################################
  ## BASE and PERTUBATED MSAs
  ##############################################


  ## Generate BASE alignment
  ###########################
  cat("Generating the base alignment \n")
  if (msa.programm == "mafft"){
    base.msa <- ips::mafft(dna, method = method)
  }
  base.msa <- as.character(base.msa)

  ## Constructing BP guide-trees for the pertubated MSAs
  #######################################################
  cat("Pertubating base alignment \n")
  pb <- txtProgressBar(min = 0, max = bootstrap, style = 3)

  base.msa.bp <- foreach(i = 1:bootstrap) %do% {
    setTxtProgressBar(pb, i)
    pertubatedMSA <- as.DNAbin(base.msa[,sample(ncol(base.msa), replace = TRUE)])
    return(pertubatedMSA)
  }
  close(pb)

  ## Generating alternative (pertubated) MSAs
  ###########################################
  cat("Generating alternative alignments \n")

  ## Compute NJ guide trees
  cat("  Generating NJ guide trees \n")
  if (parallel == TRUE){
    pb <- txtProgressBar(max = bootstrap, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    cl <- makeCluster(ncore)
    registerDoSNOW(cl)
    nj.guide.trees <- foreach(i = 1:bootstrap, .options.snow = opts) %dopar% {
      # convert to class phyDAT
      base.msa.ml <- phangorn::as.phyDat(base.msa.bp[[i]])
      # find ML distance as input to nj tree search
      ml.dist.msa <- phangorn::dist.ml(base.msa.ml)
      # NJ
      ape::nj(ml.dist.msa)
    }
    stopCluster(cl)
  }
  if (parallel == FALSE){
    nj.guide.trees <- foreach(i = 1:bootstrap) %do% {
      setTxtProgressBar(pb, i)
      # convert to class phyDAT
      base.msa.ml <- phangorn::as.phyDat(base.msa.bp[[i]])
      # find ML distance as input to nj tree search
      ml.dist.msa <- phangorn::dist.ml(base.msa.ml)
      # NJ
      ape::nj(ml.dist.msa)
    }
  }
  close(pb)


  ## Find NJ tree for base MSA
  ## 1st tip label of base.nj tree as outgroup for all guide.njs
  base.nj <- ape::nj(phangorn::dist.ml(phangorn::as.phyDat(base.msa)))

  ## Root each tree on the first tip label of the base.nj tree
  nj.guide.trees <- lapply(nj.guide.trees, root,
                           outgroup = base.nj$tip.label[1],
                           resolve.root = TRUE)

  ## Rescale branch lengths
  nj.guide.trees <- lapply(nj.guide.trees, compute.brlen)

  ## Alignment of MSA BP times with new NJ guide trees
  cat("  Alignment of pertubated MSAs using NJ guide trees \n")

  if (parallel == TRUE){
    # guide.msa <- pbmclapply(nj.guide.trees,
    #                         function(x, m, gt) ips::mafft(x, method, gt),
    #                         x = dna, m = method,
    #                         mc.cores = ncore, ignore.interactive = TRUE)
    guide.msa <- pbmclapply(nj.guide.trees, FUN = function(x) ips::mafft(dna, gt = x, method = "retree 1"), mc.cores = ncore, ignore.interactive = TRUE)
  }

  if (parallel == FALSE){

    pb <- txtProgressBar(min = 0, max = bootstrap, style = 3)

    guide.msa <- foreach(i = 1:bootstrap) %do% {
      ips::mafft(x = dna, gt = nj.guide.trees[[i]])
      setTxtProgressBar(pb, i)
    }
  }
  close(pb)

  unlist(lapply(guide.msa, function(x) dim(x)[2])) == dim(base.msa)[2]

  if(any(unlist(lapply(lapply(guide.msa, dim), is.null))) == TRUE)
    cat("one or more pertubated MSAs were NULL \n ")

  guide.msa <- guide.msa[!unlist(lapply(lapply(guide.msa, dim), is.null))]
  bootstrap <- length(guide.msa)

  mafft_created <- list.files(getwd(),
    full.names = T)[grep("tree.mafft", list.files(getwd()))]
  if(length(mafft_created)>0){
    file.remove(mafft_created)
  }
  rm(mafft_created)


  # remooving some intermediate objects not further needed
  rm(base.msa.bp)
  rm(nj.guide.trees)

  # save(base.msa, file ="out/base.msa_200.rda")
  # save(guide.msa, file ="out/guide.msa_200.rda")

  #
  load(file ="dev/GUIDANCE_in_R/out/guide.msa_200.rda")
  load(file ="dev/GUIDANCE_in_R/out/base.msa_200.rda")

  ##############################################
  ## PART II
  ##############################################
  ## Computation of GUIDANCE scores
  ##############################################
  cat("Calculating GUIDANCE scores \n")

  load(file ="dev/GUIDANCE_in_R/out/guide.msa_10.rda")
  load(file ="dev/GUIDANCE_in_R/out/base.msa_10.rda")

  guide.msa <- lapply(guide.msa, as.character)

  # ## GUIDANCE Score
  # ##############################################
  # # produce combinations factor matrix


 ####
 # HERE ADD SCORE CODE
 ####
  ## Produce output
  res <-  list(alignment_score = alignment_score,
    # GUIDANCE_residue_pair_score = rpsc,
    GUIDANCE_score = GUIDANCE_score,
    guidance_msa =guidance.msa,
    base_msa = base.msa)

  ## Return output
  return(res)

}
