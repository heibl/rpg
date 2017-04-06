#' SEMPHY
#'
#' @param msa multiple sequence alignment of \code{class} DNAbin or AAbin
#' @param dist.table Distance Table Estimation Method
#' @param model Evolutionary model (e.g. "hky")
#' @param asrv Among Site Rate Variation
#' @param categories Rate categories
#' @param exec dir where semphy is installed
#' @param bootstrap number of bootstraps, if should be performed
#' @param verbose Verbosity level of the log file (between 0 and 10, default = 8)
#' @param MoreArgs further parameters parsed to semphy
#'
#' @references Friedman et al. (2002). A Structural EM Algorithm for Phylogenetic Inference Journal of Computational Biology. 9(2):331-53
#' @references Ninio* and Privman* (2007). Phylogeny Reconstruction: Increasing the Accuracy of Pairwise Distance Estimation Using Bayesian Inference of Evolutionary Rates. Bioinformatics. 23: e136-e141
#' @references http://compbio.cs.huji.ac.il/semphy/
#'
#' @return likelihood
#' @return Ml_tree maximum likelihood tree
#' @return bootstrap_trees if wanted
#'
#' @import ips
#' @import stringr


# seq <- read.fas("dev/data/AATF.fas", type ="AA")
# msa <- mafft(seq)
# ML <- "-H"
# dist.table<- "-J"
# model <- "--aaJC"
# arsv <- "-H"
# bootstrap = 100

semphy(msa = msa,
  dist.table<- "-J",
  model <- "--aaJC",
  arsv <- "-H",
  bootstrap = 100)

semphy <- function(msa, dist.table, model, arsv,
  categories, exec= "/Applications/semphy-2.0b2-bin-mac-intel",
  bootstrap, verbose, MoreArgs){

  if (!inherits(msa, c("DNAbin", "AAbin")))
    stop("MSA not of class DNAbin or AAbin (ape)")

  out <- system(paste(exec, "--v", sep=" "), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (out == 127)
    stop("please provide exec path or install MSA program in root \n
      i.e. in mac: '/Applications/semphy-2.0b2-bin-mac-intel'")


  fns <- vector(length = 4)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "semphy", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  write.fas(msa, fns[1])

  type <- class(msa)
  type <- gsub("bin", "", type)
  if(type == "DNA")
    a <- 4
  if(type == "AA")
    a <- 20


  if(missing(verbose))
    verbose <- "-v 8"

  if(missing(dist.table)) stop("dist.table must be specified, one of:\n posteriorDTME (iterative NJ) \n homogeneousRatesDTME (homogeneous rates model) \n rate4siteDTME (ML rate for each site)")
  if(dist.table == "auto")
    dist.table <- "--posteriorDTME"


  if(missing(model)) stop("model must be specified, one of:\n JC for nucleotides (--nucjc) or amino acids (--aaJC) \n K2P (--k2p) \n HKY (--hky) \n Dayhoff (--day) \n JTT (--jtt) \n REV (--rev) \n WAG (--wag) \n cpREV (--cprev)")
  if(model == "auto")
    if(type == 4){
      model <- "--hky"
    }else{ model <- "--jtt" }

  if(missing(arsv)) stop("Among rate site variation must be specified, one of: -H (homogeneous rates), -A (Set the initial alpha parameter), or -O (Optimize the alpha parameter)")
  if(arsv == "auto")
    arsv <- "-O"

  if(missing(categories))
    categories <- ""

  if(missing(MoreArgs))
    MoreArgs <- ""

  semphy.call <- paste(paste(program.dir, "semphy", sep="/"),
    "-s", fns[1], "-o", fns[2], "-T", fns[3],
    "-l", fns[4],"-a", a, model, verbose, dist.table, arsv, MoreArgs, sep=" ")

  if(!missing(bootstrap)){
    bp <- paste("--BPrepeats", bootstrap, sep=" ")
    semphy.call <- paste(semphy.call, bp)
  }

  # CALL
  system(semphy.call)

  # read lik
  lik <- readLines(fns[2], warn = FALSE)
  likstart <- grep("likelihood", lik)+1
  lik <- lik[likstart]

  # read ml tree
  tr <- read.tree(fns[3])

  # read bootstrap trees
  if(!missing(bootstrap)){
    tbp <- read.tree(fns[4])
    tbp <- tbp[-c(1:3, length(tbp))]
    tbp <- tbp[seq(1,2*bootstrap, 2)]


    bp <- readLines(fns[2], warn = FALSE)
    bpstart <- grep("bootstrap", bp)+1
    bp <- bp[bpstart]
    library(stringr)
    bp <- str_extract_all(bp[[1]], "\\[[0-9\\.0-9*]*\\]")
    bp <- as.numeric(gsub("\\[|\\]", "", bp[[1]]))
    tr$node.label <- bp
  }

  if(!missing(bootstrap)){
    res <- list(log.likelihood = lik, ML_tree = tr, bootstrap_trees = tbp)
  }else{ res <- tr }

  files <- list.files(tempdir(), full.names = TRUE)
  files <- files[-grep("rs-graphics", files)]
  unlink(files, force = TRUE, recursive = TRUE)

  return(res)

}
