#' Compare reference MSAs with alternative MSAs
#'
#' @param ref of class data.frame, is the reference MSA ('BASE MSA') with sequences as columns
#' @param com like ref, but 1 alternative MSA
#' @param dir_path directory with multiple alternative MSA files
#'
#' @return list containing following scores:
#' @return mean_scores residue pair score and mean column score
#' @return column_score
#' @return residue_column_score
#' @return residue_pair_residue_score
#' @return residual_pair_sequence_pair_score
#' @return residual_pair_sequence_score
#' @return residue_pair_score
#'
#' @description Wrapper function for program msa_set_score v2.01 of the GUIDANCE program (see reference). Copyright: To modify the code, or use parts of it for other purposes, permission should be requested. Please contact Tal Pupko: talp@post.tau.ac.il. Please note that the use of the GUIDANCE program is for academic use only
#' @references Penn et al. (2010). An alignment confidence score capturing
#'   robustness to guide tree uncertainty. Molecular Biology and Evolution
#'   27:1759--1767
#'
#' @author Franz-Sebastian Krah
#' @import plyr
#' @export

compareMSAs <- function(ref, com, dir_path){

  if (!inherits(ref, c("AAbin", "DNAbin")))
    stop("ref not of class AAbin or DNAbin")

  # create temporary dir with temporary fasta files
  fns <- vector(length = 2)
  rn <- format(runif(1,1,100000), digits = 0, scientific = FALSE)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = paste("run", rn, sep="")
      , tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  if(is.null(com)){

    if(!dir.exists(dir_path))
      stop("Dir not found")

    write.fas(ref, fns[1])
    system(paste("src/msa_set_score_src/msa_set_score",
      fns[1],
      paste(tempdir(), "GUIDANCE", sep ="/"),
      "-d",
      dir_path), intern = TRUE, ignore.stdout = FALSE)

    ## read program putput which is in temp dir
    files <-  list.files(tempdir(), full.names = TRUE)
    read <- files[grep("GUIDANCE", files)]


  }else{

    if (!inherits(com, c("AAbin", "DNAbin")))
      stop("ref not of class AAbin or DNAbin")

    # store fasta files in temp files
    # rownames(ref) <- paste0("S", 1:dim(ref)[1])
    # rownames(ref) <- paste0("S", 1:dim(ref)[1])

    write.fas(ref, file = fns[1])
    write.fas(com, file = fns[2])

    ## call program msa_set_score in R package folder src/msa_set_score_src
    system(paste(paste("src/msa_set_score_src/msa_set_score", sep ="/"),
      fns[1],
      paste(tempdir(), rn, sep ="/"),
      "-m",
      fns[2]), intern = TRUE, ignore.stdout = FALSE)


    ## read program putput which is in temp dir
    files <-  list.files(tempdir(), full.names = TRUE)
    files <- files[grep(rn, files)]
    read <- files[-grep("\\.fas", files)]
  }

  ## Scores
  # Column score
  CS <- read.table(read[1])
  names(CS) <- c("col", "CS")

  # Mean scores
  msc <- readLines(read[2])
  msc <- gsub("#", "", msc[5])
  msc <- strsplit(msc, "  ")
  msc <- unlist(msc)
  msc <- data.frame(do.call(cbind, strsplit(msc, " ")))
  msc[] <- lapply(msc, as.character)
  names(msc) <- c(msc[1,1], msc[1,2])
  msc <- msc[-1,]


  # Residue pair column score (GUIDANCE Score)
  g.sc <- read.table(read[3])
  names(g.sc) <- c("col", "col_score")

  # Residue pair residue score
  rpr.sc <- read.table(read[4])
  names(rpr.sc) <- c("col", "residue", "score")

  # Residual pair sequence pair score
  rpsp.sc <- read.table(read[5])
  names(rpsp.sc) <- c("seq_row1", "seq_row2", "score")

  # Residual pair sequence score
  rps.sc <- read.table(read[6])
  names(rps.sc) <- c("seq", "score")

  # Residue pair score

  rp.sc <- read.table(read[7])
  names(rp.sc) <- c("col1", "row1", "row2", "score")

  ## delete temp files
  unlink(fns[file.exists(fns)], force = TRUE)
  unlink(fns, force = TRUE)
  files <- list.files(tempdir(), full.names = TRUE)
  files <- files[-grep("rs-graphics", files)]
  unlink(files, force = TRUE, recursive = TRUE)


  # output
  res <- list(score_means = msc,
    column_score = CS,
    residue_pair_column_score = g.sc,
    residue_pair_residue_score = rpr.sc,
    residual_pair_sequence_pair_score  = rpsp.sc,
    residual_pair_sequence_score = rps.sc,
    residue_pair_score = rp.sc)

  return(res)
}

