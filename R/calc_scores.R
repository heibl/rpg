#' Wrapper Function for Calculation of GUIDANCE Scores
#'
#' @param ref of class data.frame, is the reference MSA ('BASE MSA') with sequences as columns
#' @param com like ref, but 1 alternative ('perturbated MSA')
#' @param n_id is used as 'i' if function is run in parallel
#' @return list containing GUIDANCE column score, the residue score and sequence score
#'
#' @author Franz-Sebastian Krah
#' @import plyr
#' @export

calc_scores <- function(ref, com){
  # create temporary dir with temporary fasta files

  fns <- vector(length = 2)
  rn <- format(runif(1,1,100000), digits = 0, scientific = FALSE)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = paste("run", rn, sep="")
      , tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  # store fasta files in temp files
  seqinr::write.fasta(as.list(ref), names = paste0("S", 1:dim(ref)[1]),
    file = fns[1])
  seqinr::write.fasta(as.list(com), names = paste0("S", 1:dim(com)[1]),
    file = fns[2])

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
  read

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
  unlink(fns[file.exists(fns)])

  # output
  res <- list(mean_score = msc,
    column_score = CS,
    GUIDANCE_score = g.sc,
    Residue_pair_residue_score = rpr.sc,
    Residual_pair_sequence_pair_score  = rpsp.sc,
    Residual_pair_sequence_score = rps.sc,
    Residue_pair_score = rp.sc)
  return(res)
}

