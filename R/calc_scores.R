#' Wrapper Function of C Code for Calculation of the GUIDANCE Scores
#'
#' @param ref of class data.frame, is the reference MSA ('BASE MSA') with sequences as columns
#' @param com like ref, but 1 alternative ('perturbated MSA')
#' @param n_id is used as 'i' if function is run in parallel
#' @return list containing GUIDANCE column score, the residue score and sequence score
#'
#' @author Franz-Sebastian Krah

calc_scores <- function(ref, com, n_id = 1){

  locdir <- getwd()
  locdir_runs <- paste0(paste(locdir, "score_run", sep ="/"), n_id)
  dir.create(locdir_runs)
  seqinr::write.fasta(as.list(ref), names = paste0("S", 1:dim(ref)[1]),
    file =paste(locdir_runs, "ref.fasta", sep = "/"))
  seqinr::write.fasta(as.list(com), names = paste0("S", 1:dim(com)[1]),
    file = paste(locdir_runs, "com.fasta", sep = "/"))

  system(paste(paste(locdir, "src/msa_set_score", sep ="/"),
    paste(locdir_runs, "ref.fasta", sep="/"),
    paste(paste("score_run", n_id, sep =""), "score_res", sep="/"),
    "-m",
    paste(paste("score_run", n_id, sep =""), "com.fasta", sep="/")))

  files <- list.files(locdir_runs)
  read <- files[grep("score", files)]
  # GUIDANCE Score
  gsc <- read.table(paste(locdir_runs, read[3], sep ="/"), header = FALSE)
  # G residue score
  grsc <- read.table(paste(locdir_runs, read[4], sep ="/"), header = FALSE)
  # sequence score
  ssc <- read.table(paste(locdir_runs, read[6], sep ="/"), header = FALSE)

  unlink(locdir_runs, recursive = TRUE, force = TRUE)

  res <- list(GUIDANCE_score = gsc, residue_score = grsc, sequence_score = ssc)
  return(res)
}
