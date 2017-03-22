#' Wrapper Function for Calculation of GUIDANCE Scores
#'
#' @param ref of class data.frame, is the reference MSA ('BASE MSA') with sequences as columns
#' @param com like ref, but 1 alternative ('perturbated MSA')
#' @param n_id is used as 'i' if function is run in parallel
#' @return list containing GUIDANCE column score, the residue score and sequence score
#'
#' @author Franz-Sebastian Krah
#' @export

calc_scores <- function(ref, com){
  # create temporary dir with temporary fasta files
  fns <- vector(length = 2)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "mafft", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  # store fasta files in temp files
  seqinr::write.fasta(as.list(ref), names = paste0("S", 1:dim(ref)[1]),
    file = fns[1])
  seqinr::write.fasta(as.list(com), names = paste0("S", 1:dim(com)[1]),
    file = fns[2])

  ## call program msa_set_score in R package folder src/msa_set_score_src
  system(paste(paste("src/msa_set_score_src/msa_set_score", sep ="/"),
    fns[1],
    paste(tempdir(), "score_result", sep ="/"),
    "-m",
    fns[2]), intern = TRUE, ignore.stdout = FALSE)


  ## read program putput which is in temp dir
  files <-  list.files(tempdir())
  read <- files[grep("score", files)]
  # GUIDANCE Score
  gsc <- read.table(paste(tempdir(), read[3], sep ="/"), header = FALSE)
  # G residue score
  grsc <- read.table(paste(tempdir(), read[4], sep ="/"), header = FALSE)
  # sequence score
  ssc <- read.table(paste(tempdir(), read[6], sep ="/"), header = FALSE)

  ## delete temp files
  unlink(fns[file.exists(fns)])
  unlink(tempdir())

  # output
  res <- list(GUIDANCE_score = gsc, residue_score = grsc, sequence_score = ssc)
  return(res)
}

