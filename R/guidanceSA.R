#' @import stringr
#' @import useful
#' @import foreach
#' @import ips
#' @export
#'

guidanceSA <- function(sequences, msa.program, programm,
  bootstrap, gencode, outorder, msafile, cutoff= 0.93,
  moreArgs, guidance_dir,proc_num, quiet = FALSE){

  if (!inherits(sequences, c("DNAbin", "AAbin")))
    stop("sequences not of class DNAbin or AAbin (ape)")

  type <- class(sequences)
  type <- gsub("bin", "", type)
  if(type == "DNA")
    type <- "nuc"
  if(type == "AA")
    type <- "aa"

  if(missing(guidance_dir))
    guidance_dir <- "/Applications/guidance.v2.02/"

  fns <- vector(length = 1)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "guidance", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  write.fas(sequences, fns[1])

  # necessary
  seqFile <- paste("--seqFile", fns[1])
  if(lower.case(msa.program))
    msa.program <- str_to_upper(msa.program)
  msa.program <- paste("--msaProgram", msa.program, sep=" ")

  seqType <- paste("--seqType", type, sep=" ")
  # outDir <- paste("--outDir", outdir, sep=" ")
  if(programm == "hot")
    programm <- "HoT"
  if(lower.case(programm))
    programm <- str_to_upper(programm)
  programm <- paste("--program", programm, sep=" ")
  bootstraps <- paste("--bootstraps", bootstrap, sep=" ")
  outdir <- paste(tempdir(), "guidance", sep="/")
  outDir <- paste("--outDir", outdir, sep=" ")


  if(!missing(moreArgs))
    moreArgs <- paste("--MSA_Param", moreArgs, sep =" ")
  if(!missing(gencode))
    gencode <- paste("--genCode", gencode, sep=" ")
  if(!missing(outorder))
    outorder <- paste("--outOrder", outorder, sep=" ")
  if(!missing(msafile))
    msafile <- paste("--msaFile", msafile, sep=" ")
  if(!missing(cutoff))
    cutoff <- paste("--seqCutoff", cutoff, sep=" ")
  if(!missing(proc_num))
    proc_num <- paste("--proc_num", proc_num, sep=" ")


  moreArgs <- foreach(i = c("moreArgs", "gencode","outorder",
    "msafile", "cutoff", "proc_num"), .combine = 'c') %do% {
      if(is.object(i)){
        i
      }else{""}
    }
  moreArgs <- moreArgs[-grep("", moreArgs)]

  guidance.call <- paste(seqFile, msa.program, seqType,
    outDir, programm, bootstraps, moreArgs)
  perl.call <- paste("perl", paste(guidance_dir, "www/Guidance/guidance.pl" , sep=""))
  guidance.call <- paste(perl.call, guidance.call)
  ## CALL
  if(quiet){
    system(guidance.call, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }else{
    system(guidance.call)
  }

  files <- list.files(paste(tempdir(), "guidance", sep="/"), full.names = TRUE)
  # files <- list.files("../../../PhD/proj/high_priority/color_new_alignment/rpb1.fguidance/",full.names = TRUE)
  read <- files[grep("\\.scr\\b", files)]
  read <- read[-grep("csv", read)]

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

  # base.msa
  base <- files[grep("MSA.MAFFT.aln.With_Names", files)]
  base.msa <- read.fas(base, type =type)

  # guidance.msa
  guidance.msa <- files[grep("Without_low_SP_Col.With_Names", files)]
  guidance.msa <- read.fas(guidance.msa, type =type)

  ## delete temp files
  unlink(fns[file.exists(fns)], recursive = TRUE)
  unlink(files, force = TRUE)
  files <- list.files(tempdir(), full.names = TRUE)
  files <- files[-grep("rs-graphics", files)]
  unlink(files, force = TRUE, recursive = TRUE)


  # output
  res <- list(mean_score = msc,
    column_score = CS,
    residue_pair_column_scoore = g.sc,
    residue_pair_residue_score = rpr.sc,
    residual_pair_sequence_pair_score  = rpsp.sc,
    residual_pair_sequence_score = rps.sc,
    residue_pair_score = rp.sc,
    reduced_msa = guidance.msa,
    base_msa = base.msa)
  return(res)
}
