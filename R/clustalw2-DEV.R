# [FK 2017-03-27]
# Code based on 'clustal'  from ape
# NEW: alignment of PROFILE1 and PROFILE2
#' @export

clustalw2 <- function (x, y, type, pw.gapopen = 10, pw.gapext = 0.1, gapopen = 10,
  gapext = 0.2, exec = NULL, MoreArgs = "", quiet = TRUE, original.ordering = TRUE)
{
  os <- Sys.info()[1]
  if (is.null(exec)) {
    exec <- switch(os, Linux = "clustalw", Darwin = "clustalw2",
      Windows = "clustalw2.exe")
  }
  if (missing(x)) {
    out <- system(paste(exec, "-help"))
    if (out == 127)
      stop(.errorAlignment(exec, "Clustal"))
    return(invisible(NULL))
  }
  fns <- vector(length = 4)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "clustal", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  x <- as.list(x)
  labels.bak <- names(x)
  names(x) <- paste0("Id", 1:length(x))
  write.fas(x, fns[1])


  if(missing(y)){
    prefix <- c("-INFILE", "-PWGAPOPEN", "-PWGAPEXT", "-GAPOPEN","-GAPEXT", "-OUTFILE","-OUTPUT")
    suffix <- c(fns[1], pw.gapopen, pw.gapext, gapopen, gapext, fns[3], "FASTA")
    opts <- paste(prefix, suffix, sep = "=", collapse = " ")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127)
      stop(.errorAlignment(exec, "Clustal"))
    res <- read.fas(fns[3], type = type)
    if (original.ordering)
      res <- res[labels(x), ]
    rownames(res) <- labels.bak
    res
  }

  if(!missing(y)){
    y <- as.list(y)
    labels.baky <- names(y)
    names(y) <- paste0("Id", (length(x)+1):(length(x)+length(y)))
    write.fas(y, fns[2])
    # if(length(y)==1){
    #   prefix <- c("-INFILE", "-PROFILE1", "-PWGAPOPEN", "-PWGAPEXT",
    #     "-GAPOPEN","-GAPEXT", "-OUTFILE","-OUTPUT")
    #   suffix <- c(fns[1], fns[2], pw.gapopen, pw.gapext,
    #     gapopen, gapext, pns[1], "PHYLIP")
    # }else{
    prefix <- c("-PROFILE1", "-PROFILE2", "-PWGAPOPEN", "-PWGAPEXT",
      "-GAPOPEN","-GAPEXT", "-OUTFILE","-OUTPUT")
    suffix <- c(fns[1], fns[2], pw.gapopen, pw.gapext,
      gapopen, gapext, fns[3], "FASTA")
    # }
    opts <- paste(prefix, suffix, sep = "=", collapse = " ")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
  if (out == 127)
    stop(.errorAlignment(exec, "Clustal"))
  res <- read.fas(fns[3], type = type)
  if (original.ordering)
    res <- res[c(labels(x), labels(y)), ]
  rownames(res) <- c(labels.bak, labels.baky)
  res
  }
  return(res)
}
