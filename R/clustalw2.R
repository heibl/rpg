# [FK 2017-03-27]
# Code based on 'clustal'  from ape
# NEW: alignment of PROFILE1 and PROFILE2
# NEW: guide-tree alignment
#' @export

clustalw2 <- function (x, y, gt, pw.gapopen = 10, pw.gapext = 0.1, gapopen = 10,
  gapext = 0.2, exec = NULL, MoreArgs = "", quiet = TRUE, original.ordering = TRUE, file)
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

  type <- class(x)
  type <- gsub("bin", "", type)

  fns <- vector(length = 4)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "clustal", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])
  # gtt <-  tempfile(pattern = "guidetree", tmpdir = tempdir(), fileext = ".nwk")

  x <- as.list(x)
  labels.bak <- names(x)
  names(x) <- paste0("Id", 1:length(x))
  write.fas(x, fns[1])


  if(missing(y)){
    prefix <- c("-INFILE", "-PWGAPOPEN", "-PWGAPEXT", "-GAPOPEN","-GAPEXT", "-OUTFILE","-OUTPUT")
    suffix <- c(fns[1], pw.gapopen, pw.gapext, gapopen, gapext, fns[3], "FASTA")
    opts <- paste(prefix, suffix, sep = "=", collapse = " ")

    # add input guide tree
    if (!missing(gt)) {
      if (!inherits(gt, "phylo"))
        stop("object \"gt\" is not of class \"phylo\"")
      if (!all(labels.bak %in% gt$tip.label))
        stop("guide tree does not match sequence names")
      gt$tip.label <- paste0("Id", 1:length(x))
      if (!is.binary.tree(gt))
        gt <- multi2di(gt)
      if (is.null(gt$edge.length))
        gt$edge.length <- rep(1, nrow(gt$edge))
      write.tree(gt, fns[4])
      gt <- paste("-USETREE=", fns[4], sep="")
      opts <- paste(opts, gt)
    }

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
  unlink(fns[file.exists(fns)])
  # unlink(gtt[file.exists(gtt)])
  if(!missing(file)){
    write.fas(res, file)
  }else{
  return(res)
  }
}
