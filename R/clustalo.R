# [FK 2017-03-30]
# Code based on 'clustalomega'  from ape
# NEW: alignment of PROFILE1 and PROFILE2
# NEW: guide-tree alignment
#' @export

clustalo <- function (x, y, gt, exec = NULL,MoreArgs = "",
  quiet = TRUE, original.ordering = TRUE, file)
{
  os <- Sys.info()[1]
  if (is.null(exec)) {
    exec <- switch(os, Linux = "clustalo", Darwin = "clustalo",
      Windows = "clustalo.exe")
  }
  if (missing(x)) {
    out <- system(paste(exec, "-h"))
    if (out == 127)
      stop(.errorAlignment(exec, "Clustal-Omega"))
    return(invisible(NULL))
  }

  type <- class(x)
  type <- gsub("bin", "", type)

  fns <- vector(length = 4)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "clustal", tmpdir = tempdir(), fileext = ".fas")
  # gtt <-  tempfile(pattern = "guidetree", tmpdir = tempdir(), fileext = ".nwk")
  unlink(fns[file.exists(fns)])
  # unlink(gtt[file.exists(gtt)])

  x <- as.list(x)
  labels.bak <- names(x)
  names(x) <- paste0("Id", 1:length(x))
  write.fas(x, fns[1])

  if (!quiet)
    opts <- paste(opts, "-v")

  if (missing(y)){
    opts <- paste("-i", fns[1], "-o", fns[2], "--force")

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
      gt <- paste("--guidetree-in ", fns[4], sep="")
      opts <- paste(opts, gt)
    }
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127)
      stop(.errorAlignment(exec, "Clustal-Omega"))
    res <- read.fas(fns[2], type =type)
    if (original.ordering)
      res <- res[labels(x), ]
    rownames(res) <- labels.bak
    res
  }


  if(!missing(y)){
    y <- as.list(y)
    labels.baky <- names(y)
    write.fas(y, fns[2])
    names(y) <- paste0("Id", (length(x)+1):(length(x)+length(y)))
    if(length(y)==1){
      opts <- paste("-i", fns[1],"--profile1",
        fns[2], "-o", fns[3], "--force")
    }else{
      opts <- paste("--profile1", fns[1],"--profile2",
        fns[2], "-o", fns[3], "--force")
    }
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127)
      stop(.errorAlignment(exec, "Clustal-Omega"))
    res <- read.fas(fns[3], type = type)
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
