# [FK 2017-03-27]
# Code based on 'muscle' from ape
# NEW: alignment of PROFILE1 and PROFILE2
# NEW: guide-tree alignment
#' @export

muscle2 <- function (x, y,gt, exec = "muscle", MoreArgs = "",
  quiet = TRUE, original.ordering=TRUE, file)
{
  if (missing(x)) {
    out <- system(exec)
    if (out == 127)
      stop(.errorAlignment(exec, "MUSCLE"))
    return(invisible(NULL))
  }

  type <- class(x)
  type <- gsub("bin", "", type)

  # Produce TEMP files
  fns <- vector(length = 4)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "muscle", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])
  # gtt <- tempfile(pattern = "guidetree", tmpdir = tempdir(), fileext = ".nwk")

  # Write input sequences X to file
  x <- as.list(x)
  labels.bak <- names(x)
  names(x) <- paste0("Id", 1:length(x))
  write.fas(x, fns[1])

  # Run muscle for X
  if(missing(y)){
    opts <- paste("-in", fns[1], "-out", fns[3])
    if (quiet)
      opts <- paste(opts, "-quiet")
    opts <- paste(opts, MoreArgs)

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
      gt <- paste("-usetree_nowarn ", fns[4], sep="")
      opts <- paste(opts, gt)
    }

    out <- system(paste(exec, opts), ignore.stdout = TRUE,  ignore.stderr = TRUE)
    if (out == 127)
      stop(.errorAlignment(exec, "MUSCLE"))
    res <- read.fas(fns[3], type =type)
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
    opts <- paste("-profile", "-in1", fns[1],"-in2", fns[2], "-out", fns[3])
    if (quiet)
      opts <- paste(opts, "-quiet")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts),ignore.stdout = TRUE,  ignore.stderr = TRUE)
    if (out == 127)
      stop(.errorAlignment(exec, "MUSCLE"))
    res <- read.fas(fns[3], type =type)
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
