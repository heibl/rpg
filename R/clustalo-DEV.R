# [FK 2017-03-27]
# Code based on 'clustalomega'  from ape
# NEW: alignment of PROFILE1 and PROFILE2
#' @export

clustalo <- function (x, y, exec = NULL, type,MoreArgs = "", quiet = TRUE)
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

  fns <- vector(length = 4)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "clustalo", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])


  x <- as.list(x)
  labels.bak <- names(x)
  names(x) <- paste0("Id", 1:length(x))
  write.fas(x, fns[1])

  if (!quiet)
    opts <- paste(opts, "-v")

  if(missing(y)){
    opts <- paste("-i", fns[1], "-o", fns[2], "--force")
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127)
      stop(.errorAlignment(exec, "Clustal-Omega"))
    res <- read.fas(fns[2], type =type)
    rownames(res) <- labels.bak
    res
  }


  if(!missing(y)){
    y <- as.list(y)
    labels.baky <- names(y)
    write.fas(y, fns[3])
    names(y) <- paste0("Id", (length(x)+1):(length(x)+length(y)))
    if(length(y)==1){
      opts <- paste("-i", fns[1],"--profile1",
        fns[3], "-o", fns[4], "--force")
    }else{
      opts <- paste("--profile1", fns[1],"--profile2",
        fns[3], "-o", fns[4], "--force")
    }
    opts <- paste(opts, MoreArgs)
    out <- system(paste(exec, opts), ignore.stdout = quiet)
    if (out == 127)
      stop(.errorAlignment(exec, "Clustal-Omega"))
    res <- read.fas(fns[4], type = type)
    rownames(res) <- c(labels.bak, labels.baky)
    res
  }


  unlink(fns[file.exists(fns)])
  return(res)
}
