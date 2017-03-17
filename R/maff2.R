# Hi Christoph!
# wir müssten ips::mafft so anpassen das es auch mit AA sequenzen klar kommt
# ich habe es unten mal so angepasst das es läuft und die Stellen
# markiert.
# schicke auch einen beispieldatensatz mit.

## Gleiches gilt dann natürlich auch für PRANK


# suppressWarnings(
#   seq <- seqinr::read.fasta("dev/data/ref1_bb11005.fas", seqtype = "AA")
# )
#
# seq2 <- lapply(seq, as.character)
# seq2 <- rbind.fill(lapply(seq2,function(y) { as.data.frame(t(y), stringsAsFactors=FALSE) }))
# aa <- as.AAbin(as.matrix(seq2))
#
# base.msa <- mafft2(x = aa)

mafft2 <- function(x, y, add, method = "auto", maxiterate = 0, op = 1.53,
  ep = 0, gt, options, thread = -1, path, quiet)
{
  ###############################################
  if (!inherits(x, c("DNAbin", "AAbin")))
    stop("'x' is not of class 'DNAbin' or 'AAbin'")
  ###############################################
  os <- .Platform$OS
  if (missing(quiet))
    quiet <- TRUE
  qut <- ifelse(quiet, " --quiet ", " ")
  if (missing(path))
    path <- "/usr/local/bin/mafft"
  maxiterate <- match.arg(as.character(maxiterate), c("0",
    "2", "1000"))
  fns <- vector(length = 3)
  for (i in seq_along(fns)) fns[i] <- tempfile(pattern = "mafft",
    tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])
  method <- match.arg(method, c("auto", "localpair", "globalpair",
    "genafpair", "parttree", "retree 1", "retree 2"))
  if (missing(gt)) {
    gt <- ""
  }
  else {
    if (!inherits(gt, "phylo"))
      stop("object \"gt\" is not of class \"phylo\"")
    if (!all(names(x) %in% gt$tip.label))
      stop("guide tree does not match sequence names")
    gt$tip.label <- match(names(x), gt$tip.label)
    if (!is.binary.tree(gt))
      gt <- multi2di(gt)
    if (is.null(gt$edge.length))
      gt$edge.length <- rep(1, nrow(gt$edge))
    phylo2mafft(gt)
    gt <- " --treein tree.mafft "
  }
  thread <- paste("--thread", thread)
  thread <- paste(rep(" ", 2), collapse = thread)
  if (missing(options)) {
    options <- " "
  }
  else {
    options <- match.arg(options, c("--adjustdirection",
      "--adjustdirectionaccurately"))
    options <- paste(options, collapse = " ")
    options <- paste(rep(" ", 2), collapse = options)
  }
  if (missing(y)) {

    ##############################################
    if (inherits(x, "DNAbin")){ write.fas(x, fns[1])  }
    if (inherits(x, "AAbin")) {     write.phyDat(x, file = fns[1], format="fasta") }
    ###############################################

    call.mafft <- paste(path, " --", method, " --", "maxiterate ",
      maxiterate, qut, "--op ", op, " --ep ", ep, gt, options,
      thread, fns[1], " > ", fns[3], sep = "")
  }
  else {
    if (!inherits(y, "DNAbin"))
      stop("'y' is not of class 'DNAbin'")
    if (missing(add))
      add <- "addprofile"
    add <- match.arg(add, c("add", "addprofile"))
    add <- paste("--", add, sep = "")
    write.fas(x, fns[1])
    write.fas(y, fns[2])
    call.mafft <- paste(path, qut, add, fns[2], fns[1], ">",
      fns[3])
  }
  if (!quiet)
    message(call.mafft)
  if (os == "unix") {
    system(call.mafft, intern = FALSE, ignore.stdout = FALSE)
    res <- length(scan(fns[3], what = "c", quiet = TRUE))
    if (res != 0)
      ###############################################
    if (inherits(x, "DNAbin")) { res <- read.fas(fns[3]) }
    if (inherits(x, "AAbin" )) { res <- seqinr::read.fasta(fns[3], seqtype  ="AA")       ###############################################
    }
  }
  else {
    res <- system(call.mafft, intern = TRUE, ignore.stderr = FALSE)
    if (length(grep("error|ERROR", res)) > 0) {
      res <- 0
    }
    else {
      ###############################################
      if (inherits(x, "DNAbin")) { res <- read.fas(fns[3]) }
      if (inherits(x, "AAbin" )) { res <- seqinr::read.fasta(fns[3], seqtype  ="AA")
      ###############################################
    }
  }
  unlink(fns[file.exists(fns)])
  return(res)
}
