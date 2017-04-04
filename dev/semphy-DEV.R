
# for DNA
msa <- base.msa
semphy <- function(msa, ML, dist.table, model, arsv, categories, program.dir)

  fns <- vector(length = 4)
  for (i in seq_along(fns))
    fns[i] <- tempfile(pattern = "semphy", tmpdir = tempdir(), fileext = ".fas")
  unlink(fns[file.exists(fns)])

  write.fas(msa, fns[1])

  type <- class(sequences)
  type <- gsub("bin", "", type)
  if(type == "DNA")
    a <- 4
  if(type == "AA")
    a <- 20

  ML <- ifelse(!missing(ML), "-S", "")

  if(missing(dist.table)) stop("dist.table must be specified, one of:\n posteriorDTME (iterative NJ) \n homogeneousRatesDTME (homogeneous rates model) \n rate4siteDTME (ML rate for each site)")
  if(dist.table == "auto")
    dist.table <- "--posteriorDTME"


  if(missing(model)) stop("model must be specified, one of:\n JC for nucleotides (--nucjc) or amino acids (--aaJC) \n K2P (--k2p) \n HKY (--hky) \n Dayhoff (--day) \n JTT (--jtt) \n REV (--rev) \n WAG (--wag) \n cpREV (--cprev)")
  if(model == "auto")
    if(type == 4){
      model <- "--hky"
    }else{ model <- "--jtt" }

  if(missing(arsv)) stop("Among rate site variation must be specified, one of: -H (homogeneous rates), -A (Set the initial alpha parameter), or -O (Optimize the alpha parameter)")
  if(arsv == "auto")
    arsv <- "-O"

  program.dir <- "/Applications/semphy-2.0b2-bin-mac-intel"

  if(missing(MoreArgs))
    MoreArgs <- ""

  semphy.call <- paste(paste(program.dir, "semphy", sep="/"),
    "-s", fns[1], "-o", fns[2], "-T", fns[3],
    "-l", fns[4],"-a", a, model, ML, dist.table, arsv, MoreArgs, sep=" ")

  system(semphy.call)
