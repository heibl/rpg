library(ips); library(foreach);
library(parallel)
library(doSNOW)
library(pbmcapply)

# 1.	Read sequences
ms <- read.fas("dev/data/cortinarius_28s_ms.fas")
# ms <- read.FASTA("data/_28s.fas")
set.seed(100)
ms_test <- sample(ms, 10)
# seqinr::write.fasta(as.character(ms_test), names = names(ms_test),
#   file = "dev/out/ms_test.fas")

# ms_test <- read.FASTA("out/msa_test.fas")


dna = ms_test
ncore = 24
bootstrap = 100
msa.programm = "mafft"
method = "auto"
parallel  = TRUE
G.CS.th = 0.93

# system.time(test <- guidance(dna = ms_test,
#                              ncore = 12,
#                              bootstrap = 100,
#                              msa.programm = "mafft",
#                              parallel  = TRUE,
#                              G.CS.th = 0.93))
#
