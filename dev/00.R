library("ips")
library("foreach")
library("parallel")
library("doSNOW")
library(pbmcapply)
# 1.	Read sequences
ms <- read.fas("dev/data/cortinarius_28s_ms.fas")
# ms <- read.FASTA("data/_28s.fas")
set.seed(100)
ms_test <- sample(ms, 10)


dna = ms_test
ncore = 2
bootstrap = 100
msa.program = "mafft"
method = "auto"
parallel  = TRUE
cutoff = 0.93

system.time(test <- guidance(dna = ms_test,
                                ncore = 2,
                                bootstrap = 100,
                                msa.program = "mafft",
                                method = "auto",
                                parallel  = TRUE,
                                cutoff = 0.93))
# 15 min with 2 cores
