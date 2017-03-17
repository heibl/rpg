## Development script for function guidance
## Last update 2017-03-15

library("ips")
library("foreach")
library("parallel")
library("doSNOW")
library("pbmcapply")

library(rpg)

# 1.	Read sequences
ms <- read.fas("dev/data/cortinarius_28s_ms.fas")
set.seed(100)
seqs <- sample(ms, 10) ## choose 10 sequences

ncore = 24
bootstrap = 100
msa.program = "mafft"
method = "auto"
exec = "/usr/local/bin/mafft"
parallel  = TRUE
cutoff = 0.93

system.time(test <- guidance(seqs = seqs,
                             ncore = 24,
                             bootstrap = 100,
                             msa.program = "mafft",
                             method = "auto",
                             parallel = TRUE,
                             cutoff = 0.93))
# 15 min with 2 cores
