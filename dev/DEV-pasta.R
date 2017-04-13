library("ips")
library("foreach")
library("parallel")
library("doSNOW")
library("pbmcapply")

library(rpg)

# 1.	Read sequences
seqs <- read.fas("dev/data/cortinarius_28s_ms.fas")
# set.seed(100)
# seq <- sample(ms, 10) ## choose 10 sequences

k = 50
ncore = 24
bootstrap = 100
msa.program = "mafft"
method = "auto"
parallel  = TRUE
exec = "/usr/local/bin/mafft"

system.time(test <- pasta(seqs = ms,
                          k = 50,
                          ncore = 24,
                          bootstrap = 100,
                          msa.program = "mafft",
                          method = "auto",
                          parallel = TRUE,
                          cutoff = 0.93))
# 15 min with 2 cores


seqs <- read.fas(x = "dev/data/AATF.fas")
