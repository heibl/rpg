library("ips")
library("foreach")
library("parallel")
library("doSNOW")
library("pbmcapply")
library("phangorn")
library("plyr")
library("ggplot2")

# 1.	Read example sequences
## DNA
ms <- read.fas("dev/data/cortinarius_28s_ms.fas")
set.seed(100)
seq <- sample(ms, 200)

## Amino Acids
suppressWarnings(
  seq <- seqinr::read.fasta("dev/data/AATF.fas", seqtype = "AA")
)


seq = seq
ncore = 2
bootstrap = 100
msa.program = "mafft"
method = "auto"
parallel  = TRUE
cutoff = 0.93
mask = FALSE

# system.time(test <- guidance(seq = seq,
#                                 ncore = 2,
#                                 bootstrap = 100,
#                                 msa.program = "mafft",
#                                 method = "auto",
#                                 parallel  = TRUE,
#                                 cutoff = 0.93,
#                                 mask = FALSE))

# guidance_heatmap(test, file ="dev/data/test_msa_plot.pdf")
# guidance_heatmap(test)
