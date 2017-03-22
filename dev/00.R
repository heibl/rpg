library("ips")
library("foreach")
library("parallel")
library("doSNOW")
library("pbmcapply")
library("phangorn")
library("plyr")
library("ggplot2")
library("scales")
library("adephylo")

# 1.	Read example sequences
## DNA
ms <- read.fas("dev/data/cortinarius_28s_ms.fas", type ="DNA")
set.seed(100)
seq_dna <- sample(ms, 10)
## Amino Acids
seq_aa <- read.fas("dev/data/AATF.fas", type ="AA")


# seq = seq_dna
seq = seq_aa
ncore = 2
bootstrap = 100
msa.program = "mafft"
method = "auto"
parallel  = TRUE
cutoff = 0.93
mask = FALSE
mafft_exec <- "/usr/local/bin/mafft"
plot_guide = TRUE

# system.time(g_msa <- guidance(seq = seq,
#                                 ncore = 2,
#                                 bootstrap = 100,
#                                 msa.program = "mafft",
#                                 method = "auto",
#                                 parallel  = TRUE,
#                                 cutoff = 0.93,
#                                 mask = FALSE))

# guidance_heatmap(g_msa, file ="dev/data/test_msa_plot.pdf")
# system("open dev/data/test_msa_plot.pdf")
# guidance_heatmap(guidance_obj = g_msa)

library("rpg")
system.time(hot_msa <- HoT(seq = seq,
  ncore = 2,
  msa.program = "mafft",
  method = "auto",
  parallel  = FALSE,
  cutoff = 0.93,
  mask = FALSE,
  plot_guide = TRUE))
