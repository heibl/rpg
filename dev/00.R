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
library("rpg")

# 1.	Read example sequences
## DNA
ms <- read.fas("dev/data/cortinarius_28s_ms.fas", type ="DNA")
set.seed(100)
seq_dna <- sample(ms, 10)
## Amino Acids
seq_aa <- read.fas("dev/data/AATF.fas", type ="AA")


# seq = seq_dna
sequences = seq_aa
ncore = 4
bootstrap = 2
msa.program = "mafft"
method = "auto"
parallel  = TRUE
cutoff = 0.93
mask = FALSE
plot_guide = TRUE
n.part = "auto"
exec <- "/usr/local/bin/mafft"
exec <- "/Applications/clustalo"
exec <- "/Applications/clustalw2"
exec <- "/Applications/muscle"

seq_dna <- read.fas("/Users/krah/Dropbox/color_gradient/data/new_phylo_march17/selected_msa/rpb1.fas", type ="DNA")
system.time(g_msa <- guidance(seq = seq_dna,
                                bootstrap = 100,
                                msa.program = "mafft",
                                exec = "/usr/local/bin/mafft",
                                method = "auto",
                                parallel  = TRUE,
                                ncore = 4,
                                cutoff = 0.93,
                                mask = FALSE))

# heatmap.msa(g_msa, file ="dev/data/test_msa_plot.pdf")
# system("open dev/data/test_msa_plot.pdf")
heatmap.msa(g_msa)

msa.program = "mafft"
exec <- "/usr/local/bin/mafft"
n.part = "auto"
col.cutoff = "auto"
seq.cutoff = "auto"
mask.cutoff = "auto"
parallel = FALSE
ncore = 4
method = "auto"
plot_guide = TRUE

system.time(hot_msa <- HoT_dev(sequences = seq_aa, # MUSCLE stops after Scores
  msa.program = "muscle",
  exec = "/Applications/muscle",
  n.part = "auto",
  col.cutoff = "auto",
  seq.cutoff = "auto",
  mask.cutoff = "auto",
  parallel = TRUE, ncore = 4,
  method = "auto",
  plot_guide = TRUE))

system.time(hot_msa <- HoT_dev(sequences = seq_dna,
  ncore = 4,
  msa.program = "prank",
  parallel  = TRUE))

heatmap.msa(obj = hot_msa)
