library("rpg")
library("ips")
library("parallel")
library("foreach")
library("phangorn")
library("doSNOW")
library("adephylo")
library("useful")
library("stringr")
library("scales")
library("ggplot2")
library("zoo")
library(cowplot)
# 1.	Read example sequences
## DNA
seq_dna <- read.fas("dev/data/cortinarius_28s_ms.fas", type ="DNA")
set.seed(100)
seq_dna <- sample(seq_dna,10)
## Amino Acids
seq_aa <- read.fas("dev/data/AATF.fas", type ="AA")


# msa.program <- "mafft"
# exec <- "/usr/local/bin/mafft"
# msa.program <- "clustalo"
# exec <- "/Applications/clustalo"
# msa.program <- "clustalw2"
# exec <- "/Applications/clustalw2"
# msa.program <- "muscle"
# exec <- "/Applications/muscle"

system.time(g_r <- guidance(sequences = seq_aa,
  msa.program = "mafft",
  # exec = exec,
  bootstrap = 100,
  col.cutoff = "auto",
  seq.cutoff = "auto",
  mask.cutoff = "auto",
  parallel = TRUE, ncore = "auto",
  method = "retree 1"))

system.time(
  g_sa <- guidanceSA(sequences = seq_dna,
  msa.program = "mafft",
  programm = "guidance",
  bootstrap = 100,
  proc_num = 4)
)
g_r.h <- confidence.heatmap(g_r, title = "GUIDANCE R", legend = FALSE,
  guidance_score = FALSE)
g_sa.h <- confidence.heatmap(g_sa, title = "GUIDANCE SA", legend = FALSE,
  guidance_score = FALSE)


sequences = seq_aa
msa.program = "mafft"
parallel = FALSE
ncore = 4
method = "retree 1"
plot_guide = TRUE

system.time(hot_r <- HoT(sequences = seq_aa,
  msa.program = "mafft",
  parallel = FALSE, ncore = 4,
  method = "retree 1",
  plot_guide = TRUE))

system.time(
  hot_sa <- guidanceSA(sequences = seq_aa,
    msa.program = "mafft",
    programm = "hot",
    bootstrap = 100,
    proc_num = 4,
    quiet = FALSE)
)
hot_r.h <- confidence.heatmap(hot_msa, title = "HoT SA", legend = FALSE,
  guidance_score = FALSE)
hot_sa.h <- confidence.heatmap(hot_sa, title = "HoT SA", legend = FALSE,
  guidance_score = FALSE)
plot_grid(hot_r.h, hot_sa.h, nrow = 2, ncol = 1)







system.time(
g2_r <- guidance2(sequences = seq_aa,
  msa.program = "mafft",
  # exec <- "/Applications/clustalw2",
  bootstrap = 100,
  parallel = TRUE, ncore ="auto",
  method = "auto",
  n.coopt = "auto")
)
heatmap.msa(obj = g2, file =paste(getwd(), "test.pdf", sep="/"))


system.time(g2_sa <- guidanceSA(sequences = seq_aa,
  msa.program = "mafft",
  outdir = "test.guidance",
  programm = "guidance2",
  bootstrap = 100,
  proc_num = 4))
