library("rpg")
library("ips")
library("parallel")
library("foreach")
library("phangorn")
library("doSNOW")
library("adephylo")

# 1.	Read example sequences
## DNA
ms <- read.fas("dev/data/cortinarius_28s_ms.fas", type ="DNA")
set.seed(100)
seq_dna <- sample(ms, 10)
## Amino Acids
seq_aa <- read.fas("dev/data/AATF.fas", type ="AA")


sequences = seq_aa
parallel = TRUE
ncore = 4
bootstrap = 10
col.cutoff = "auto"
seq.cutoff = "auto"
mask.cutoff = "auto"
msa.program <- "clustalo"
exec <- "/Applications/clustalo"
msa.program <- "mafft"
exec <- "/usr/local/bin/mafft"
msa.program <- "clustalw2"
exec <- "/Applications/clustalw2"
msa.program <- "muscle"
exec <- "/Applications/muscle"

system.time(g_msa <- guidance(sequences = seq_dna,
  msa.program = "muscle",
  exec = "/Applications/muscle",
  bootstrap = 10,
  col.cutoff = "auto",
  seq.cutoff = "auto",
  mask.cutoff = "auto",
  parallel = TRUE, ncore = "auto",
  method = "auto"))

heatmap.msa(g_msa)


##
sequences <- seq_dna
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
n.coopt = "auto"

system.time(hot_msa <- HoT_dev(sequences = seq_dna, # MUSCLE stops after Scores
  msa.program = "muscle",
  exec = "/Applications/muscle",
  n.coopt = "auto",
  col.cutoff = "auto",
  seq.cutoff = "auto",
  mask.cutoff = "auto",
  parallel = TRUE, ncore = 4,
  method = "auto",
  plot_guide = TRUE,
  alt.msas.file))

system.time(hot_msa <- HoT_dev(sequences = seq_dna,
  ncore = 4,
  msa.program = "prank",
  parallel  = TRUE))

heatmap.msa(obj = hot_msa)






##
sequences = seq_dna
msa.program = "muscle"
# exec <- "/usr/local/bin/mafft"
exec <- "/Applications/muscle"
bootstrap = 10
n.part="auto"
col.cutoff = "auto"
seq.cutoff = "auto"
mask.cutoff = "auto"
parallel = TRUE
ncore =4
method = "auto"
# alt.msas.file
n.coopt = "auto"

g2 <- guidance2(sequences = seq_dna,
  msa.program = "clustalw2",
  exec <- "/Applications/clustalw2",
  bootstrap = 10,
  n.part="auto",
  col.cutoff = "auto",
  seq.cutoff = "auto",
  mask.cutoff = "auto",
  parallel = TRUE, ncore ="auto",
  method = "auto",
  n.coopt = "auto",
  alt.msas.file = paste(getwd(), "R", "test", sep="/"))

heatmap.msa(obj = g2)
