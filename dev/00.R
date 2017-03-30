library("rpg")

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

system.time(g_msa <- guidance(sequences = seq_aa,
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
