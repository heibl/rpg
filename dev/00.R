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
# 1.	Read example sequences
## DNA
seq_dna <- read.fas("dev/data/cortinarius_28s_ms.fas", type ="DNA")
set.seed(100)
seq_dna <- sample(seq_dna,10)
## Amino Acids
seq_aa <- read.fas("dev/data/AATF.fas", type ="AA")

file <- system.file("extdata", "BB50009.fasta", package = "rpg")
aa_seq<- read.fas(file, type ="AA")

g_msa <- guidance(sequences = aa_seq,
  msa.program = "mafft",
  exec = "/usr/local/bin/mafft",
  bootstrap = 100,
  parallel = FALSE,
  method = "retree 1")
confidence.heatmap(g_msa, title = "GUIDANCE BB30015 benchmark", legend = TRUE,
  guidance_score = TRUE)

sequences = seq_dna
parallel = TRUE
ncore = 4
bootstrap = 10
col.cutoff = "auto"
seq.cutoff = "auto"
mask.cutoff = "auto"
method <- "retree 1"
msa.program <- "mafft"
exec <- "/usr/local/bin/mafft"

msa.program <- "clustalo"
exec <- "/Applications/clustalo"
msa.program <- "clustalw2"
exec <- "/Applications/clustalw2"
msa.program <- "muscle"
exec <- "/Applications/muscle"

system.time(g_msa <- guidance(sequences = seq_aa,
  msa.program = "mafft",
  # exec = exec,
  bootstrap = 100,
  col.cutoff = "auto",
  seq.cutoff = "auto",
  mask.cutoff = "auto",
  parallel = FALSE, ncore = "auto",
  method = "retree 1"))

confidence.heatmap(g_msa, title = "GUIDANCE R", legend = TRUE,
  guidance_score = TRUE)

system.time(
guidanceSA(sequences = seq_dna,
  msa.program = "mafft",
  programm = "guidance",
  bootstrap = 100,
  proc_num = 4)
)

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

system.time(hot_msa <- HoT_dev(sequences = seq_aa, # MUSCLE stops after Scores
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

system.time(
  hot_sa <- guidanceSA(sequences = seq_dna,
    msa.program = "mafft",
    programm = "hot",
    bootstrap = 100,
    proc_num = 4)
)

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

system.time(
g2 <- guidance2(sequences = seq_aa,
  msa.program = "mafft",
  # exec <- "/Applications/clustalw2",
  bootstrap = 10,
  n.part="auto",
  col.cutoff = "auto",
  seq.cutoff = "auto",
  mask.cutoff = "auto",
  parallel = TRUE, ncore ="auto",
  method = "auto",
  n.coopt = "auto")
)
heatmap.msa(obj = g2, file =paste(getwd(), "test.pdf", sep="/"))




##

# seq_dna <- read.fas("Dropbox/cortinarius_28s_ms.fas", type ="DNA")
# # set.seed(100)
# # seq_dna <- sample(sea_dna, 10)
# sequences = seq_dna
# msa.program = "mafft"
# outdir = "test.guidance"
# programm = "GUIDANCE"
# bootstrap = 100
# proc_num = 24
# gencode
# outorder
# msafile

# system.time(gSA <- guidanceSA(sequences = seq_dna,
#   msa.program = "mafft",
#   outdir = "test.guidance",
#   programm = "GUIDANCE",
#   bootstrap = 100,
#   proc_num = 12))
