library("rpg")
# seq_aa <- read.fas("dev/data/AATF.fas", type ="AA")
seq_dna <- read.fas("dev/data/cortinarius_28s_ms.fas", type ="DNA")
set.seed(100)
seq_dna <- sample(seq_dna, 200)

seq_aa <- read.fas("dev/data/bb3_release/RV40/BB40002.tfa", type ="aa")


bootstrap <- 100
msa.program <- "mafft"
exec <- "/usr/local/bin/mafft"
# exec <- "/Applications/clustalo"
# exec <- "/Applications/clustalw2"
# exec <- "/Applications/muscle"


## GUIDANCE
g_time <- system.time(g <- guidance(sequences = seq_dna,
  msa.program = msa.program,
  exec = exec,
  bootstrap = bootstrap,
  parallel = TRUE, ncore = "auto",
  method = "retree 1"))
g_h <- heatmap.msa(g)

gSA_time <- system.time(
  gSA <- guidanceSA(sequences = seq_dna,
    msa.program = "mafft",
    programm = "guidance",
    bootstrap = bootstrap,
    proc_num = 4)
)
gSA_h <- heatmap.msa(gSA)

## HoT
h_time <- system.time(h <- HoT_dev(sequences = seq_dna,
  msa.program = msa.program,
  exec = exec,
  parallel = TRUE,
  ncore = "auto",
  plot_guide = FALSE))
h_h <- heatmap.msa(h)

hSA_time <- system.time(
  hSA <- guidanceSA(sequences = seq_dna,
    msa.program = "mafft",
    programm = "hot",
    bootstrap = bootstrap,
    proc_num = 4)
)
hSA_h <- heatmap.msa(hSA)

## GUIDANCE2
g2_time <- system.time(
g2 <- guidance2(sequences = seq_dna,
  msa.program = msa.program,
  exec = exec,
  bootstrap = bootstrap,
  parallel = TRUE,
  method = "retree 1",
  ncore ="auto"))
g2_h <- heatmap.msa(obj = g2)


g2SA_time <- system.time(
  g2SA <- guidanceSA(sequences = seq_dna,
    msa.program = "mafft",
    programm = "guidance2",
    bootstrap = bootstrap,
    proc_num = 4)
)
g2SA_h <- heatmap.msa(g2SA)



library("cowplot")
plot_grid(g_h, gSA_h, nrow = 2, ncol = 1)
plot_grid(g2_h, g2SA_h, nrow = 1, ncol = 2)

time_comp <- c(gtime[3], htime[3], g2time[3])
barplot(time_comp, names.arg = c("GUIDANCE", "HoT", "GUIDANCE2"))

