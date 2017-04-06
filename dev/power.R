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
library(bios2mds)
library(XML)



#### SOME CORE FUNCTIONS
## extract core blocks from XML files
read.coreBlocks <- function(file){
  doc <- xmlParse(file)
  doc <- xmlToList(doc)
  core_block <- doc$alignment$`column-score`$`colsco-data`
  core_block <- gsub("\n", "", core_block)
  core_block <- unlist(str_split(core_block, " "))
  core_block <- core_block[-which(nchar(core_block)==0)]
  return(core_block)
}


## Generate residue pairs residue scores from GUIDANCE outputs
# based on core blocks
rocr_mat <- function(true_msa, guidance_out, core_block){
  # REDUCTION based on GUIDANCE2 paper:
  # reduce dataset to contain only columns of
  # the MSA which have more than 2 residues
  # then select those column which are core blocks
  bmsa <- guidance_out$base_msa

  if(length(grep("scores", names(guidance_out)))>0)
    guidance_out <- guidance_out$scores

  true_msa <- as.character(true_msa)
    gaps <- apply(data.frame(true_msa), 2, grep, pattern = "-")
    true_msa <- true_msa[,!lapply(gaps, length)==nrow(true_msa)]
    true_msa <-  as.AAbin(true_msa)
  true_g <- compareMSAs(com = bmsa, ref = true_msa)

  # table with resides scores for true and test
  tbase <- as.vector(t(as.character(true_msa)))
  rbase <- as.vector(t(as.character(bmsa)))

  d_true <- cbind(true_g$residue_pair_residue_score, t.base = tbase)
  d_test <- cbind(guidance_out$residue_pair_residue_score, r.base = rbase)

  # remove gaps and add core blocks score
  d_true_red <- d_true[-grep("-", d_true$t.base),]
  d_true_red <- cbind(d_true_red, cb = core_block[match(d_true_red$col, core_block$col),]$cb)
  d_test_red <- d_test[-grep("-", d_test$r.base),]

  # combine residue score tables of true and test and core blocks score
  d <- data.frame(d_true_red, d_test_red)
  # remoove NAs
  d <- d[!is.na(d$score),]
  # keep only columns of the true alignment which have more than 2 residues
  keep <- table(d$col.1)
  keep <- as.numeric(names(keep[keep>1]))
  d <- d[d$col.1 %in% keep,]
  # and now keep only column of the test alignment which have more than 2 residues
  keep <- table(d$col)
  keep <- as.numeric(names(keep[keep>1]))
  d <- d[d$col %in% keep,]
  # reduce table to columns which are core blocks
  cols_for_analysis <- unique(d[d$cb==1,]$col)

  # generate residue pair score comparision once again
  # and reduce to cols selected as containing residues of core block column
  # of the true MSA
  true_g <- compareMSAs(ref = bmsa, com = true_msa)
  test_mat <- data.frame(true_g$residue_pair_score, guidance_out$residue_pair_score)
  test_mat <- test_mat[test_mat$col1 %in% cols_for_analysis,]
  test_mat <- test_mat[,c(1:4, 8)]
  names(test_mat)[4:5] <- c("true_score", "test_score")
  return(test_mat)
}



## PREPARE BENCHMARK DATASET
#------------------------------
files_bench <- list.files("dev/data/bb3_release/", full.names = TRUE, recursive = TRUE)
files_bench <- files_bench[-grep("README", files_bench)]

files_bench_f <- files_bench[grep("\\.tfa", files_bench)]
files_bench_xml <- files_bench[grep("\\.xml", files_bench)]
files_bench_msf <- files_bench[grep("\\.msf", files_bench)]

fas <- lapply(files_bench_f, read.fas, type = "AA")
g7 <- lapply(fas, length)>7

fas <- fas[g7]
files_bench_xml <- files_bench_xml[g7]
files_bench_msf <- files_bench_msf[g7]

msf <- lapply(files_bench_msf, import.msf)
for(z in 1:length(msf)){
  class(msf[[z]]) <- "list"
  msf[[z]] <- as.AAbin(do.call(rbind, msf[[z]]))
}
# cb <- lapply(files_bench_xml, read.coreBlocks)
# save(cb, file ="dev/res_benchmark/cb_>7.rda")
load(file ="dev/res_benchmark/cb_>7.rda")



#######################################################
#################   BENCHMARKING   ####################
#######################################################
#######################################################


#######################################################
######## COMPARISON BETWEEN GUIDANCE R vs. SA #########
#######################################################

## GUIDANCE BENCHMARKING START
#------------------------------------------------------
#######################################################

## GUIDANCE R
##################################
time <- list()
for(i in 1:length(fas)){
  g_time <- system.time(g <- guidance(sequences = fas[[i]],
    msa.program = "mafft",
    exec = "/usr/local/bin/mafft",
    bootstrap = 10,
    parallel = TRUE, ncore = "auto",
    method = "retree 1"))
  time[[i]] <- g_time

  coreblock <- data.frame(col = 1:length(cb[[i]]), cb = cb[[i]])
  g_corc <- rocr_mat(true_msa = msf[[i]], guidance_out = g, core_block = coreblock)
  if(dim(g_corc)[1]>0)
    write.csv(g_corc,
      file =paste("dev/res_balibase_1-5/GUIDANCE/gR_corc",
        str_extract(files_bench_f[g7], "BB[A-Z]*[0-9]*")[i],
        ".csv", sep=""))
}
time_g <- do.call(rbind, time)
save(time_g, file ="dev/res_balibase_1-5/GUIDANCE/time_g.rda")


## GUIDANCE SA
#######################

time <- list()
for(i in 1:length(fas)){
  gSA_time <- system.time(
    gSA <- guidanceSA(sequences = fas[[i]],
      msa.program = "mafft",
      programm = "guidance",
      bootstrap = 100,
      proc_num = 4)
  )
  time[[i]] <- gSA_time

  coreblock <- data.frame(col = 1:length(cb[[i]]), cb = cb[[i]])
  gSA_corc <- rocr_mat(true_msa = msf[[i]], guidance_out = gSA, core_block = coreblock)
  if(dim(gSA_corc)[1]>0)
    write.csv(gSA_corc,
      file =paste("dev/res_balibase_1-5/GUIDANCE/gSA_corc",
        str_extract(files_bench_f[g7], "BB[A-Z]*[0-9]*")[i],
        ".csv", sep=""))
}
time_gSA <- do.call(rbind, time)
save(time_gSA, file ="dev/res_balibase_1-5/GUIDANCE/time_gSA.rda")


####
library("ROCR")
g_rocr <- read.csv("dev/res_balibase_1-5/GUIDANCE/gR_balibase/gR_benchmark.csv")
gSA_rocr <- read.csv("dev/res_balibase_1-5/GUIDANCE/gSA_balibase/gSA_benchmark.csv")

pred_g <- prediction(as.numeric(as.character(g_rocr$test_score)),
  as.numeric(as.character(g_rocr$true_score)))
perf_g <- performance(pred_g,"tpr","fpr")
pred_gSA <- prediction(as.numeric(as.character(gSA_rocr$test_score)),
  as.numeric(as.character(gSA_rocr$true_score)))
perf_gSA <- performance(pred_gSA,"tpr","fpr")
pdf("dev/res_balibase_1-5/GUIDANCE_R_vs_SA_balibase_1-5.pdf", width = 6, height = 6)
plot(perf_g, col = "red", main = "Benchmark test \nBALIBASE 1-5, based on SPS & core blocks", lwd = 1.2)
plot(perf_gSA, add = TRUE, col = "blue", lwd = 1.2)
legend("bottomright", c("GUIDANCE R", "GUIDANCE SA"), text.col = c("red", "blue"),
  bty = "n")
dev.off()

## Time elapsed plot
load("dev/res_balibase_1-5/GUIDANCE/time_g.rda")
load("dev/res_balibase_1-5/GUIDANCE/time_gSA.rda")
time <- data.frame(time_g, time_gSA)
se <- function(x) sd(x)/length(x)
df<- cbind.data.frame(
  mean = colMeans(time),
  se = apply(time, 2, se),
  Implementation = c(rep("GUIDANCE R", 5), rep("GUIDANCE SA", 5)))
df <- df[grep("elapsed", rownames(df)),]

p <- ggplot(df, aes(x = Implementation, y = mean, fill = Implementation)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2,
    position = position_dodge(.9))
p <- p + ylab("Time (sek.)")
p <- p  +  theme_bw() + theme(legend.position = "none")
p

pdf("dev/res_balibase_1-5/GUIDANCE_R_vs_SA_balibase_1-5_time.pdf", width = 6, height = 6)
p
dev.off()

## GUIDANCE BENCHMARKING END
#---------------------------------------------------
####################################################
