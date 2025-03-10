library(readr)
library(plyr)
library(dplyr)
library(egg)
library(tidyverse)
library(igraph)
library(future)
library(future.apply)
library(magrittr)
library(readbulk)
library(ggpubr)
library(scales)
library(ggplot2)
library(ggh4x)
library(cowplot)
library(broom)

future::plan(multisession)

setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm/')

###reading 2-, 3-, and 4-way qpAdm models, additional processing and generation of intermediate files for plotting:

files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low-high2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low-high2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low-high2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low-high2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low-high2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low-high2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_standard",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_standard",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_standard",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_standard",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_standard",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_standard",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_standard",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_standard",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_high",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_high",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_high",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_high",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_high",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_high",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_high",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_high")

###read and process *p3* files separately:
#files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2",
#           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2")

###modify the lines below:
#df2 <- read_delim(paste0(f,"_2way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:7,16,19:56,63)]
#df3 <- read_delim(paste0(f,"_3way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:8,26,30:67,76)]
#df4 <- read_delim(paste0(f,"_4way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:9,44,49:86,97)]

for (f in files) {
  df2 <- read_delim(paste0(f,"_2way_p2_correct.txt.gz"), delim="\t", trim_ws=T)[,c(1:7,16,19:56,63)]
  df3 <- read_delim(paste0(f,"_3way_p2_correct.txt.gz"), delim="\t", trim_ws=T)[,c(1:8,26,30:67,76)]
  df4 <- read_delim(paste0(f,"_4way_p2_correct.txt.gz"), delim="\t", trim_ws=T)[,c(1:9,44,49:86,97)]
  
  ###rearranging columns
  df2$model <- "2-way"
  df2 <- df2 %>% relocate(model, .after=target)
  df3$model <- "3-way"
  df3 <- df3 %>% relocate(model, .after=target)
  df4$model <- "4-way"
  df4 <- df4 %>% relocate(model, .after=target)
  
  ###merging all the tables
  dfa <- rbind.fill(df2, df3, df4)
  rm(df2)
  rm(df3)
  rm(df4)
  
  columns <-c("simulation replicate", "max. distance", "avg. distance", "min_S1_T_S2_angle")
  dfa[, columns] <- lapply(columns, function(x) as.numeric(as.character(dfa[[x]])))
  
  write.table(dfa, file=paste0(f,"_for_2Doptimality_metrics_p2_correct.txt"), row.names=F, sep='\t')
}


###reading the intermediate files###
setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm/model_optimality_metrics_2-4way')
df <- read_bulk(directory = "./Isolation_by_distance/stepping-stone_simulations/qpAdm/model_optimality_metrics_2-4way/", name_contains = "_for_2Doptimality_metrics", fun = readr::read_delim, delim="\t")[,c(1:4,6,10:48)]

df<-filter(df, !is.na(min_S1_T_S2_angle))
df$min_S1_T_S2_angle_bin = cut_interval(df$min_S1_T_S2_angle, length=10)

df$min_S1_T_S2_angle_bin <- gsub("[0,10]", "10", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(10,20]", "20", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(20,30]", "30", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(30,40]", "40", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(40,50]", "50", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(50,60]", "60", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(60,70]", "70", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(70,80]", "80", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(80,90]", "90", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(90,100]", "100", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(100,110]", "110", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(110,120]", "120", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(120,130]", "130", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(130,140]", "140", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(140,150]", "150", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(150,160]", "160", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(160,170]", "170", df$min_S1_T_S2_angle_bin, fixed=T)
df$min_S1_T_S2_angle_bin <- gsub("(170,180]", "180", df$min_S1_T_S2_angle_bin, fixed=T)

df <- df %>% relocate(min_S1_T_S2_angle, .after=avg..distance)
df <- df %>% relocate(min_S1_T_S2_angle_bin, .after=min_S1_T_S2_angle)
columns <-c("min_S1_T_S2_angle_bin")
df[, columns] <- lapply(columns, function(x) as.numeric(df[[x]]))
df$avg..distance <- round_any(df$avg..distance, 0.5, f = round) #!!!!!!!

dfr1 <- filter(df, qpAdm.setup=="dist. rot.")
dfr2 <- filter(df, qpAdm.setup=="prox. rot.")
dfnr1 <- filter(df, qpAdm.setup=="dist. non-rot.")
dfnr2 <- filter(df, qpAdm.setup=="prox. non-rot.")
dfr <- rbind(dfr1, dfr2)
dfnr <- rbind(dfnr1, dfnr2)
rm(dfr1)
rm(dfr2)
rm(dfnr1)
rm(dfnr2)


###PLOTTING FIGURES###

###the results are stratified by model complexity and "sparsity" of stepping-stone landscapes
###("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2": these are approximate intervals for random sampling of gene-flow intensities from uniform or normal distributions)
###the results are shown for seven selected composite model feasibility criteria, with all the qpAdm protocols combined:

###counting fitting and non-fitting models in 2D bins###
fit_criteria <- read_delim("fit_criteria2.txt", delim="\t", trim_ws=T)
#crs = c(32,28,24,15,9,3,7) #v. 1
crs = c(32,28,22,15,9,3,7) #v. 2

###model feasibility criterion no. 1###
cr <- crs[1]
dfa <- df[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df01 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df01) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df01 <- df01[with(df01, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df01f <- filter(df01, positive=="fitting")
df01nf <- filter(df01, positive=="non-fitting")
df01f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01fin <- cbind(df01f[,1:4], df01f$model_count/(df01f$model_count+df01nf$model_count))
colnames(df01fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[df01fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df01fin[[x]])))
df01fin <- df01fin[with(df01fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df01fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts01.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts01.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts01.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts01.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df01)

####model feasibility criterion no. 2###
cr <- crs[2]
dfa <- df[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df02 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df02) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df02 <- df02[with(df02, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df02f <- filter(df02, positive=="fitting")
df02nf <- filter(df02, positive=="non-fitting")
df02f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02fin <- cbind(df02f[,1:4], df02f$model_count/(df02f$model_count+df02nf$model_count))
colnames(df02fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[df02fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df02fin[[x]])))
df02fin <- df02fin[with(df02fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df02fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts02.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts02.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts02.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts02.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df02)

####model feasibility criterion no. 3###
cr <- crs[3]
dfa <- df[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df03 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df03) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df03 <- df03[with(df03, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df03f <- filter(df03, positive=="fitting")
df03nf <- filter(df03, positive=="non-fitting")
df03f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03fin <- cbind(df03f[,1:4], df03f$model_count/(df03f$model_count+df03nf$model_count))
colnames(df03fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[df03fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df03fin[[x]])))
df03fin <- df03fin[with(df03fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df03fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts03.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts03.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts03.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts03.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df03)

####model feasibility criterion no. 4###
cr <- crs[4]
dfa <- df[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df04 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df04) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df04 <- df04[with(df04, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df04f <- filter(df04, positive=="fitting")
df04nf <- filter(df04, positive=="non-fitting")
df04f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04fin <- cbind(df04f[,1:4], df04f$model_count/(df04f$model_count+df04nf$model_count))
colnames(df04fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[df04fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df04fin[[x]])))
df04fin <- df04fin[with(df04fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df04fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts04.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts04.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts04.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts04.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df04)

####model feasibility criterion no. 5###
cr <- crs[5]
dfa <- df[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df05 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df05) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df05 <- df05[with(df05, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df05f <- filter(df05, positive=="fitting")
df05nf <- filter(df05, positive=="non-fitting")
df05f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05fin <- cbind(df05f[,1:4], df05f$model_count/(df05f$model_count+df05nf$model_count))
colnames(df05fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[df05fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df05fin[[x]])))
df05fin <- df05fin[with(df05fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df05fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts05.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts05.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts05.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts05.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df05)

####model feasibility criterion no. 6###
cr <- crs[6]
dfa <- df[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df06 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df06) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df06 <- df06[with(df06, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df06f <- filter(df06, positive=="fitting")
df06nf <- filter(df06, positive=="non-fitting")
df06f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06fin <- cbind(df06f[,1:4], df06f$model_count/(df06f$model_count+df06nf$model_count))
colnames(df06fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[df06fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df06fin[[x]])))
df06fin <- df06fin[with(df06fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df06fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts06.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts06.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts06.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts06.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df06)

####model feasibility criterion no. 7###
cr <- crs[7]
dfa <- df[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df07 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df07) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df07 <- df07[with(df07, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df07f <- filter(df07, positive=="fitting")
df07nf <- filter(df07, positive=="non-fitting")
df07f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df07nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df07fin <- cbind(df07f[,1:4], df07f$model_count/(df07f$model_count+df07nf$model_count))
colnames(df07fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df07fin[df07fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df07fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df07fin[[x]])))
df07fin <- df07fin[with(df07fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df07fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts07.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts07.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts07.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts07.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts07.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts07.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts07.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts07.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts07.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts07.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df07)

dff <- rbind(df01f, df02f, df03f, df04f, df05f, df06f, df07f)
dfnf <- rbind(df01nf, df02nf, df03nf, df04nf, df05nf, df06nf, df07nf)
dffin <- rbind(df01fin, df02fin, df03fin, df04fin, df05fin, df06fin, df07fin)
counts.f <- rbind(counts01.f, counts02.f, counts03.f, counts04.f, counts05.f, counts06.f, counts07.f)
counts.nf <- rbind(counts01.nf, counts02.nf, counts03.nf, counts04.nf, counts05.nf, counts06.nf, counts07.nf)
rm(df01f)
rm(df02f)
rm(df03f)
rm(df04f)
rm(df05f)
rm(df06f)
rm(df07f)
rm(df01nf)
rm(df02nf)
rm(df03nf)
rm(df04nf)
rm(df05nf)
rm(df06nf)
rm(df07nf)
rm(df01fin)
rm(df02fin)
rm(df03fin)
rm(df04fin)
rm(df05fin)
rm(df06fin)
rm(df07fin)
rm(counts01.f)
rm(counts02.f)
rm(counts03.f)
rm(counts04.f)
rm(counts05.f)
rm(counts06.f)
rm(counts07.f)
rm(counts01.nf)
rm(counts02.nf)
rm(counts03.nf)
rm(counts04.nf)
rm(counts05.nf)
rm(counts06.nf)
rm(counts07.nf)

dffin$nf <- dff$model_count
dffin$nt <- dff$model_count+dfnf$model_count
dffin$nf <- as.double(dffin$nf)
dffin$nt <- as.double(dffin$nt)
dffin$nt[dffin$nt==0] <- NA ##IMPORTANT!
corr <- dffin  %>% group_by(gene.flow.intensity.per.generation, model, fit_criteria1, fit_criteria2, fit_criteria3)  %>% do(tidy(cor.test(.$nt, .$nf, alternative="two.sided", method = "spearman")))
#counts.f$PR_global <- counts.f$model_count/(counts.f$model_count+counts.nf$model_count)
#counts.f$n_total <- counts.f$model_count+counts.nf$model_count
#counts.f$n_total[13:84] <- NA

###plotting:
p00<-ggplot(df, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_bin2d(aes(fill = after_stat(count)), binwidth = c(1, 10), origin = c(-0.5, -5))+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  ggtitle("model count, logarithmic color scale")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,180), expand=T, default=F, clip="on")+
  ylab("")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=32, face='bold'),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'grey95'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(10,110,10,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p01<-ggplot(dffin, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_tile(aes(fill = nf))+
  geom_text(data=counts.f, aes(x=2.5, y=200, label=scales::comma(model_count)), colour="darkgreen", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=counts.nf, aes(x=7, y=200, label=scales::comma(model_count)), colour="red3", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=corr, aes(x=9, y=150, label=round(estimate, 2)), colour="tomato4", size=8, fontface="bold", check_overlap=T, angle=90)+
  geom_hline(aes(yintercept=125), colour="red3", linewidth=1, alpha=0.5)+
  geom_hline(aes(yintercept=95), colour="red3", linewidth=1, alpha=0.5)+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,200), expand=T, default=F, clip="on")+
  ggtitle('counts of feasible models for 2-4-way models (23,439,127 in total) in the space of optimality metrics,\nlogarithmic color scale, progressively more stringent feasibility criteria')+
  xlab("max. source-target distance")+
  ylab("min. S1-target-S2 angle")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=24))+
  theme(plot.title = element_text(color="darkmagenta", face="bold", size=32),
        axis.text.x=element_text(size=24, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_text(size=36, face='bold'),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,10,5,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(fit_criteria3 + factor(fit_criteria1, level=c("% (-0.3,1.3)", "% (0,1)", "%+/-2SE (0,1)")) + fit_criteria2 ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

pall <- plot_grid(
  plot_grid(p00+theme(legend.position="none"), p01+theme(legend.position="none"), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  plot_grid(get_legend(p00), get_legend(p01), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  rel_widths = c(33,3))

###This is Figure 6###
pdf("distributions_maxdist_angles_2D_counts_2-4way_simplified_p3.pdf", width=36, height=28)
print(pall)
dev.off()




###the results are stratified by qpAdm protocols, by model complexity and "sparsity" of stepping-stone landscapes
###("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2": these are approximate intervals for random sampling of gene-flow intensities from uniform or normal distributions)
###the results are shown for a selected composite model feasibility criterion no. 22:

###counting fitting and non-fitting models in 2D bins###
fit_criteria <- read_delim("fit_criteria2.txt", delim="\t", trim_ws=T)
cr <- 22
dfa <- df[,c(1:6,9,9+cr)]
dfa$qpAdm.setup <- gsub("dist. rot.", "1.dist. rot.", dfa$qpAdm.setup)
dfa$qpAdm.setup <- gsub("prox. rot.", "2.prox. rot.", dfa$qpAdm.setup)
dfa$qpAdm.setup <- gsub("dist. non-rot.", "3.dist. non-rot.", dfa$qpAdm.setup)
dfa$qpAdm.setup <- gsub("prox. non-rot.", "4.prox. non-rot.", dfa$qpAdm.setup)
colnames(dfa) <- c("pop_no", "qpAdm.setup", "gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df01 <- as.data.frame(xtabs(~ pop_no+qpAdm.setup+gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df01) <- c("pop_no", "qpAdm.setup", "gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df01 <- df01[with(df01, order(pop_no, qpAdm.setup, gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df01f <- filter(df01, positive=="fitting")
df01nf <- filter(df01, positive=="non-fitting")
df01fin <- cbind(df01f[,1:6], df01f$model_count/(df01f$model_count+df01nf$model_count))
colnames(df01fin) <- c("pop_no", "qpAdm.setup", "gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[df01fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df01fin[[x]])))
df01fin <- df01fin[with(df01fin, order(pop_no, qpAdm.setup, gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]

counts.f <- as.data.frame(xtabs(~ pop_no+qpAdm.setup+gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts.nf <- as.data.frame(xtabs(~ pop_no+qpAdm.setup+gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts.f) <- c("pop_no", "qpAdm.setup", "gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts.nf) <- c("pop_no", "qpAdm.setup", "gene.flow.intensity.per.generation", "model", "model_count")
counts.f$PR_global <- counts.f$model_count/(counts.f$model_count+counts.nf$model_count)
counts.f$n_total <- counts.f$model_count+counts.nf$model_count
#counts.f$n_total[13:84] <- NA

df01fin$nf <- df01f$model_count
df01fin$nt <- df01f$model_count+df01nf$model_count
df01fin$nf <- as.double(df01fin$nf)
df01fin$nt <- as.double(df01fin$nt)
df01fin$nt[df01fin$nt==0] <- NA ##IMPORTANT!
corr1 <- df01fin  %>% group_by(pop_no, qpAdm.setup, gene.flow.intensity.per.generation, model)  %>% do(tidy(cor.test(.$nt, .$PR, alternative="two.sided", method = "spearman"))) 
corr2 <- df01fin  %>% group_by(pop_no, qpAdm.setup, gene.flow.intensity.per.generation, model)  %>% do(tidy(cor.test(.$nt, .$nf, alternative="two.sided", method = "spearman"))) 


###plotting:
p00<-ggplot(df01fin, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_tile(aes(fill = nt))+
  geom_text(data=counts.f, aes(x=2.5, y=200, label=scales::comma(model_count)), colour="darkgreen", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=counts.nf, aes(x=7, y=200, label=scales::comma(model_count)), colour="red3", size=7, fontface="bold", check_overlap=T)+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  ggtitle(paste0(fit_criteria[cr,2], "; ", fit_criteria[cr,3], "; ", fit_criteria[cr,4]))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,200), expand=T, default=F, clip="on")+
  xlab("max. source-target distance")+
  ylab("min. S1-target-S2 angle")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=24))+
  theme(plot.title = element_text(size=32, face='bold'),
        axis.text.x=element_text(size=24, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_text(size=36, face='bold'),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'grey95'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,10,5,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(pop_no + qpAdm.setup ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p01<-ggplot(df01fin, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_tile(aes(fill = PR))+
  geom_text(data=counts.f, aes(x=2.5, y=200, label=scales::scientific(PR_global, digits=2)), colour="darkgreen", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=counts.f, aes(x=7, y=200, label=scales::comma(n_total)), colour="navy", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=corr1, aes(x=9, y=150, label=round(estimate, 2)), colour="tomato4", size=8, fontface="bold", check_overlap=T, angle=90)+
  scale_fill_gradientn(colours = pals::parula(20), trans = "sqrt", breaks = c(0,0.05,0.15,0.3,0.6,0.9), labels = c(0,0.05,0.15,0.3,0.6,0.9), limits = c(0,1), na.value='grey95')+
  ggtitle(paste0(fit_criteria[cr,2], "; ", fit_criteria[cr,3], "; ", fit_criteria[cr,4]))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,200), expand=T, default=F, clip="on")+
  xlab("max. source-target distance")+
  ylab("min. S1-target-S2 angle")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=24))+
  theme(plot.title = element_text(size=32, face='bold'),
        axis.text.x=element_text(size=24, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_text(size=36, face='bold'),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,10,5,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(pop_no + qpAdm.setup ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p02<-ggplot(df01fin, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_tile(aes(fill = nf))+
  geom_text(data=counts.f, aes(x=2.5, y=200, label=scales::scientific(PR_global, digits=2)), colour="darkgreen", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=counts.f, aes(x=7, y=200, label=scales::comma(n_total)), colour="navy", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=corr2, aes(x=9, y=150, label=round(estimate, 2)), colour="tomato4", size=8, fontface="bold", check_overlap=T, angle=90)+
  geom_hline(aes(yintercept=125), colour="red3", linewidth=1, alpha=0.5)+
  geom_hline(aes(yintercept=95), colour="red3", linewidth=1, alpha=0.5)+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  ggtitle(paste0(fit_criteria[cr,2], "; ", fit_criteria[cr,3], "; ", fit_criteria[cr,4]))+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,200), expand=T, default=F, clip="on")+
  xlab("max. source-target distance")+
  ylab("min. S1-target-S2 angle")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=24))+
  theme(plot.title = element_text(size=32, face='bold'),
        axis.text.x=element_text(size=24, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_text(size=36, face='bold'),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,10,5,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(pop_no + qpAdm.setup ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

pall00 <- plot_grid(p00+theme(legend.position="none"), get_legend(p00), rel_widths = c(33,3))

pall01 <- plot_grid(p01+theme(legend.position="none"), get_legend(p01), rel_widths = c(33,3))

pall02 <- plot_grid(p02+theme(legend.position="none"), get_legend(p02), rel_widths = c(33,3))

###This is Figure S25a###
pdf(paste0("distributions_maxdist_angles_2D_density_2-4way_p",cr,".pdf"), width=36, height=27)
print(annotate_figure(pall00, top=text_grob('density of 2-4-way models tested (23,439,127 in total) in the space of optimality metrics, logarithmic color scale\n', color="darkmagenta", face="bold", size=32)))
dev.off()

###This is Figure S25b###
pdf(paste0("distributions_maxdist_angles_2D_PR_2-4way_p",cr,".pdf"), width=36, height=27)
print(annotate_figure(pall01, top=text_grob('positive rates for 2-4-way models (23,439,127 in total) in the space of optimality metrics, square root color scale\n', color="darkmagenta", face="bold", size=32)))
dev.off()

###This is Figure S25c###
pdf(paste0("distributions_maxdist_angles_2D_counts_2-4way_p",cr,".pdf"), width=36, height=27)
print(annotate_figure(pall02, top=text_grob('counts of feasible models for 2-4-way models (23,439,127 in total) in the space of optimality metrics, logarithmic color scale\n', color="darkmagenta", face="bold", size=32)))
dev.off()




###the results are stratified by model complexity and "sparsity" of stepping-stone landscapes
###("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2": these are approximate intervals for random sampling of gene-flow intensities from uniform or normal distributions)
###the results are shown for seven selected composite model feasibility criteria, with all rotating or all non-rotating qpAdm protocols combined:

###counting fitting and non-fitting models in 2D bins###
###model feasibility criterion no. 1###
fit_criteria <- read_delim("fit_criteria2.txt", delim="\t", trim_ws=T)
crs = c(32,28,22,15,9,3,7) #v. 2
cr <- crs[1]
dfa <- dfr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df01 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df01) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df01 <- df01[with(df01, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df01f <- filter(df01, positive=="fitting")
df01nf <- filter(df01, positive=="non-fitting")
df01f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01fin <- cbind(df01f[,1:4], df01f$model_count/(df01f$model_count+df01nf$model_count))
colnames(df01fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[df01fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df01fin[[x]])))
df01fin <- df01fin[with(df01fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df01fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts01.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts01.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts01.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts01.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df01)

###model feasibility criterion no. 2###
cr <- crs[2]
dfa <- dfr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df02 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df02) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df02 <- df02[with(df02, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df02f <- filter(df02, positive=="fitting")
df02nf <- filter(df02, positive=="non-fitting")
df02f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02fin <- cbind(df02f[,1:4], df02f$model_count/(df02f$model_count+df02nf$model_count))
colnames(df02fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[df02fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df02fin[[x]])))
df02fin <- df02fin[with(df02fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df02fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts02.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts02.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts02.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts02.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df02)

###model feasibility criterion no. 3###
cr <- crs[3]
dfa <- dfr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df03 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df03) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df03 <- df03[with(df03, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df03f <- filter(df03, positive=="fitting")
df03nf <- filter(df03, positive=="non-fitting")
df03f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03fin <- cbind(df03f[,1:4], df03f$model_count/(df03f$model_count+df03nf$model_count))
colnames(df03fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[df03fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df03fin[[x]])))
df03fin <- df03fin[with(df03fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df03fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts03.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts03.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts03.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts03.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df03)

###model feasibility criterion no. 4###
cr <- crs[4]
dfa <- dfr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df04 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df04) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df04 <- df04[with(df04, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df04f <- filter(df04, positive=="fitting")
df04nf <- filter(df04, positive=="non-fitting")
df04f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04fin <- cbind(df04f[,1:4], df04f$model_count/(df04f$model_count+df04nf$model_count))
colnames(df04fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[df04fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df04fin[[x]])))
df04fin <- df04fin[with(df04fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df04fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts04.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts04.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts04.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts04.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df04)

###model feasibility criterion no. 5###
cr <- crs[5]
dfa <- dfr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df05 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df05) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df05 <- df05[with(df05, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df05f <- filter(df05, positive=="fitting")
df05nf <- filter(df05, positive=="non-fitting")
df05f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05fin <- cbind(df05f[,1:4], df05f$model_count/(df05f$model_count+df05nf$model_count))
colnames(df05fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[df05fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df05fin[[x]])))
df05fin <- df05fin[with(df05fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df05fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts05.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts05.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts05.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts05.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df05)

###model feasibility criterion no. 6###
cr <- crs[6]
dfa <- dfr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df06 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df06) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df06 <- df06[with(df06, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df06f <- filter(df06, positive=="fitting")
df06nf <- filter(df06, positive=="non-fitting")
df06f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06fin <- cbind(df06f[,1:4], df06f$model_count/(df06f$model_count+df06nf$model_count))
colnames(df06fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[df06fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df06fin[[x]])))
df06fin <- df06fin[with(df06fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df06fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts06.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts06.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts06.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts06.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df06)

###model feasibility criterion no. 7###
cr <- crs[7]
dfa <- dfr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df07 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df07) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df07 <- df07[with(df07, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df07f <- filter(df07, positive=="fitting")
df07nf <- filter(df07, positive=="non-fitting")
df07f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df07nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df07fin <- cbind(df07f[,1:4], df07f$model_count/(df07f$model_count+df07nf$model_count))
colnames(df07fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df07fin[df07fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df07fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df07fin[[x]])))
df07fin <- df07fin[with(df07fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df07fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts07.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts07.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts07.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts07.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts07.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts07.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts07.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts07.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts07.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts07.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df07)

dff <- rbind(df01f, df02f, df03f, df04f, df05f, df06f, df07f)
dfnf <- rbind(df01nf, df02nf, df03nf, df04nf, df05nf, df06nf, df07nf)
dffin <- rbind(df01fin, df02fin, df03fin, df04fin, df05fin, df06fin, df07fin)
counts.f <- rbind(counts01.f, counts02.f, counts03.f, counts04.f, counts05.f, counts06.f, counts07.f)
counts.nf <- rbind(counts01.nf, counts02.nf, counts03.nf, counts04.nf, counts05.nf, counts06.nf, counts07.nf)
rm(df01f)
rm(df02f)
rm(df03f)
rm(df04f)
rm(df05f)
rm(df06f)
rm(df07f)
rm(df01nf)
rm(df02nf)
rm(df03nf)
rm(df04nf)
rm(df05nf)
rm(df06nf)
rm(df07nf)
rm(df01fin)
rm(df02fin)
rm(df03fin)
rm(df04fin)
rm(df05fin)
rm(df06fin)
rm(df07fin)
rm(counts01.f)
rm(counts02.f)
rm(counts03.f)
rm(counts04.f)
rm(counts05.f)
rm(counts06.f)
rm(counts07.f)
rm(counts01.nf)
rm(counts02.nf)
rm(counts03.nf)
rm(counts04.nf)
rm(counts05.nf)
rm(counts06.nf)
rm(counts07.nf)

dffin$nf <- dff$model_count
dffin$nt <- dff$model_count+dfnf$model_count
dffin$nf <- as.double(dffin$nf)
dffin$nt <- as.double(dffin$nt)
dffin$nt[dffin$nt==0] <- NA ##IMPORTANT!
corr <- dffin  %>% group_by(gene.flow.intensity.per.generation, model, fit_criteria1, fit_criteria2, fit_criteria3)  %>% do(tidy(cor.test(.$nt, .$nf, alternative="two.sided", method = "spearman"))) 
#counts.f$PR_global <- counts.f$model_count/(counts.f$model_count+counts.nf$model_count)
#counts.f$n_total <- counts.f$model_count+counts.nf$model_count
#counts.f$n_total[13:84] <- NA

###plotting:
p00<-ggplot(dfr, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_bin2d(aes(fill = after_stat(count)), binwidth = c(1, 10), origin = c(-0.5, -5))+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  ggtitle("model count, logarithmic color scale")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,180), expand=T, default=F, clip="on")+
  ylab("")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=32, face='bold'),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'grey95'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(10,110,10,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p01<-ggplot(dffin, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_tile(aes(fill = nf))+
  geom_text(data=counts.f, aes(x=2.5, y=200, label=scales::comma(model_count)), colour="darkgreen", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=counts.nf, aes(x=7, y=200, label=scales::comma(model_count)), colour="red3", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=corr, aes(x=9, y=150, label=round(estimate, 2)), colour="tomato4", size=8, fontface="bold", check_overlap=T, angle=90)+
  geom_hline(aes(yintercept=125), colour="red3", linewidth=1, alpha=0.5)+
  geom_hline(aes(yintercept=95), colour="red3", linewidth=1, alpha=0.5)+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,200), expand=T, default=F, clip="on")+
  ggtitle('counts of feasible models for 2-4-way "rotating" models (16,390,119 in total) in the space of optimality metrics,\nlogarithmic color scale, progressively more stringent feasibility criteria')+
  xlab("max. source-target distance")+
  ylab("min. S1-target-S2 angle")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=24))+
  theme(plot.title = element_text(color="darkmagenta", face="bold", size=32),
        axis.text.x=element_text(size=24, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_text(size=36, face='bold'),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,10,5,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(fit_criteria3 + factor(fit_criteria1, level=c("% (-0.3,1.3)", "% (0,1)", "%+/-2SE (0,1)")) + fit_criteria2 ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

pall <- plot_grid(
  plot_grid(p00+theme(legend.position="none"), p01+theme(legend.position="none"), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  plot_grid(get_legend(p00), get_legend(p01), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  rel_widths = c(33,3))

###This is Figure S18c###
pdf("distributions_maxdist_angles_2D_counts_2-4way_simplified_p3_rot.pdf", width=36, height=28)
print(pall)
dev.off()


###model feasibility criterion no. 1###
cr <- crs[1]
dfa <- dfnr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df01 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df01) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df01 <- df01[with(df01, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df01f <- filter(df01, positive=="fitting")
df01nf <- filter(df01, positive=="non-fitting")
df01f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df01fin <- cbind(df01f[,1:4], df01f$model_count/(df01f$model_count+df01nf$model_count))
colnames(df01fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[df01fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df01fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df01fin[[x]])))
df01fin <- df01fin[with(df01fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df01fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df01fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df01fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts01.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts01.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts01.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts01.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts01.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts01.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts01.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df01)

###model feasibility criterion no. 2###
cr <- crs[2]
dfa <- dfnr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df02 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df02) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df02 <- df02[with(df02, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df02f <- filter(df02, positive=="fitting")
df02nf <- filter(df02, positive=="non-fitting")
df02f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df02fin <- cbind(df02f[,1:4], df02f$model_count/(df02f$model_count+df02nf$model_count))
colnames(df02fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[df02fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df02fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df02fin[[x]])))
df02fin <- df02fin[with(df02fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df02fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df02fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df02fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts02.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts02.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts02.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts02.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts02.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts02.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts02.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df02)

###model feasibility criterion no. 3###
cr <- crs[3]
dfa <- dfnr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df03 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df03) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df03 <- df03[with(df03, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df03f <- filter(df03, positive=="fitting")
df03nf <- filter(df03, positive=="non-fitting")
df03f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df03fin <- cbind(df03f[,1:4], df03f$model_count/(df03f$model_count+df03nf$model_count))
colnames(df03fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[df03fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df03fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df03fin[[x]])))
df03fin <- df03fin[with(df03fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df03fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df03fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df03fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts03.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts03.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts03.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts03.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts03.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts03.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts03.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df03)

###model feasibility criterion no. 4###
cr <- crs[4]
dfa <- dfnr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df04 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df04) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df04 <- df04[with(df04, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df04f <- filter(df04, positive=="fitting")
df04nf <- filter(df04, positive=="non-fitting")
df04f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df04fin <- cbind(df04f[,1:4], df04f$model_count/(df04f$model_count+df04nf$model_count))
colnames(df04fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[df04fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df04fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df04fin[[x]])))
df04fin <- df04fin[with(df04fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df04fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df04fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df04fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts04.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts04.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts04.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts04.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts04.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts04.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts04.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df04)

###model feasibility criterion no. 5###
cr <- crs[5]
dfa <- dfnr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df05 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df05) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df05 <- df05[with(df05, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df05f <- filter(df05, positive=="fitting")
df05nf <- filter(df05, positive=="non-fitting")
df05f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df05fin <- cbind(df05f[,1:4], df05f$model_count/(df05f$model_count+df05nf$model_count))
colnames(df05fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[df05fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df05fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df05fin[[x]])))
df05fin <- df05fin[with(df05fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df05fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df05fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df05fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts05.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts05.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts05.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts05.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts05.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts05.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts05.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df05)

###model feasibility criterion no. 6###
cr <- crs[6]
dfa <- dfnr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df06 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df06) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df06 <- df06[with(df06, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df06f <- filter(df06, positive=="fitting")
df06nf <- filter(df06, positive=="non-fitting")
df06f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df06fin <- cbind(df06f[,1:4], df06f$model_count/(df06f$model_count+df06nf$model_count))
colnames(df06fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[df06fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df06fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df06fin[[x]])))
df06fin <- df06fin[with(df06fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df06fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df06fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df06fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts06.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts06.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts06.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts06.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts06.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts06.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts06.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df06)

###model feasibility criterion no. 7###
cr <- crs[7]
dfa <- dfnr[,c(3:6,9,9+cr)]
colnames(dfa) <- c("gene.flow.intensity.per.generation", "simulation.replicate", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive")
df07 <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model+max..distance+min_S1_T_S2_angle_bin+positive, data = dfa))
colnames(df07) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "positive", "model_count")
df07 <- df07[with(df07, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin, positive)), ]
df07f <- filter(df07, positive=="fitting")
df07nf <- filter(df07, positive=="non-fitting")
df07f$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07f$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07f$fit_criteria3 <- as.character(fit_criteria[cr,4])
df07nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
df07fin <- cbind(df07f[,1:4], df07f$model_count/(df07f$model_count+df07nf$model_count))
colnames(df07fin) <- c("gene.flow.intensity.per.generation", "model", "max..distance", "min_S1_T_S2_angle_bin", "PR")
df07fin[df07fin == "NaN"] <- NA
columns <-c("max..distance", "min_S1_T_S2_angle_bin", "PR")
df07fin[, columns] <- lapply(columns, function(x) as.numeric(as.character(df07fin[[x]])))
df07fin <- df07fin[with(df07fin, order(gene.flow.intensity.per.generation, model, max..distance, min_S1_T_S2_angle_bin)), ]
df07fin$fit_criteria1 <- as.character(fit_criteria[cr,2])
df07fin$fit_criteria2 <- as.character(fit_criteria[cr,3])
df07fin$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts07.f <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="fitting")))
counts07.nf <- as.data.frame(xtabs(~ gene.flow.intensity.per.generation+model, data = filter(dfa, positive=="non-fitting")))
colnames(counts07.f) <- c("gene.flow.intensity.per.generation", "model", "model_count")
colnames(counts07.nf) <- c("gene.flow.intensity.per.generation", "model", "model_count")
counts07.f$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts07.f$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts07.f$fit_criteria3 <- as.character(fit_criteria[cr,4])
counts07.nf$fit_criteria1 <- as.character(fit_criteria[cr,2])
counts07.nf$fit_criteria2 <- as.character(fit_criteria[cr,3])
counts07.nf$fit_criteria3 <- as.character(fit_criteria[cr,4])
rm(dfa)
rm(df07)

dff <- rbind(df01f, df02f, df03f, df04f, df05f, df06f, df07f)
dfnf <- rbind(df01nf, df02nf, df03nf, df04nf, df05nf, df06nf, df07nf)
dffin <- rbind(df01fin, df02fin, df03fin, df04fin, df05fin, df06fin, df07fin)
counts.f <- rbind(counts01.f, counts02.f, counts03.f, counts04.f, counts05.f, counts06.f, counts07.f)
counts.nf <- rbind(counts01.nf, counts02.nf, counts03.nf, counts04.nf, counts05.nf, counts06.nf, counts07.nf)
rm(df01f)
rm(df02f)
rm(df03f)
rm(df04f)
rm(df05f)
rm(df06f)
rm(df07f)
rm(df01nf)
rm(df02nf)
rm(df03nf)
rm(df04nf)
rm(df05nf)
rm(df06nf)
rm(df07nf)
rm(df01fin)
rm(df02fin)
rm(df03fin)
rm(df04fin)
rm(df05fin)
rm(df06fin)
rm(df07fin)
rm(counts01.f)
rm(counts02.f)
rm(counts03.f)
rm(counts04.f)
rm(counts05.f)
rm(counts06.f)
rm(counts07.f)
rm(counts01.nf)
rm(counts02.nf)
rm(counts03.nf)
rm(counts04.nf)
rm(counts05.nf)
rm(counts06.nf)
rm(counts07.nf)

dffin$nf <- dff$model_count
dffin$nt <- dff$model_count+dfnf$model_count
dffin$nf <- as.double(dffin$nf)
dffin$nt <- as.double(dffin$nt)
dffin$nt[dffin$nt==0] <- NA ##IMPORTANT!
corr <- dffin  %>% group_by(gene.flow.intensity.per.generation, model, fit_criteria1, fit_criteria2, fit_criteria3)  %>% do(tidy(cor.test(.$nt, .$nf, alternative="two.sided", method = "spearman"))) 
#counts.f$PR_global <- counts.f$model_count/(counts.f$model_count+counts.nf$model_count)
#counts.f$n_total <- counts.f$model_count+counts.nf$model_count
#counts.f$n_total[13:84] <- NA

###plotting:
p00<-ggplot(dfnr, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_bin2d(aes(fill = after_stat(count)), binwidth = c(1, 10), origin = c(-0.5, -5))+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  ggtitle("model count, logarithmic color scale")+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,180), expand=T, default=F, clip="on")+
  ylab("")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=32, face='bold'),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'grey95'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(10,110,10,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p01<-ggplot(dffin, aes(x=max..distance, y=min_S1_T_S2_angle_bin))+
  geom_tile(aes(fill = nf))+
  geom_text(data=counts.f, aes(x=2.5, y=200, label=scales::comma(model_count)), colour="darkgreen", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=counts.nf, aes(x=7, y=200, label=scales::comma(model_count)), colour="red3", size=7, fontface="bold", check_overlap=T)+
  geom_text(data=corr, aes(x=9, y=150, label=round(estimate, 2)), colour="tomato4", size=8, fontface="bold", check_overlap=T, angle=90)+
  geom_hline(aes(yintercept=125), colour="red3", linewidth=1, alpha=0.5)+
  geom_hline(aes(yintercept=95), colour="red3", linewidth=1, alpha=0.5)+
  scale_fill_gradientn(colours = pals::parula(20), trans = "log", breaks = c(1,10,100,1000,10000,100000), labels = c(1,10,100,1000,10000,100000), na.value='grey95')+
  scale_x_continuous(breaks=c(1,2,3,4,5,6,7,8,9), expand=c(0.05,0.05))+
  scale_y_continuous(breaks=c(30,60,90,120,150,180), expand=c(0.05,0.05))+
  coord_cartesian(xlim=c(1,9), ylim=c(10,200), expand=T, default=F, clip="on")+
  ggtitle('counts of feasible models for 2-4-way "non-rotating" models (7,049,008 in total) in the space of optimality metrics,\nlogarithmic color scale, progressively more stringent feasibility criteria')+
  xlab("max. source-target distance")+
  ylab("min. S1-target-S2 angle")+
  guides(fill = guide_colorbar(title = '', order=1, barwidth=4, barheight=24))+
  theme(plot.title = element_text(color="darkmagenta", face="bold", size=32),
        axis.text.x=element_text(size=24, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=24, face='bold'),
        axis.title.x=element_text(size=36, face='bold'),
        axis.title.y=element_text(size=36, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=30),
        panel.background=element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5,10,5,10), "pt"),
        strip.text=element_text(size=26, face='bold'))+
  facet_nested(fit_criteria3 + factor(fit_criteria1, level=c("% (-0.3,1.3)", "% (0,1)", "%+/-2SE (0,1)")) + fit_criteria2 ~ model + factor(gene.flow.intensity.per.generation, level=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

pall <- plot_grid(
  plot_grid(p00+theme(legend.position="none"), p01+theme(legend.position="none"), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  plot_grid(get_legend(p00), get_legend(p01), ncol = 1, axis="l", align = "h", rel_heights = c(1.2,6.8)),
  rel_widths = c(33,3))

###This is Figure S18e###
pdf("distributions_maxdist_angles_2D_counts_2-4way_simplified_p3_non-rot.pdf", width=36, height=28)
print(pall)
dev.off()

