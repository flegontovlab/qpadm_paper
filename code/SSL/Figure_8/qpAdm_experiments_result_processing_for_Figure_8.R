library(readr)
library(plyr)
library(dplyr)
library(egg)
library(tidyverse)
library(igraph)
library(future)
library(future.apply)
library(magrittr)
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
#df2 <- read_delim(paste0(f,"_2way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:7,13:56,64,65)]
#df3 <- read_delim(paste0(f,"_3way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:8,26,30:67,77,78)]
#df4 <- read_delim(paste0(f,"_4way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:9,44,49:86,98,99)]


for (f in files) {
df2 <- read_delim(paste0(f,"_2way_p2_correct.txt.gz"), delim="\t", trim_ws=T)[,c(1:7,13:56,64,65)]
df3 <- read_delim(paste0(f,"_3way_p2_correct.txt.gz"), delim="\t", trim_ws=T)[,c(1:8,26,30:67,77,78)]
df4 <- read_delim(paste0(f,"_4way_p2_correct.txt.gz"), delim="\t", trim_ws=T)[,c(1:9,44,49:86,98,99)]

###adding s1-s2 distances###
g = read.table('filledcircle_graph.txt', header = T) %>% 
  transmute(from = pop1, to = pop2) %>% 
  igraph::graph_from_data_frame()

s1 = as.data.frame(t(as.data.frame(strsplit(df2$s1, split="_"))))
s2 = as.data.frame(t(as.data.frame(strsplit(df2$s2, split="_"))))
models <- as.data.frame(cbind(s1$V1, s2$V1))
colnames(models) <- c("s1", "s2")

res1 = future_apply(models, 1, function(x) {
  out<-shortest_paths(g, from=x[1], to=x[2], mode="all", output="vpath", algorithm="unweighted")$vpath
})

combs <- c(1:length(models$s1))
outtab2 = sapply(combs, function(comb){
  s1_s2_dist = length(res1[[comb]][[1]])-1
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab2) <- c("distance_s1.s2")

df2 <- cbind(df2, outtab2)
df2 <- df2 %>% relocate(distance_s1.s2, .before=`max. distance`)
rm(combs)
rm(s1)
rm(s2)
rm(res1)
rm(models)
rm(outtab2)

###generating a 1-way table
df5a <- df2[,c(1:4,5,6,8,11,12)]
df5b <- df2[,c(1:4,5,7,9,11,13)]
df5c <- df2[,c(1:4,6,7,10,11,14)] #s1*s2 pairs
colnames(df5a) <- c("pop_no", "qpAdm setup", "gene flow intensity per generation", "simulation replicate", "s1", "s2", "p_1way_s1.s2", "right", "distance_t.s1")
colnames(df5b) <- c("pop_no", "qpAdm setup", "gene flow intensity per generation", "simulation replicate", "s1", "s2", "p_1way_s1.s2", "right", "distance_t.s1")
colnames(df5c) <- c("pop_no", "qpAdm setup", "gene flow intensity per generation", "simulation replicate", "s1", "s2", "p_1way_s1.s2", "right", "distance_t.s1") #s1*s2 pairs

dff <- rbind(df5a, df5b, df5c)
dfff <- dff[with(dff, order(`gene flow intensity per generation`, `simulation replicate`, right, s1, s2)), ]
df1 <- dfff[!duplicated(dfff[, c("gene flow intensity per generation", "simulation replicate", "right", "s1", "s2")]), ]
rm(df5a)
rm(df5b)
rm(df5c)
rm(dff)
rm(dfff)

df1$positive01 <- ifelse(df1$p_1way_s1.s2>=0.001, "fitting", "non-fitting")
df1$positive02 <- ifelse(df1$p_1way_s1.s2>=0.01, "fitting", "non-fitting")
df1$positive03 <- ifelse(df1$p_1way_s1.s2>=0.05, "fitting", "non-fitting")
df1$positive04 <- ifelse(df1$p_1way_s1.s2>=0.1, "fitting", "non-fitting")
df1$positive05 <- ifelse(df1$p_1way_s1.s2>=0.5, "fitting", "non-fitting")
df1$positive06 <- ifelse(df1$p_1way_s1.s2>=0.05, "fitting", "non-fitting")
df1$positive07 <- ifelse(df1$p_1way_s1.s2>=0.5, "fitting", "non-fitting")
df1$positive08 <- df1$positive01
df1$positive09 <- df1$positive02
df1$positive10 <- df1$positive03
df1$positive11 <- df1$positive04
df1$positive12 <- df1$positive05
df1$positive13 <- df1$positive06
df1$positive14 <- df1$positive07
df1$positive15 <- df1$positive01
df1$positive16 <- df1$positive02
df1$positive17 <- df1$positive03
df1$positive18 <- df1$positive04
df1$positive19 <- df1$positive05
df1$positive20 <- df1$positive06
df1$positive21 <- df1$positive07
df1$positive22 <- df1$positive01
df1$positive23 <- df1$positive02
df1$positive24 <- df1$positive03
df1$positive25 <- df1$positive04
df1$positive26 <- df1$positive05
df1$positive27 <- df1$positive01
df1$positive28 <- df1$positive02
df1$positive29 <- df1$positive03
df1$positive30 <- df1$positive04
df1$positive31 <- df1$positive05
df1$positive32 <- df1$positive01
df1$positive33 <- df1$positive02
df1$positive34 <- df1$positive03
df1$positive35 <- df1$positive04
df1$positive36 <- df1$positive05

###renaming and rearranging columns
df2$model <- "2-way"
df2 <- df2 %>% relocate(model, .after=target)
df3$model <- "3-way"
df3 <- df3 %>% relocate(model, .after=target)
df4$model <- "4-way"
df4 <- df4 %>% relocate(model, .after=target)

df1$`max. distance` <- df1$distance_t.s1
df1$`avg. distance` <- df1$distance_t.s1
df1_1 <- df1
df1_2 <- df1
colnames(df1_1)[colnames(df1_1) == "s1"] <- "target"
colnames(df1_1)[colnames(df1_1) == "s2"] <- "s1"
colnames(df1_2)[colnames(df1_2) == "s2"] <- "target"
df1_2 <- df1_2 %>% relocate(target, .before=s1)
df1 <- rbind(df1_1, df1_2)
df1$model <- "1-way"
df1 <- df1 %>% relocate(model, .after=target)
df1 <- df1 %>% relocate(`max. distance`, .after=distance_t.s1)
df1 <- df1 %>% relocate(`avg. distance`, .after=`max. distance`)
rm(df1_1)
rm(df1_2)

###merging all the tables and selecting targets
df1 <- df1[,c(1:7,9,11:48)]
df2 <- df2[,c(1:8,12,16:55)]
dfa <- rbind.fill(df1, df2, df3, df4)
rm(df1)
rm(df2)
rm(df3)
rm(df4)

###subdividing the results into experiments:
dfa <- dfa[with(dfa, order(pop_no, `gene flow intensity per generation`, `simulation replicate`, right, target, model, s1, s2, s3, s4)), ]
columns <-c("simulation replicate", "max. distance", "avg. distance", "dist_to_corner_maxdist", "dist_to_corner_avgdist")
dfa[, columns] <- lapply(columns, function(x) as.numeric(as.character(dfa[[x]])))
dfa234 <- filter(dfa, model!="1-way")
targets <- (dfa234[!duplicated(dfa234[, c("pop_no", "qpAdm setup", "gene flow intensity per generation", "simulation replicate", "right", "target")]), ])[,c(1:5,8)] %>% set_rownames(NULL)

###counting fitting and non-fitting, optimal and non-optimal models per experiment:
combs <- c(1:length(targets$target))
outtab = sapply(combs, function(x){
test <- filter(dfa, pop_no==targets[x,1], `qpAdm setup`==targets[x,2], `gene flow intensity per generation`==targets[x,3], `simulation replicate`==targets[x,4], target==targets[x,5], right==targets[x,6])
test2f01 <- filter(test, model=="2-way", positive01=="fitting")
test2f02 <- filter(test, model=="2-way", positive02=="fitting")
test2f03 <- filter(test, model=="2-way", positive03=="fitting")
test2f04 <- filter(test, model=="2-way", positive04=="fitting")
test2f05 <- filter(test, model=="2-way", positive05=="fitting")
test2f06 <- filter(test, model=="2-way", positive06=="fitting")
test2f07 <- filter(test, model=="2-way", positive07=="fitting")
test2f08 <- filter(test, model=="2-way", positive08=="fitting")
test2f09 <- filter(test, model=="2-way", positive09=="fitting")
test2f10 <- filter(test, model=="2-way", positive10=="fitting")
test2f11 <- filter(test, model=="2-way", positive11=="fitting")
test2f12 <- filter(test, model=="2-way", positive12=="fitting")
test2f13 <- filter(test, model=="2-way", positive13=="fitting")
test2f14 <- filter(test, model=="2-way", positive14=="fitting")
test2f15 <- filter(test, model=="2-way", positive15=="fitting")
test2f16 <- filter(test, model=="2-way", positive16=="fitting")
test2f17 <- filter(test, model=="2-way", positive17=="fitting")
test2f18 <- filter(test, model=="2-way", positive18=="fitting")
test2f19 <- filter(test, model=="2-way", positive19=="fitting")
test2f20 <- filter(test, model=="2-way", positive20=="fitting")
test2f21 <- filter(test, model=="2-way", positive21=="fitting")
test2f22 <- filter(test, model=="2-way", positive22=="fitting")
test2f23 <- filter(test, model=="2-way", positive23=="fitting")
test2f24 <- filter(test, model=="2-way", positive24=="fitting")
test2f25 <- filter(test, model=="2-way", positive25=="fitting")
test2f26 <- filter(test, model=="2-way", positive26=="fitting")
test2f27 <- filter(test, model=="2-way", positive27=="fitting")
test2f28 <- filter(test, model=="2-way", positive28=="fitting")
test2f29 <- filter(test, model=="2-way", positive29=="fitting")
test2f30 <- filter(test, model=="2-way", positive30=="fitting")
test2f31 <- filter(test, model=="2-way", positive31=="fitting")
test2f32 <- filter(test, model=="2-way", positive32=="fitting")
test2f33 <- filter(test, model=="2-way", positive33=="fitting")
test2f34 <- filter(test, model=="2-way", positive34=="fitting")
test2f35 <- filter(test, model=="2-way", positive35=="fitting")
test2f36 <- filter(test, model=="2-way", positive36=="fitting")

test3f01 <- filter(test, model=="3-way", positive01=="fitting")
test3f02 <- filter(test, model=="3-way", positive02=="fitting")
test3f03 <- filter(test, model=="3-way", positive03=="fitting")
test3f04 <- filter(test, model=="3-way", positive04=="fitting")
test3f05 <- filter(test, model=="3-way", positive05=="fitting")
test3f06 <- filter(test, model=="3-way", positive06=="fitting")
test3f07 <- filter(test, model=="3-way", positive07=="fitting")
test3f08 <- filter(test, model=="3-way", positive08=="fitting")
test3f09 <- filter(test, model=="3-way", positive09=="fitting")
test3f10 <- filter(test, model=="3-way", positive10=="fitting")
test3f11 <- filter(test, model=="3-way", positive11=="fitting")
test3f12 <- filter(test, model=="3-way", positive12=="fitting")
test3f13 <- filter(test, model=="3-way", positive13=="fitting")
test3f14 <- filter(test, model=="3-way", positive14=="fitting")
test3f15 <- filter(test, model=="3-way", positive15=="fitting")
test3f16 <- filter(test, model=="3-way", positive16=="fitting")
test3f17 <- filter(test, model=="3-way", positive17=="fitting")
test3f18 <- filter(test, model=="3-way", positive18=="fitting")
test3f19 <- filter(test, model=="3-way", positive19=="fitting")
test3f20 <- filter(test, model=="3-way", positive20=="fitting")
test3f21 <- filter(test, model=="3-way", positive21=="fitting")
test3f22 <- filter(test, model=="3-way", positive22=="fitting")
test3f23 <- filter(test, model=="3-way", positive23=="fitting")
test3f24 <- filter(test, model=="3-way", positive24=="fitting")
test3f25 <- filter(test, model=="3-way", positive25=="fitting")
test3f26 <- filter(test, model=="3-way", positive26=="fitting")
test3f27 <- filter(test, model=="3-way", positive27=="fitting")
test3f28 <- filter(test, model=="3-way", positive28=="fitting")
test3f29 <- filter(test, model=="3-way", positive29=="fitting")
test3f30 <- filter(test, model=="3-way", positive30=="fitting")
test3f31 <- filter(test, model=="3-way", positive31=="fitting")
test3f32 <- filter(test, model=="3-way", positive32=="fitting")
test3f33 <- filter(test, model=="3-way", positive33=="fitting")
test3f34 <- filter(test, model=="3-way", positive34=="fitting")
test3f35 <- filter(test, model=="3-way", positive35=="fitting")
test3f36 <- filter(test, model=="3-way", positive36=="fitting")

test4f01 <- filter(test, model=="4-way", positive01=="fitting")
test4f02 <- filter(test, model=="4-way", positive02=="fitting")
test4f03 <- filter(test, model=="4-way", positive03=="fitting")
test4f04 <- filter(test, model=="4-way", positive04=="fitting")
test4f05 <- filter(test, model=="4-way", positive05=="fitting")
test4f06 <- filter(test, model=="4-way", positive06=="fitting")
test4f07 <- filter(test, model=="4-way", positive07=="fitting")
test4f08 <- filter(test, model=="4-way", positive08=="fitting")
test4f09 <- filter(test, model=="4-way", positive09=="fitting")
test4f10 <- filter(test, model=="4-way", positive10=="fitting")
test4f11 <- filter(test, model=="4-way", positive11=="fitting")
test4f12 <- filter(test, model=="4-way", positive12=="fitting")
test4f13 <- filter(test, model=="4-way", positive13=="fitting")
test4f14 <- filter(test, model=="4-way", positive14=="fitting")
test4f15 <- filter(test, model=="4-way", positive15=="fitting")
test4f16 <- filter(test, model=="4-way", positive16=="fitting")
test4f17 <- filter(test, model=="4-way", positive17=="fitting")
test4f18 <- filter(test, model=="4-way", positive18=="fitting")
test4f19 <- filter(test, model=="4-way", positive19=="fitting")
test4f20 <- filter(test, model=="4-way", positive20=="fitting")
test4f21 <- filter(test, model=="4-way", positive21=="fitting")
test4f22 <- filter(test, model=="4-way", positive22=="fitting")
test4f23 <- filter(test, model=="4-way", positive23=="fitting")
test4f24 <- filter(test, model=="4-way", positive24=="fitting")
test4f25 <- filter(test, model=="4-way", positive25=="fitting")
test4f26 <- filter(test, model=="4-way", positive26=="fitting")
test4f27 <- filter(test, model=="4-way", positive27=="fitting")
test4f28 <- filter(test, model=="4-way", positive28=="fitting")
test4f29 <- filter(test, model=="4-way", positive29=="fitting")
test4f30 <- filter(test, model=="4-way", positive30=="fitting")
test4f31 <- filter(test, model=="4-way", positive31=="fitting")
test4f32 <- filter(test, model=="4-way", positive32=="fitting")
test4f33 <- filter(test, model=="4-way", positive33=="fitting")
test4f34 <- filter(test, model=="4-way", positive34=="fitting")
test4f35 <- filter(test, model=="4-way", positive35=="fitting")
test4f36 <- filter(test, model=="4-way", positive36=="fitting")

test2nf01 <- filter(test, model=="2-way", positive01=="non-fitting")
test2nf02 <- filter(test, model=="2-way", positive02=="non-fitting")
test2nf03 <- filter(test, model=="2-way", positive03=="non-fitting")
test2nf04 <- filter(test, model=="2-way", positive04=="non-fitting")
test2nf05 <- filter(test, model=="2-way", positive05=="non-fitting")
test2nf06 <- filter(test, model=="2-way", positive06=="non-fitting")
test2nf07 <- filter(test, model=="2-way", positive07=="non-fitting")
test2nf08 <- filter(test, model=="2-way", positive08=="non-fitting")
test2nf09 <- filter(test, model=="2-way", positive09=="non-fitting")
test2nf10 <- filter(test, model=="2-way", positive10=="non-fitting")
test2nf11 <- filter(test, model=="2-way", positive11=="non-fitting")
test2nf12 <- filter(test, model=="2-way", positive12=="non-fitting")
test2nf13 <- filter(test, model=="2-way", positive13=="non-fitting")
test2nf14 <- filter(test, model=="2-way", positive14=="non-fitting")
test2nf15 <- filter(test, model=="2-way", positive15=="non-fitting")
test2nf16 <- filter(test, model=="2-way", positive16=="non-fitting")
test2nf17 <- filter(test, model=="2-way", positive17=="non-fitting")
test2nf18 <- filter(test, model=="2-way", positive18=="non-fitting")
test2nf19 <- filter(test, model=="2-way", positive19=="non-fitting")
test2nf20 <- filter(test, model=="2-way", positive20=="non-fitting")
test2nf21 <- filter(test, model=="2-way", positive21=="non-fitting")
test2nf22 <- filter(test, model=="2-way", positive22=="non-fitting")
test2nf23 <- filter(test, model=="2-way", positive23=="non-fitting")
test2nf24 <- filter(test, model=="2-way", positive24=="non-fitting")
test2nf25 <- filter(test, model=="2-way", positive25=="non-fitting")
test2nf26 <- filter(test, model=="2-way", positive26=="non-fitting")
test2nf27 <- filter(test, model=="2-way", positive27=="non-fitting")
test2nf28 <- filter(test, model=="2-way", positive28=="non-fitting")
test2nf29 <- filter(test, model=="2-way", positive29=="non-fitting")
test2nf30 <- filter(test, model=="2-way", positive30=="non-fitting")
test2nf31 <- filter(test, model=="2-way", positive31=="non-fitting")
test2nf32 <- filter(test, model=="2-way", positive32=="non-fitting")
test2nf33 <- filter(test, model=="2-way", positive33=="non-fitting")
test2nf34 <- filter(test, model=="2-way", positive34=="non-fitting")
test2nf35 <- filter(test, model=="2-way", positive35=="non-fitting")
test2nf36 <- filter(test, model=="2-way", positive36=="non-fitting")

test3nf01 <- filter(test, model=="3-way", positive01=="non-fitting")
test3nf02 <- filter(test, model=="3-way", positive02=="non-fitting")
test3nf03 <- filter(test, model=="3-way", positive03=="non-fitting")
test3nf04 <- filter(test, model=="3-way", positive04=="non-fitting")
test3nf05 <- filter(test, model=="3-way", positive05=="non-fitting")
test3nf06 <- filter(test, model=="3-way", positive06=="non-fitting")
test3nf07 <- filter(test, model=="3-way", positive07=="non-fitting")
test3nf08 <- filter(test, model=="3-way", positive08=="non-fitting")
test3nf09 <- filter(test, model=="3-way", positive09=="non-fitting")
test3nf10 <- filter(test, model=="3-way", positive10=="non-fitting")
test3nf11 <- filter(test, model=="3-way", positive11=="non-fitting")
test3nf12 <- filter(test, model=="3-way", positive12=="non-fitting")
test3nf13 <- filter(test, model=="3-way", positive13=="non-fitting")
test3nf14 <- filter(test, model=="3-way", positive14=="non-fitting")
test3nf15 <- filter(test, model=="3-way", positive15=="non-fitting")
test3nf16 <- filter(test, model=="3-way", positive16=="non-fitting")
test3nf17 <- filter(test, model=="3-way", positive17=="non-fitting")
test3nf18 <- filter(test, model=="3-way", positive18=="non-fitting")
test3nf19 <- filter(test, model=="3-way", positive19=="non-fitting")
test3nf20 <- filter(test, model=="3-way", positive20=="non-fitting")
test3nf21 <- filter(test, model=="3-way", positive21=="non-fitting")
test3nf22 <- filter(test, model=="3-way", positive22=="non-fitting")
test3nf23 <- filter(test, model=="3-way", positive23=="non-fitting")
test3nf24 <- filter(test, model=="3-way", positive24=="non-fitting")
test3nf25 <- filter(test, model=="3-way", positive25=="non-fitting")
test3nf26 <- filter(test, model=="3-way", positive26=="non-fitting")
test3nf27 <- filter(test, model=="3-way", positive27=="non-fitting")
test3nf28 <- filter(test, model=="3-way", positive28=="non-fitting")
test3nf29 <- filter(test, model=="3-way", positive29=="non-fitting")
test3nf30 <- filter(test, model=="3-way", positive30=="non-fitting")
test3nf31 <- filter(test, model=="3-way", positive31=="non-fitting")
test3nf32 <- filter(test, model=="3-way", positive32=="non-fitting")
test3nf33 <- filter(test, model=="3-way", positive33=="non-fitting")
test3nf34 <- filter(test, model=="3-way", positive34=="non-fitting")
test3nf35 <- filter(test, model=="3-way", positive35=="non-fitting")
test3nf36 <- filter(test, model=="3-way", positive36=="non-fitting")

test4nf01 <- filter(test, model=="4-way", positive01=="non-fitting")
test4nf02 <- filter(test, model=="4-way", positive02=="non-fitting")
test4nf03 <- filter(test, model=="4-way", positive03=="non-fitting")
test4nf04 <- filter(test, model=="4-way", positive04=="non-fitting")
test4nf05 <- filter(test, model=="4-way", positive05=="non-fitting")
test4nf06 <- filter(test, model=="4-way", positive06=="non-fitting")
test4nf07 <- filter(test, model=="4-way", positive07=="non-fitting")
test4nf08 <- filter(test, model=="4-way", positive08=="non-fitting")
test4nf09 <- filter(test, model=="4-way", positive09=="non-fitting")
test4nf10 <- filter(test, model=="4-way", positive10=="non-fitting")
test4nf11 <- filter(test, model=="4-way", positive11=="non-fitting")
test4nf12 <- filter(test, model=="4-way", positive12=="non-fitting")
test4nf13 <- filter(test, model=="4-way", positive13=="non-fitting")
test4nf14 <- filter(test, model=="4-way", positive14=="non-fitting")
test4nf15 <- filter(test, model=="4-way", positive15=="non-fitting")
test4nf16 <- filter(test, model=="4-way", positive16=="non-fitting")
test4nf17 <- filter(test, model=="4-way", positive17=="non-fitting")
test4nf18 <- filter(test, model=="4-way", positive18=="non-fitting")
test4nf19 <- filter(test, model=="4-way", positive19=="non-fitting")
test4nf20 <- filter(test, model=="4-way", positive20=="non-fitting")
test4nf21 <- filter(test, model=="4-way", positive21=="non-fitting")
test4nf22 <- filter(test, model=="4-way", positive22=="non-fitting")
test4nf23 <- filter(test, model=="4-way", positive23=="non-fitting")
test4nf24 <- filter(test, model=="4-way", positive24=="non-fitting")
test4nf25 <- filter(test, model=="4-way", positive25=="non-fitting")
test4nf26 <- filter(test, model=="4-way", positive26=="non-fitting")
test4nf27 <- filter(test, model=="4-way", positive27=="non-fitting")
test4nf28 <- filter(test, model=="4-way", positive28=="non-fitting")
test4nf29 <- filter(test, model=="4-way", positive29=="non-fitting")
test4nf30 <- filter(test, model=="4-way", positive30=="non-fitting")
test4nf31 <- filter(test, model=="4-way", positive31=="non-fitting")
test4nf32 <- filter(test, model=="4-way", positive32=="non-fitting")
test4nf33 <- filter(test, model=="4-way", positive33=="non-fitting")
test4nf34 <- filter(test, model=="4-way", positive34=="non-fitting")
test4nf35 <- filter(test, model=="4-way", positive35=="non-fitting")
test4nf36 <- filter(test, model=="4-way", positive36=="non-fitting")

test2nf01$diff <- min(test2f01$dist_to_corner_maxdist, na.rm=T)-test2nf01$dist_to_corner_maxdist
test2nf02$diff <- min(test2f02$dist_to_corner_maxdist, na.rm=T)-test2nf02$dist_to_corner_maxdist
test2nf03$diff <- min(test2f03$dist_to_corner_maxdist, na.rm=T)-test2nf03$dist_to_corner_maxdist
test2nf04$diff <- min(test2f04$dist_to_corner_maxdist, na.rm=T)-test2nf04$dist_to_corner_maxdist
test2nf05$diff <- min(test2f05$dist_to_corner_maxdist, na.rm=T)-test2nf05$dist_to_corner_maxdist
test2nf06$diff <- min(test2f06$dist_to_corner_maxdist, na.rm=T)-test2nf06$dist_to_corner_maxdist
test2nf07$diff <- min(test2f07$dist_to_corner_maxdist, na.rm=T)-test2nf07$dist_to_corner_maxdist
test2nf08$diff <- min(test2f08$dist_to_corner_maxdist, na.rm=T)-test2nf08$dist_to_corner_maxdist
test2nf09$diff <- min(test2f09$dist_to_corner_maxdist, na.rm=T)-test2nf09$dist_to_corner_maxdist
test2nf10$diff <- min(test2f10$dist_to_corner_maxdist, na.rm=T)-test2nf10$dist_to_corner_maxdist
test2nf11$diff <- min(test2f11$dist_to_corner_maxdist, na.rm=T)-test2nf11$dist_to_corner_maxdist
test2nf12$diff <- min(test2f12$dist_to_corner_maxdist, na.rm=T)-test2nf12$dist_to_corner_maxdist
test2nf13$diff <- min(test2f13$dist_to_corner_maxdist, na.rm=T)-test2nf13$dist_to_corner_maxdist
test2nf14$diff <- min(test2f14$dist_to_corner_maxdist, na.rm=T)-test2nf14$dist_to_corner_maxdist
test2nf15$diff <- min(test2f15$dist_to_corner_maxdist, na.rm=T)-test2nf15$dist_to_corner_maxdist
test2nf16$diff <- min(test2f16$dist_to_corner_maxdist, na.rm=T)-test2nf16$dist_to_corner_maxdist
test2nf17$diff <- min(test2f17$dist_to_corner_maxdist, na.rm=T)-test2nf17$dist_to_corner_maxdist
test2nf18$diff <- min(test2f18$dist_to_corner_maxdist, na.rm=T)-test2nf18$dist_to_corner_maxdist
test2nf19$diff <- min(test2f19$dist_to_corner_maxdist, na.rm=T)-test2nf19$dist_to_corner_maxdist
test2nf20$diff <- min(test2f20$dist_to_corner_maxdist, na.rm=T)-test2nf20$dist_to_corner_maxdist
test2nf21$diff <- min(test2f21$dist_to_corner_maxdist, na.rm=T)-test2nf21$dist_to_corner_maxdist
test2nf22$diff <- min(test2f22$dist_to_corner_maxdist, na.rm=T)-test2nf22$dist_to_corner_maxdist
test2nf23$diff <- min(test2f23$dist_to_corner_maxdist, na.rm=T)-test2nf23$dist_to_corner_maxdist
test2nf24$diff <- min(test2f24$dist_to_corner_maxdist, na.rm=T)-test2nf24$dist_to_corner_maxdist
test2nf25$diff <- min(test2f25$dist_to_corner_maxdist, na.rm=T)-test2nf25$dist_to_corner_maxdist
test2nf26$diff <- min(test2f26$dist_to_corner_maxdist, na.rm=T)-test2nf26$dist_to_corner_maxdist
test2nf27$diff <- min(test2f27$dist_to_corner_maxdist, na.rm=T)-test2nf27$dist_to_corner_maxdist
test2nf28$diff <- min(test2f28$dist_to_corner_maxdist, na.rm=T)-test2nf28$dist_to_corner_maxdist
test2nf29$diff <- min(test2f29$dist_to_corner_maxdist, na.rm=T)-test2nf29$dist_to_corner_maxdist
test2nf30$diff <- min(test2f30$dist_to_corner_maxdist, na.rm=T)-test2nf30$dist_to_corner_maxdist
test2nf31$diff <- min(test2f31$dist_to_corner_maxdist, na.rm=T)-test2nf31$dist_to_corner_maxdist
test2nf32$diff <- min(test2f32$dist_to_corner_maxdist, na.rm=T)-test2nf32$dist_to_corner_maxdist
test2nf33$diff <- min(test2f33$dist_to_corner_maxdist, na.rm=T)-test2nf33$dist_to_corner_maxdist
test2nf34$diff <- min(test2f34$dist_to_corner_maxdist, na.rm=T)-test2nf34$dist_to_corner_maxdist
test2nf35$diff <- min(test2f35$dist_to_corner_maxdist, na.rm=T)-test2nf35$dist_to_corner_maxdist
test2nf36$diff <- min(test2f36$dist_to_corner_maxdist, na.rm=T)-test2nf36$dist_to_corner_maxdist

test3nf01$diff <- min(test3f01$dist_to_corner_maxdist, na.rm=T)-test3nf01$dist_to_corner_maxdist
test3nf02$diff <- min(test3f02$dist_to_corner_maxdist, na.rm=T)-test3nf02$dist_to_corner_maxdist
test3nf03$diff <- min(test3f03$dist_to_corner_maxdist, na.rm=T)-test3nf03$dist_to_corner_maxdist
test3nf04$diff <- min(test3f04$dist_to_corner_maxdist, na.rm=T)-test3nf04$dist_to_corner_maxdist
test3nf05$diff <- min(test3f05$dist_to_corner_maxdist, na.rm=T)-test3nf05$dist_to_corner_maxdist
test3nf06$diff <- min(test3f06$dist_to_corner_maxdist, na.rm=T)-test3nf06$dist_to_corner_maxdist
test3nf07$diff <- min(test3f07$dist_to_corner_maxdist, na.rm=T)-test3nf07$dist_to_corner_maxdist
test3nf08$diff <- min(test3f08$dist_to_corner_maxdist, na.rm=T)-test3nf08$dist_to_corner_maxdist
test3nf09$diff <- min(test3f09$dist_to_corner_maxdist, na.rm=T)-test3nf09$dist_to_corner_maxdist
test3nf10$diff <- min(test3f10$dist_to_corner_maxdist, na.rm=T)-test3nf10$dist_to_corner_maxdist
test3nf11$diff <- min(test3f11$dist_to_corner_maxdist, na.rm=T)-test3nf11$dist_to_corner_maxdist
test3nf12$diff <- min(test3f12$dist_to_corner_maxdist, na.rm=T)-test3nf12$dist_to_corner_maxdist
test3nf13$diff <- min(test3f13$dist_to_corner_maxdist, na.rm=T)-test3nf13$dist_to_corner_maxdist
test3nf14$diff <- min(test3f14$dist_to_corner_maxdist, na.rm=T)-test3nf14$dist_to_corner_maxdist
test3nf15$diff <- min(test3f15$dist_to_corner_maxdist, na.rm=T)-test3nf15$dist_to_corner_maxdist
test3nf16$diff <- min(test3f16$dist_to_corner_maxdist, na.rm=T)-test3nf16$dist_to_corner_maxdist
test3nf17$diff <- min(test3f17$dist_to_corner_maxdist, na.rm=T)-test3nf17$dist_to_corner_maxdist
test3nf18$diff <- min(test3f18$dist_to_corner_maxdist, na.rm=T)-test3nf18$dist_to_corner_maxdist
test3nf19$diff <- min(test3f19$dist_to_corner_maxdist, na.rm=T)-test3nf19$dist_to_corner_maxdist
test3nf20$diff <- min(test3f20$dist_to_corner_maxdist, na.rm=T)-test3nf20$dist_to_corner_maxdist
test3nf21$diff <- min(test3f21$dist_to_corner_maxdist, na.rm=T)-test3nf21$dist_to_corner_maxdist
test3nf22$diff <- min(test3f22$dist_to_corner_maxdist, na.rm=T)-test3nf22$dist_to_corner_maxdist
test3nf23$diff <- min(test3f23$dist_to_corner_maxdist, na.rm=T)-test3nf23$dist_to_corner_maxdist
test3nf24$diff <- min(test3f24$dist_to_corner_maxdist, na.rm=T)-test3nf24$dist_to_corner_maxdist
test3nf25$diff <- min(test3f25$dist_to_corner_maxdist, na.rm=T)-test3nf25$dist_to_corner_maxdist
test3nf26$diff <- min(test3f26$dist_to_corner_maxdist, na.rm=T)-test3nf26$dist_to_corner_maxdist
test3nf27$diff <- min(test3f27$dist_to_corner_maxdist, na.rm=T)-test3nf27$dist_to_corner_maxdist
test3nf28$diff <- min(test3f28$dist_to_corner_maxdist, na.rm=T)-test3nf28$dist_to_corner_maxdist
test3nf29$diff <- min(test3f29$dist_to_corner_maxdist, na.rm=T)-test3nf29$dist_to_corner_maxdist
test3nf30$diff <- min(test3f30$dist_to_corner_maxdist, na.rm=T)-test3nf30$dist_to_corner_maxdist
test3nf31$diff <- min(test3f31$dist_to_corner_maxdist, na.rm=T)-test3nf31$dist_to_corner_maxdist
test3nf32$diff <- min(test3f32$dist_to_corner_maxdist, na.rm=T)-test3nf32$dist_to_corner_maxdist
test3nf33$diff <- min(test3f33$dist_to_corner_maxdist, na.rm=T)-test3nf33$dist_to_corner_maxdist
test3nf34$diff <- min(test3f34$dist_to_corner_maxdist, na.rm=T)-test3nf34$dist_to_corner_maxdist
test3nf35$diff <- min(test3f35$dist_to_corner_maxdist, na.rm=T)-test3nf35$dist_to_corner_maxdist
test3nf36$diff <- min(test3f36$dist_to_corner_maxdist, na.rm=T)-test3nf36$dist_to_corner_maxdist

test4nf01$diff <- min(test4f01$dist_to_corner_maxdist, na.rm=T)-test4nf01$dist_to_corner_maxdist
test4nf02$diff <- min(test4f02$dist_to_corner_maxdist, na.rm=T)-test4nf02$dist_to_corner_maxdist
test4nf03$diff <- min(test4f03$dist_to_corner_maxdist, na.rm=T)-test4nf03$dist_to_corner_maxdist
test4nf04$diff <- min(test4f04$dist_to_corner_maxdist, na.rm=T)-test4nf04$dist_to_corner_maxdist
test4nf05$diff <- min(test4f05$dist_to_corner_maxdist, na.rm=T)-test4nf05$dist_to_corner_maxdist
test4nf06$diff <- min(test4f06$dist_to_corner_maxdist, na.rm=T)-test4nf06$dist_to_corner_maxdist
test4nf07$diff <- min(test4f07$dist_to_corner_maxdist, na.rm=T)-test4nf07$dist_to_corner_maxdist
test4nf08$diff <- min(test4f08$dist_to_corner_maxdist, na.rm=T)-test4nf08$dist_to_corner_maxdist
test4nf09$diff <- min(test4f09$dist_to_corner_maxdist, na.rm=T)-test4nf09$dist_to_corner_maxdist
test4nf10$diff <- min(test4f10$dist_to_corner_maxdist, na.rm=T)-test4nf10$dist_to_corner_maxdist
test4nf11$diff <- min(test4f11$dist_to_corner_maxdist, na.rm=T)-test4nf11$dist_to_corner_maxdist
test4nf12$diff <- min(test4f12$dist_to_corner_maxdist, na.rm=T)-test4nf12$dist_to_corner_maxdist
test4nf13$diff <- min(test4f13$dist_to_corner_maxdist, na.rm=T)-test4nf13$dist_to_corner_maxdist
test4nf14$diff <- min(test4f14$dist_to_corner_maxdist, na.rm=T)-test4nf14$dist_to_corner_maxdist
test4nf15$diff <- min(test4f15$dist_to_corner_maxdist, na.rm=T)-test4nf15$dist_to_corner_maxdist
test4nf16$diff <- min(test4f16$dist_to_corner_maxdist, na.rm=T)-test4nf16$dist_to_corner_maxdist
test4nf17$diff <- min(test4f17$dist_to_corner_maxdist, na.rm=T)-test4nf17$dist_to_corner_maxdist
test4nf18$diff <- min(test4f18$dist_to_corner_maxdist, na.rm=T)-test4nf18$dist_to_corner_maxdist
test4nf19$diff <- min(test4f19$dist_to_corner_maxdist, na.rm=T)-test4nf19$dist_to_corner_maxdist
test4nf20$diff <- min(test4f20$dist_to_corner_maxdist, na.rm=T)-test4nf20$dist_to_corner_maxdist
test4nf21$diff <- min(test4f21$dist_to_corner_maxdist, na.rm=T)-test4nf21$dist_to_corner_maxdist
test4nf22$diff <- min(test4f22$dist_to_corner_maxdist, na.rm=T)-test4nf22$dist_to_corner_maxdist
test4nf23$diff <- min(test4f23$dist_to_corner_maxdist, na.rm=T)-test4nf23$dist_to_corner_maxdist
test4nf24$diff <- min(test4f24$dist_to_corner_maxdist, na.rm=T)-test4nf24$dist_to_corner_maxdist
test4nf25$diff <- min(test4f25$dist_to_corner_maxdist, na.rm=T)-test4nf25$dist_to_corner_maxdist
test4nf26$diff <- min(test4f26$dist_to_corner_maxdist, na.rm=T)-test4nf26$dist_to_corner_maxdist
test4nf27$diff <- min(test4f27$dist_to_corner_maxdist, na.rm=T)-test4nf27$dist_to_corner_maxdist
test4nf28$diff <- min(test4f28$dist_to_corner_maxdist, na.rm=T)-test4nf28$dist_to_corner_maxdist
test4nf29$diff <- min(test4f29$dist_to_corner_maxdist, na.rm=T)-test4nf29$dist_to_corner_maxdist
test4nf30$diff <- min(test4f30$dist_to_corner_maxdist, na.rm=T)-test4nf30$dist_to_corner_maxdist
test4nf31$diff <- min(test4f31$dist_to_corner_maxdist, na.rm=T)-test4nf31$dist_to_corner_maxdist
test4nf32$diff <- min(test4f32$dist_to_corner_maxdist, na.rm=T)-test4nf32$dist_to_corner_maxdist
test4nf33$diff <- min(test4f33$dist_to_corner_maxdist, na.rm=T)-test4nf33$dist_to_corner_maxdist
test4nf34$diff <- min(test4f34$dist_to_corner_maxdist, na.rm=T)-test4nf34$dist_to_corner_maxdist
test4nf35$diff <- min(test4f35$dist_to_corner_maxdist, na.rm=T)-test4nf35$dist_to_corner_maxdist
test4nf36$diff <- min(test4f36$dist_to_corner_maxdist, na.rm=T)-test4nf36$dist_to_corner_maxdist

c("pop_no" = targets[x,1],
  "qpAdm setup" = targets[x,2],
  "gene flow intensity per generation" = targets[x,3],
  "simulation replicate" = targets[x,4],
  "target" = targets[x,5],
  "right" = targets[x,6],
  "N1" = nrow(test[test$model == '1-way', ]),
  "PC1_01" = nrow(test[test$model == '1-way' & test$positive01 == 'fitting', ]),
  "PC1_02" = nrow(test[test$model == '1-way' & test$positive02 == 'fitting', ]),
  "PC1_03" = nrow(test[test$model == '1-way' & test$positive03 == 'fitting', ]),
  "PC1_04" = nrow(test[test$model == '1-way' & test$positive04 == 'fitting', ]),
  "PC1_05" = nrow(test[test$model == '1-way' & test$positive05 == 'fitting', ]),
  "PC1_06" = nrow(test[test$model == '1-way' & test$positive06 == 'fitting', ]),
  "PC1_07" = nrow(test[test$model == '1-way' & test$positive07 == 'fitting', ]),
  "PC1_08" = nrow(test[test$model == '1-way' & test$positive08 == 'fitting', ]),
  "PC1_09" = nrow(test[test$model == '1-way' & test$positive09 == 'fitting', ]),
  "PC1_10" = nrow(test[test$model == '1-way' & test$positive10 == 'fitting', ]),
  "PC1_11" = nrow(test[test$model == '1-way' & test$positive11 == 'fitting', ]),
  "PC1_12" = nrow(test[test$model == '1-way' & test$positive12 == 'fitting', ]),
  "PC1_13" = nrow(test[test$model == '1-way' & test$positive13 == 'fitting', ]),
  "PC1_14" = nrow(test[test$model == '1-way' & test$positive14 == 'fitting', ]),
  "PC1_15" = nrow(test[test$model == '1-way' & test$positive15 == 'fitting', ]),
  "PC1_16" = nrow(test[test$model == '1-way' & test$positive16 == 'fitting', ]),
  "PC1_17" = nrow(test[test$model == '1-way' & test$positive17 == 'fitting', ]),
  "PC1_18" = nrow(test[test$model == '1-way' & test$positive18 == 'fitting', ]),
  "PC1_19" = nrow(test[test$model == '1-way' & test$positive19 == 'fitting', ]),
  "PC1_20" = nrow(test[test$model == '1-way' & test$positive20 == 'fitting', ]),
  "PC1_21" = nrow(test[test$model == '1-way' & test$positive21 == 'fitting', ]),
  "PC1_22" = nrow(test[test$model == '1-way' & test$positive22 == 'fitting', ]),
  "PC1_23" = nrow(test[test$model == '1-way' & test$positive23 == 'fitting', ]),
  "PC1_24" = nrow(test[test$model == '1-way' & test$positive24 == 'fitting', ]),
  "PC1_25" = nrow(test[test$model == '1-way' & test$positive25 == 'fitting', ]),
  "PC1_26" = nrow(test[test$model == '1-way' & test$positive26 == 'fitting', ]),
  "PC1_27" = nrow(test[test$model == '1-way' & test$positive27 == 'fitting', ]),
  "PC1_28" = nrow(test[test$model == '1-way' & test$positive28 == 'fitting', ]),
  "PC1_29" = nrow(test[test$model == '1-way' & test$positive29 == 'fitting', ]),
  "PC1_30" = nrow(test[test$model == '1-way' & test$positive30 == 'fitting', ]),
  "PC1_31" = nrow(test[test$model == '1-way' & test$positive31 == 'fitting', ]),
  "PC1_32" = nrow(test[test$model == '1-way' & test$positive32 == 'fitting', ]),
  "PC1_33" = nrow(test[test$model == '1-way' & test$positive33 == 'fitting', ]),
  "PC1_34" = nrow(test[test$model == '1-way' & test$positive34 == 'fitting', ]),
  "PC1_35" = nrow(test[test$model == '1-way' & test$positive35 == 'fitting', ]),
  "PC1_36" = nrow(test[test$model == '1-way' & test$positive36 == 'fitting', ]),

  "N2" = nrow(test[test$model == '2-way', ]),
  "PC2_01" = nrow(test[test$model == '2-way' & test$positive01 == 'fitting', ]),
  "PC2_02" = nrow(test[test$model == '2-way' & test$positive02 == 'fitting', ]),
  "PC2_03" = nrow(test[test$model == '2-way' & test$positive03 == 'fitting', ]),
  "PC2_04" = nrow(test[test$model == '2-way' & test$positive04 == 'fitting', ]),
  "PC2_05" = nrow(test[test$model == '2-way' & test$positive05 == 'fitting', ]),
  "PC2_06" = nrow(test[test$model == '2-way' & test$positive06 == 'fitting', ]),
  "PC2_07" = nrow(test[test$model == '2-way' & test$positive07 == 'fitting', ]),
  "PC2_08" = nrow(test[test$model == '2-way' & test$positive08 == 'fitting', ]),
  "PC2_09" = nrow(test[test$model == '2-way' & test$positive09 == 'fitting', ]),
  "PC2_10" = nrow(test[test$model == '2-way' & test$positive10 == 'fitting', ]),
  "PC2_11" = nrow(test[test$model == '2-way' & test$positive11 == 'fitting', ]),
  "PC2_12" = nrow(test[test$model == '2-way' & test$positive12 == 'fitting', ]),
  "PC2_13" = nrow(test[test$model == '2-way' & test$positive13 == 'fitting', ]),
  "PC2_14" = nrow(test[test$model == '2-way' & test$positive14 == 'fitting', ]),
  "PC2_15" = nrow(test[test$model == '2-way' & test$positive15 == 'fitting', ]),
  "PC2_16" = nrow(test[test$model == '2-way' & test$positive16 == 'fitting', ]),
  "PC2_17" = nrow(test[test$model == '2-way' & test$positive17 == 'fitting', ]),
  "PC2_18" = nrow(test[test$model == '2-way' & test$positive18 == 'fitting', ]),
  "PC2_19" = nrow(test[test$model == '2-way' & test$positive19 == 'fitting', ]),
  "PC2_20" = nrow(test[test$model == '2-way' & test$positive20 == 'fitting', ]),
  "PC2_21" = nrow(test[test$model == '2-way' & test$positive21 == 'fitting', ]),
  "PC2_22" = nrow(test[test$model == '2-way' & test$positive22 == 'fitting', ]),
  "PC2_23" = nrow(test[test$model == '2-way' & test$positive23 == 'fitting', ]),
  "PC2_24" = nrow(test[test$model == '2-way' & test$positive24 == 'fitting', ]),
  "PC2_25" = nrow(test[test$model == '2-way' & test$positive25 == 'fitting', ]),
  "PC2_26" = nrow(test[test$model == '2-way' & test$positive26 == 'fitting', ]),
  "PC2_27" = nrow(test[test$model == '2-way' & test$positive27 == 'fitting', ]),
  "PC2_28" = nrow(test[test$model == '2-way' & test$positive28 == 'fitting', ]),
  "PC2_29" = nrow(test[test$model == '2-way' & test$positive29 == 'fitting', ]),
  "PC2_30" = nrow(test[test$model == '2-way' & test$positive30 == 'fitting', ]),
  "PC2_31" = nrow(test[test$model == '2-way' & test$positive31 == 'fitting', ]),
  "PC2_32" = nrow(test[test$model == '2-way' & test$positive32 == 'fitting', ]),
  "PC2_33" = nrow(test[test$model == '2-way' & test$positive33 == 'fitting', ]),
  "PC2_34" = nrow(test[test$model == '2-way' & test$positive34 == 'fitting', ]),
  "PC2_35" = nrow(test[test$model == '2-way' & test$positive35 == 'fitting', ]),
  "PC2_36" = nrow(test[test$model == '2-way' & test$positive36 == 'fitting', ]),

  "N3" = nrow(test[test$model == '3-way', ]),
  "PC3_01" = nrow(test[test$model == '3-way' & test$positive01 == 'fitting', ]),
  "PC3_02" = nrow(test[test$model == '3-way' & test$positive02 == 'fitting', ]),
  "PC3_03" = nrow(test[test$model == '3-way' & test$positive03 == 'fitting', ]),
  "PC3_04" = nrow(test[test$model == '3-way' & test$positive04 == 'fitting', ]),
  "PC3_05" = nrow(test[test$model == '3-way' & test$positive05 == 'fitting', ]),
  "PC3_06" = nrow(test[test$model == '3-way' & test$positive06 == 'fitting', ]),
  "PC3_07" = nrow(test[test$model == '3-way' & test$positive07 == 'fitting', ]),
  "PC3_08" = nrow(test[test$model == '3-way' & test$positive08 == 'fitting', ]),
  "PC3_09" = nrow(test[test$model == '3-way' & test$positive09 == 'fitting', ]),
  "PC3_10" = nrow(test[test$model == '3-way' & test$positive10 == 'fitting', ]),
  "PC3_11" = nrow(test[test$model == '3-way' & test$positive11 == 'fitting', ]),
  "PC3_12" = nrow(test[test$model == '3-way' & test$positive12 == 'fitting', ]),
  "PC3_13" = nrow(test[test$model == '3-way' & test$positive13 == 'fitting', ]),
  "PC3_14" = nrow(test[test$model == '3-way' & test$positive14 == 'fitting', ]),
  "PC3_15" = nrow(test[test$model == '3-way' & test$positive15 == 'fitting', ]),
  "PC3_16" = nrow(test[test$model == '3-way' & test$positive16 == 'fitting', ]),
  "PC3_17" = nrow(test[test$model == '3-way' & test$positive17 == 'fitting', ]),
  "PC3_18" = nrow(test[test$model == '3-way' & test$positive18 == 'fitting', ]),
  "PC3_19" = nrow(test[test$model == '3-way' & test$positive19 == 'fitting', ]),
  "PC3_20" = nrow(test[test$model == '3-way' & test$positive20 == 'fitting', ]),
  "PC3_21" = nrow(test[test$model == '3-way' & test$positive21 == 'fitting', ]),
  "PC3_22" = nrow(test[test$model == '3-way' & test$positive22 == 'fitting', ]),
  "PC3_23" = nrow(test[test$model == '3-way' & test$positive23 == 'fitting', ]),
  "PC3_24" = nrow(test[test$model == '3-way' & test$positive24 == 'fitting', ]),
  "PC3_25" = nrow(test[test$model == '3-way' & test$positive25 == 'fitting', ]),
  "PC3_26" = nrow(test[test$model == '3-way' & test$positive26 == 'fitting', ]),
  "PC3_27" = nrow(test[test$model == '3-way' & test$positive27 == 'fitting', ]),
  "PC3_28" = nrow(test[test$model == '3-way' & test$positive28 == 'fitting', ]),
  "PC3_29" = nrow(test[test$model == '3-way' & test$positive29 == 'fitting', ]),
  "PC3_30" = nrow(test[test$model == '3-way' & test$positive30 == 'fitting', ]),
  "PC3_31" = nrow(test[test$model == '3-way' & test$positive31 == 'fitting', ]),
  "PC3_32" = nrow(test[test$model == '3-way' & test$positive32 == 'fitting', ]),
  "PC3_33" = nrow(test[test$model == '3-way' & test$positive33 == 'fitting', ]),
  "PC3_34" = nrow(test[test$model == '3-way' & test$positive34 == 'fitting', ]),
  "PC3_35" = nrow(test[test$model == '3-way' & test$positive35 == 'fitting', ]),
  "PC3_36" = nrow(test[test$model == '3-way' & test$positive36 == 'fitting', ]),

  "N4" = nrow(test[test$model == '4-way', ]),
  "PC4_01" = nrow(test[test$model == '4-way' & test$positive01 == 'fitting', ]),
  "PC4_02" = nrow(test[test$model == '4-way' & test$positive02 == 'fitting', ]),
  "PC4_03" = nrow(test[test$model == '4-way' & test$positive03 == 'fitting', ]),
  "PC4_04" = nrow(test[test$model == '4-way' & test$positive04 == 'fitting', ]),
  "PC4_05" = nrow(test[test$model == '4-way' & test$positive05 == 'fitting', ]),
  "PC4_06" = nrow(test[test$model == '4-way' & test$positive06 == 'fitting', ]),
  "PC4_07" = nrow(test[test$model == '4-way' & test$positive07 == 'fitting', ]),
  "PC4_08" = nrow(test[test$model == '4-way' & test$positive08 == 'fitting', ]),
  "PC4_09" = nrow(test[test$model == '4-way' & test$positive09 == 'fitting', ]),
  "PC4_10" = nrow(test[test$model == '4-way' & test$positive10 == 'fitting', ]),
  "PC4_11" = nrow(test[test$model == '4-way' & test$positive11 == 'fitting', ]),
  "PC4_12" = nrow(test[test$model == '4-way' & test$positive12 == 'fitting', ]),
  "PC4_13" = nrow(test[test$model == '4-way' & test$positive13 == 'fitting', ]),
  "PC4_14" = nrow(test[test$model == '4-way' & test$positive14 == 'fitting', ]),
  "PC4_15" = nrow(test[test$model == '4-way' & test$positive15 == 'fitting', ]),
  "PC4_16" = nrow(test[test$model == '4-way' & test$positive16 == 'fitting', ]),
  "PC4_17" = nrow(test[test$model == '4-way' & test$positive17 == 'fitting', ]),
  "PC4_18" = nrow(test[test$model == '4-way' & test$positive18 == 'fitting', ]),
  "PC4_19" = nrow(test[test$model == '4-way' & test$positive19 == 'fitting', ]),
  "PC4_20" = nrow(test[test$model == '4-way' & test$positive20 == 'fitting', ]),
  "PC4_21" = nrow(test[test$model == '4-way' & test$positive21 == 'fitting', ]),
  "PC4_22" = nrow(test[test$model == '4-way' & test$positive22 == 'fitting', ]),
  "PC4_23" = nrow(test[test$model == '4-way' & test$positive23 == 'fitting', ]),
  "PC4_24" = nrow(test[test$model == '4-way' & test$positive24 == 'fitting', ]),
  "PC4_25" = nrow(test[test$model == '4-way' & test$positive25 == 'fitting', ]),
  "PC4_26" = nrow(test[test$model == '4-way' & test$positive26 == 'fitting', ]),
  "PC4_27" = nrow(test[test$model == '4-way' & test$positive27 == 'fitting', ]),
  "PC4_28" = nrow(test[test$model == '4-way' & test$positive28 == 'fitting', ]),
  "PC4_29" = nrow(test[test$model == '4-way' & test$positive29 == 'fitting', ]),
  "PC4_30" = nrow(test[test$model == '4-way' & test$positive30 == 'fitting', ]),
  "PC4_31" = nrow(test[test$model == '4-way' & test$positive31 == 'fitting', ]),
  "PC4_32" = nrow(test[test$model == '4-way' & test$positive32 == 'fitting', ]),
  "PC4_33" = nrow(test[test$model == '4-way' & test$positive33 == 'fitting', ]),
  "PC4_34" = nrow(test[test$model == '4-way' & test$positive34 == 'fitting', ]),
  "PC4_35" = nrow(test[test$model == '4-way' & test$positive35 == 'fitting', ]),
  "PC4_36" = nrow(test[test$model == '4-way' & test$positive36 == 'fitting', ]),

  "FNC2_01" = nrow(test2nf01[test2nf01$diff>sqrt(8) & !is.na(test2nf01$diff) & test2nf01$diff!=Inf & test2nf01$diff!=-Inf, ]),
  "FNC2_02" = nrow(test2nf02[test2nf02$diff>sqrt(8) & !is.na(test2nf02$diff) & test2nf02$diff!=Inf & test2nf02$diff!=-Inf, ]),
  "FNC2_03" = nrow(test2nf03[test2nf03$diff>sqrt(8) & !is.na(test2nf03$diff) & test2nf03$diff!=Inf & test2nf03$diff!=-Inf, ]),
  "FNC2_04" = nrow(test2nf04[test2nf04$diff>sqrt(8) & !is.na(test2nf04$diff) & test2nf04$diff!=Inf & test2nf04$diff!=-Inf, ]),
  "FNC2_05" = nrow(test2nf05[test2nf05$diff>sqrt(8) & !is.na(test2nf05$diff) & test2nf05$diff!=Inf & test2nf05$diff!=-Inf, ]),
  "FNC2_06" = nrow(test2nf06[test2nf06$diff>sqrt(8) & !is.na(test2nf06$diff) & test2nf06$diff!=Inf & test2nf06$diff!=-Inf, ]),
  "FNC2_07" = nrow(test2nf07[test2nf07$diff>sqrt(8) & !is.na(test2nf07$diff) & test2nf07$diff!=Inf & test2nf07$diff!=-Inf, ]),
  "FNC2_08" = nrow(test2nf08[test2nf08$diff>sqrt(8) & !is.na(test2nf08$diff) & test2nf08$diff!=Inf & test2nf08$diff!=-Inf, ]),
  "FNC2_09" = nrow(test2nf09[test2nf09$diff>sqrt(8) & !is.na(test2nf09$diff) & test2nf09$diff!=Inf & test2nf09$diff!=-Inf, ]),
  "FNC2_10" = nrow(test2nf10[test2nf10$diff>sqrt(8) & !is.na(test2nf10$diff) & test2nf10$diff!=Inf & test2nf10$diff!=-Inf, ]),
  "FNC2_11" = nrow(test2nf11[test2nf11$diff>sqrt(8) & !is.na(test2nf11$diff) & test2nf11$diff!=Inf & test2nf11$diff!=-Inf, ]),
  "FNC2_12" = nrow(test2nf12[test2nf12$diff>sqrt(8) & !is.na(test2nf12$diff) & test2nf12$diff!=Inf & test2nf12$diff!=-Inf, ]),
  "FNC2_13" = nrow(test2nf13[test2nf13$diff>sqrt(8) & !is.na(test2nf13$diff) & test2nf13$diff!=Inf & test2nf13$diff!=-Inf, ]),
  "FNC2_14" = nrow(test2nf14[test2nf14$diff>sqrt(8) & !is.na(test2nf14$diff) & test2nf14$diff!=Inf & test2nf14$diff!=-Inf, ]),
  "FNC2_15" = nrow(test2nf15[test2nf15$diff>sqrt(8) & !is.na(test2nf15$diff) & test2nf15$diff!=Inf & test2nf15$diff!=-Inf, ]),
  "FNC2_16" = nrow(test2nf16[test2nf16$diff>sqrt(8) & !is.na(test2nf16$diff) & test2nf16$diff!=Inf & test2nf16$diff!=-Inf, ]),
  "FNC2_17" = nrow(test2nf17[test2nf17$diff>sqrt(8) & !is.na(test2nf17$diff) & test2nf17$diff!=Inf & test2nf17$diff!=-Inf, ]),
  "FNC2_18" = nrow(test2nf18[test2nf18$diff>sqrt(8) & !is.na(test2nf18$diff) & test2nf18$diff!=Inf & test2nf18$diff!=-Inf, ]),
  "FNC2_19" = nrow(test2nf19[test2nf19$diff>sqrt(8) & !is.na(test2nf19$diff) & test2nf19$diff!=Inf & test2nf19$diff!=-Inf, ]),
  "FNC2_20" = nrow(test2nf20[test2nf20$diff>sqrt(8) & !is.na(test2nf20$diff) & test2nf20$diff!=Inf & test2nf20$diff!=-Inf, ]),
  "FNC2_21" = nrow(test2nf21[test2nf21$diff>sqrt(8) & !is.na(test2nf21$diff) & test2nf21$diff!=Inf & test2nf21$diff!=-Inf, ]),
  "FNC2_22" = nrow(test2nf22[test2nf22$diff>sqrt(8) & !is.na(test2nf22$diff) & test2nf22$diff!=Inf & test2nf22$diff!=-Inf, ]),
  "FNC2_23" = nrow(test2nf23[test2nf23$diff>sqrt(8) & !is.na(test2nf23$diff) & test2nf23$diff!=Inf & test2nf23$diff!=-Inf, ]),
  "FNC2_24" = nrow(test2nf24[test2nf24$diff>sqrt(8) & !is.na(test2nf24$diff) & test2nf24$diff!=Inf & test2nf24$diff!=-Inf, ]),
  "FNC2_25" = nrow(test2nf25[test2nf25$diff>sqrt(8) & !is.na(test2nf25$diff) & test2nf25$diff!=Inf & test2nf25$diff!=-Inf, ]),
  "FNC2_26" = nrow(test2nf26[test2nf26$diff>sqrt(8) & !is.na(test2nf26$diff) & test2nf26$diff!=Inf & test2nf26$diff!=-Inf, ]),
  "FNC2_27" = nrow(test2nf27[test2nf27$diff>sqrt(8) & !is.na(test2nf27$diff) & test2nf27$diff!=Inf & test2nf27$diff!=-Inf, ]),
  "FNC2_28" = nrow(test2nf28[test2nf28$diff>sqrt(8) & !is.na(test2nf28$diff) & test2nf28$diff!=Inf & test2nf28$diff!=-Inf, ]),
  "FNC2_29" = nrow(test2nf29[test2nf29$diff>sqrt(8) & !is.na(test2nf29$diff) & test2nf29$diff!=Inf & test2nf29$diff!=-Inf, ]),
  "FNC2_30" = nrow(test2nf30[test2nf30$diff>sqrt(8) & !is.na(test2nf30$diff) & test2nf30$diff!=Inf & test2nf30$diff!=-Inf, ]),
  "FNC2_31" = nrow(test2nf31[test2nf31$diff>sqrt(8) & !is.na(test2nf31$diff) & test2nf31$diff!=Inf & test2nf31$diff!=-Inf, ]),
  "FNC2_32" = nrow(test2nf32[test2nf32$diff>sqrt(8) & !is.na(test2nf32$diff) & test2nf32$diff!=Inf & test2nf32$diff!=-Inf, ]),
  "FNC2_33" = nrow(test2nf33[test2nf33$diff>sqrt(8) & !is.na(test2nf33$diff) & test2nf33$diff!=Inf & test2nf33$diff!=-Inf, ]),
  "FNC2_34" = nrow(test2nf34[test2nf34$diff>sqrt(8) & !is.na(test2nf34$diff) & test2nf34$diff!=Inf & test2nf34$diff!=-Inf, ]),
  "FNC2_35" = nrow(test2nf35[test2nf35$diff>sqrt(8) & !is.na(test2nf35$diff) & test2nf35$diff!=Inf & test2nf35$diff!=-Inf, ]),
  "FNC2_36" = nrow(test2nf36[test2nf36$diff>sqrt(8) & !is.na(test2nf36$diff) & test2nf36$diff!=Inf & test2nf36$diff!=-Inf, ]),
  
  "FNC3_01" = nrow(test3nf01[test3nf01$diff>sqrt(8) & !is.na(test3nf01$diff) & test3nf01$diff!=Inf & test3nf01$diff!=-Inf, ]),
  "FNC3_02" = nrow(test3nf02[test3nf02$diff>sqrt(8) & !is.na(test3nf02$diff) & test3nf02$diff!=Inf & test3nf02$diff!=-Inf, ]),
  "FNC3_03" = nrow(test3nf03[test3nf03$diff>sqrt(8) & !is.na(test3nf03$diff) & test3nf03$diff!=Inf & test3nf03$diff!=-Inf, ]),
  "FNC3_04" = nrow(test3nf04[test3nf04$diff>sqrt(8) & !is.na(test3nf04$diff) & test3nf04$diff!=Inf & test3nf04$diff!=-Inf, ]),
  "FNC3_05" = nrow(test3nf05[test3nf05$diff>sqrt(8) & !is.na(test3nf05$diff) & test3nf05$diff!=Inf & test3nf05$diff!=-Inf, ]),
  "FNC3_06" = nrow(test3nf06[test3nf06$diff>sqrt(8) & !is.na(test3nf06$diff) & test3nf06$diff!=Inf & test3nf06$diff!=-Inf, ]),
  "FNC3_07" = nrow(test3nf07[test3nf07$diff>sqrt(8) & !is.na(test3nf07$diff) & test3nf07$diff!=Inf & test3nf07$diff!=-Inf, ]),
  "FNC3_08" = nrow(test3nf08[test3nf08$diff>sqrt(8) & !is.na(test3nf08$diff) & test3nf08$diff!=Inf & test3nf08$diff!=-Inf, ]),
  "FNC3_09" = nrow(test3nf09[test3nf09$diff>sqrt(8) & !is.na(test3nf09$diff) & test3nf09$diff!=Inf & test3nf09$diff!=-Inf, ]),
  "FNC3_10" = nrow(test3nf10[test3nf10$diff>sqrt(8) & !is.na(test3nf10$diff) & test3nf10$diff!=Inf & test3nf10$diff!=-Inf, ]),
  "FNC3_11" = nrow(test3nf11[test3nf11$diff>sqrt(8) & !is.na(test3nf11$diff) & test3nf11$diff!=Inf & test3nf11$diff!=-Inf, ]),
  "FNC3_12" = nrow(test3nf12[test3nf12$diff>sqrt(8) & !is.na(test3nf12$diff) & test3nf12$diff!=Inf & test3nf12$diff!=-Inf, ]),
  "FNC3_13" = nrow(test3nf13[test3nf13$diff>sqrt(8) & !is.na(test3nf13$diff) & test3nf13$diff!=Inf & test3nf13$diff!=-Inf, ]),
  "FNC3_14" = nrow(test3nf14[test3nf14$diff>sqrt(8) & !is.na(test3nf14$diff) & test3nf14$diff!=Inf & test3nf14$diff!=-Inf, ]),
  "FNC3_15" = nrow(test3nf15[test3nf15$diff>sqrt(8) & !is.na(test3nf15$diff) & test3nf15$diff!=Inf & test3nf15$diff!=-Inf, ]),
  "FNC3_16" = nrow(test3nf16[test3nf16$diff>sqrt(8) & !is.na(test3nf16$diff) & test3nf16$diff!=Inf & test3nf16$diff!=-Inf, ]),
  "FNC3_17" = nrow(test3nf17[test3nf17$diff>sqrt(8) & !is.na(test3nf17$diff) & test3nf17$diff!=Inf & test3nf17$diff!=-Inf, ]),
  "FNC3_18" = nrow(test3nf18[test3nf18$diff>sqrt(8) & !is.na(test3nf18$diff) & test3nf18$diff!=Inf & test3nf18$diff!=-Inf, ]),
  "FNC3_19" = nrow(test3nf19[test3nf19$diff>sqrt(8) & !is.na(test3nf19$diff) & test3nf19$diff!=Inf & test3nf19$diff!=-Inf, ]),
  "FNC3_20" = nrow(test3nf20[test3nf20$diff>sqrt(8) & !is.na(test3nf20$diff) & test3nf20$diff!=Inf & test3nf20$diff!=-Inf, ]),
  "FNC3_21" = nrow(test3nf21[test3nf21$diff>sqrt(8) & !is.na(test3nf21$diff) & test3nf21$diff!=Inf & test3nf21$diff!=-Inf, ]),
  "FNC3_22" = nrow(test3nf22[test3nf22$diff>sqrt(8) & !is.na(test3nf22$diff) & test3nf22$diff!=Inf & test3nf22$diff!=-Inf, ]),
  "FNC3_23" = nrow(test3nf23[test3nf23$diff>sqrt(8) & !is.na(test3nf23$diff) & test3nf23$diff!=Inf & test3nf23$diff!=-Inf, ]),
  "FNC3_24" = nrow(test3nf24[test3nf24$diff>sqrt(8) & !is.na(test3nf24$diff) & test3nf24$diff!=Inf & test3nf24$diff!=-Inf, ]),
  "FNC3_25" = nrow(test3nf25[test3nf25$diff>sqrt(8) & !is.na(test3nf25$diff) & test3nf25$diff!=Inf & test3nf25$diff!=-Inf, ]),
  "FNC3_26" = nrow(test3nf26[test3nf26$diff>sqrt(8) & !is.na(test3nf26$diff) & test3nf26$diff!=Inf & test3nf26$diff!=-Inf, ]),
  "FNC3_27" = nrow(test3nf27[test3nf27$diff>sqrt(8) & !is.na(test3nf27$diff) & test3nf27$diff!=Inf & test3nf27$diff!=-Inf, ]),
  "FNC3_28" = nrow(test3nf28[test3nf28$diff>sqrt(8) & !is.na(test3nf28$diff) & test3nf28$diff!=Inf & test3nf28$diff!=-Inf, ]),
  "FNC3_29" = nrow(test3nf29[test3nf29$diff>sqrt(8) & !is.na(test3nf29$diff) & test3nf29$diff!=Inf & test3nf29$diff!=-Inf, ]),
  "FNC3_30" = nrow(test3nf30[test3nf30$diff>sqrt(8) & !is.na(test3nf30$diff) & test3nf30$diff!=Inf & test3nf30$diff!=-Inf, ]),
  "FNC3_31" = nrow(test3nf31[test3nf31$diff>sqrt(8) & !is.na(test3nf31$diff) & test3nf31$diff!=Inf & test3nf31$diff!=-Inf, ]),
  "FNC3_32" = nrow(test3nf32[test3nf32$diff>sqrt(8) & !is.na(test3nf32$diff) & test3nf32$diff!=Inf & test3nf32$diff!=-Inf, ]),
  "FNC3_33" = nrow(test3nf33[test3nf33$diff>sqrt(8) & !is.na(test3nf33$diff) & test3nf33$diff!=Inf & test3nf33$diff!=-Inf, ]),
  "FNC3_34" = nrow(test3nf34[test3nf34$diff>sqrt(8) & !is.na(test3nf34$diff) & test3nf34$diff!=Inf & test3nf34$diff!=-Inf, ]),
  "FNC3_35" = nrow(test3nf35[test3nf35$diff>sqrt(8) & !is.na(test3nf35$diff) & test3nf35$diff!=Inf & test3nf35$diff!=-Inf, ]),
  "FNC3_36" = nrow(test3nf36[test3nf36$diff>sqrt(8) & !is.na(test3nf36$diff) & test3nf36$diff!=Inf & test3nf36$diff!=-Inf, ]),
  
  "FNC4_01" = nrow(test4nf01[test4nf01$diff>sqrt(8) & !is.na(test4nf01$diff) & test4nf01$diff!=Inf & test4nf01$diff!=-Inf, ]),
  "FNC4_02" = nrow(test4nf02[test4nf02$diff>sqrt(8) & !is.na(test4nf02$diff) & test4nf02$diff!=Inf & test4nf02$diff!=-Inf, ]),
  "FNC4_03" = nrow(test4nf03[test4nf03$diff>sqrt(8) & !is.na(test4nf03$diff) & test4nf03$diff!=Inf & test4nf03$diff!=-Inf, ]),
  "FNC4_04" = nrow(test4nf04[test4nf04$diff>sqrt(8) & !is.na(test4nf04$diff) & test4nf04$diff!=Inf & test4nf04$diff!=-Inf, ]),
  "FNC4_05" = nrow(test4nf05[test4nf05$diff>sqrt(8) & !is.na(test4nf05$diff) & test4nf05$diff!=Inf & test4nf05$diff!=-Inf, ]),
  "FNC4_06" = nrow(test4nf06[test4nf06$diff>sqrt(8) & !is.na(test4nf06$diff) & test4nf06$diff!=Inf & test4nf06$diff!=-Inf, ]),
  "FNC4_07" = nrow(test4nf07[test4nf07$diff>sqrt(8) & !is.na(test4nf07$diff) & test4nf07$diff!=Inf & test4nf07$diff!=-Inf, ]),
  "FNC4_08" = nrow(test4nf08[test4nf08$diff>sqrt(8) & !is.na(test4nf08$diff) & test4nf08$diff!=Inf & test4nf08$diff!=-Inf, ]),
  "FNC4_09" = nrow(test4nf09[test4nf09$diff>sqrt(8) & !is.na(test4nf09$diff) & test4nf09$diff!=Inf & test4nf09$diff!=-Inf, ]),
  "FNC4_10" = nrow(test4nf10[test4nf10$diff>sqrt(8) & !is.na(test4nf10$diff) & test4nf10$diff!=Inf & test4nf10$diff!=-Inf, ]),
  "FNC4_11" = nrow(test4nf11[test4nf11$diff>sqrt(8) & !is.na(test4nf11$diff) & test4nf11$diff!=Inf & test4nf11$diff!=-Inf, ]),
  "FNC4_12" = nrow(test4nf12[test4nf12$diff>sqrt(8) & !is.na(test4nf12$diff) & test4nf12$diff!=Inf & test4nf12$diff!=-Inf, ]),
  "FNC4_13" = nrow(test4nf13[test4nf13$diff>sqrt(8) & !is.na(test4nf13$diff) & test4nf13$diff!=Inf & test4nf13$diff!=-Inf, ]),
  "FNC4_14" = nrow(test4nf14[test4nf14$diff>sqrt(8) & !is.na(test4nf14$diff) & test4nf14$diff!=Inf & test4nf14$diff!=-Inf, ]),
  "FNC4_15" = nrow(test4nf15[test4nf15$diff>sqrt(8) & !is.na(test4nf15$diff) & test4nf15$diff!=Inf & test4nf15$diff!=-Inf, ]),
  "FNC4_16" = nrow(test4nf16[test4nf16$diff>sqrt(8) & !is.na(test4nf16$diff) & test4nf16$diff!=Inf & test4nf16$diff!=-Inf, ]),
  "FNC4_17" = nrow(test4nf17[test4nf17$diff>sqrt(8) & !is.na(test4nf17$diff) & test4nf17$diff!=Inf & test4nf17$diff!=-Inf, ]),
  "FNC4_18" = nrow(test4nf18[test4nf18$diff>sqrt(8) & !is.na(test4nf18$diff) & test4nf18$diff!=Inf & test4nf18$diff!=-Inf, ]),
  "FNC4_19" = nrow(test4nf19[test4nf19$diff>sqrt(8) & !is.na(test4nf19$diff) & test4nf19$diff!=Inf & test4nf19$diff!=-Inf, ]),
  "FNC4_20" = nrow(test4nf20[test4nf20$diff>sqrt(8) & !is.na(test4nf20$diff) & test4nf20$diff!=Inf & test4nf20$diff!=-Inf, ]),
  "FNC4_21" = nrow(test4nf21[test4nf21$diff>sqrt(8) & !is.na(test4nf21$diff) & test4nf21$diff!=Inf & test4nf21$diff!=-Inf, ]),
  "FNC4_22" = nrow(test4nf22[test4nf22$diff>sqrt(8) & !is.na(test4nf22$diff) & test4nf22$diff!=Inf & test4nf22$diff!=-Inf, ]),
  "FNC4_23" = nrow(test4nf23[test4nf23$diff>sqrt(8) & !is.na(test4nf23$diff) & test4nf23$diff!=Inf & test4nf23$diff!=-Inf, ]),
  "FNC4_24" = nrow(test4nf24[test4nf24$diff>sqrt(8) & !is.na(test4nf24$diff) & test4nf24$diff!=Inf & test4nf24$diff!=-Inf, ]),
  "FNC4_25" = nrow(test4nf25[test4nf25$diff>sqrt(8) & !is.na(test4nf25$diff) & test4nf25$diff!=Inf & test4nf25$diff!=-Inf, ]),
  "FNC4_26" = nrow(test4nf26[test4nf26$diff>sqrt(8) & !is.na(test4nf26$diff) & test4nf26$diff!=Inf & test4nf26$diff!=-Inf, ]),
  "FNC4_27" = nrow(test4nf27[test4nf27$diff>sqrt(8) & !is.na(test4nf27$diff) & test4nf27$diff!=Inf & test4nf27$diff!=-Inf, ]),
  "FNC4_28" = nrow(test4nf28[test4nf28$diff>sqrt(8) & !is.na(test4nf28$diff) & test4nf28$diff!=Inf & test4nf28$diff!=-Inf, ]),
  "FNC4_29" = nrow(test4nf29[test4nf29$diff>sqrt(8) & !is.na(test4nf29$diff) & test4nf29$diff!=Inf & test4nf29$diff!=-Inf, ]),
  "FNC4_30" = nrow(test4nf30[test4nf30$diff>sqrt(8) & !is.na(test4nf30$diff) & test4nf30$diff!=Inf & test4nf30$diff!=-Inf, ]),
  "FNC4_31" = nrow(test4nf31[test4nf31$diff>sqrt(8) & !is.na(test4nf31$diff) & test4nf31$diff!=Inf & test4nf31$diff!=-Inf, ]),
  "FNC4_32" = nrow(test4nf32[test4nf32$diff>sqrt(8) & !is.na(test4nf32$diff) & test4nf32$diff!=Inf & test4nf32$diff!=-Inf, ]),
  "FNC4_33" = nrow(test4nf33[test4nf33$diff>sqrt(8) & !is.na(test4nf33$diff) & test4nf33$diff!=Inf & test4nf33$diff!=-Inf, ]),
  "FNC4_34" = nrow(test4nf34[test4nf34$diff>sqrt(8) & !is.na(test4nf34$diff) & test4nf34$diff!=Inf & test4nf34$diff!=-Inf, ]),
  "FNC4_35" = nrow(test4nf35[test4nf35$diff>sqrt(8) & !is.na(test4nf35$diff) & test4nf35$diff!=Inf & test4nf35$diff!=-Inf, ]),
  "FNC4_36" = nrow(test4nf36[test4nf36$diff>sqrt(8) & !is.na(test4nf36$diff) & test4nf36$diff!=Inf & test4nf36$diff!=-Inf, ]))
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)

write.table(outtab, file=paste0(f,"_experiments_dist_to_corner_maxdist_diffsqrt8.txt"), row.names=F, sep='\t')
}


###generating summaries for plotting:
files <- c("./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low-high2",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low-high2",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low-high2",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low-high2",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low-high2",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low-high2",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_standard",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_standard",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_standard",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_standard",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_standard",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_standard",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_standard",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_standard",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_high",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_high",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_high",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_high",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_high",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_high",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_high",
           "./experiment_summaries_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_high")


for (f in files) {
df <- read_delim(paste0(f,"_experiments_dist_to_corner_maxdist_diffsqrt8.txt"), delim="\t", trim_ws=T)
  
columns <- c("simulation replicate",
             "N1","PC1_01","PC1_02","PC1_03","PC1_04","PC1_05","PC1_06","PC1_07","PC1_08","PC1_09","PC1_10","PC1_11","PC1_12","PC1_13","PC1_14",
             "PC1_15","PC1_16","PC1_17","PC1_18","PC1_19","PC1_20","PC1_21","PC1_22","PC1_23","PC1_24","PC1_25","PC1_26","PC1_27","PC1_28","PC1_29",
             "PC1_30","PC1_31","PC1_32","PC1_33","PC1_34","PC1_35","PC1_36",
             "N2","PC2_01","PC2_02","PC2_03","PC2_04","PC2_05","PC2_06","PC2_07","PC2_08","PC2_09","PC2_10","PC2_11","PC2_12","PC2_13","PC2_14",
             "PC2_15","PC2_16","PC2_17","PC2_18","PC2_19","PC2_20","PC2_21","PC2_22","PC2_23","PC2_24","PC2_25","PC2_26","PC2_27","PC2_28","PC2_29",
             "PC2_30","PC2_31","PC2_32","PC2_33","PC2_34","PC2_35","PC2_36",
             "N3","PC3_01","PC3_02","PC3_03","PC3_04","PC3_05","PC3_06","PC3_07","PC3_08","PC3_09","PC3_10","PC3_11","PC3_12","PC3_13","PC3_14",
             "PC3_15","PC3_16","PC3_17","PC3_18","PC3_19","PC3_20","PC3_21","PC3_22","PC3_23","PC3_24","PC3_25","PC3_26","PC3_27","PC3_28","PC3_29",
             "PC3_30","PC3_31","PC3_32","PC3_33","PC3_34","PC3_35","PC3_36",
             "N4","PC4_01","PC4_02","PC4_03","PC4_04","PC4_05","PC4_06","PC4_07","PC4_08","PC4_09","PC4_10","PC4_11","PC4_12","PC4_13","PC4_14",
             "PC4_15","PC4_16","PC4_17","PC4_18","PC4_19","PC4_20","PC4_21","PC4_22","PC4_23","PC4_24","PC4_25","PC4_26","PC4_27","PC4_28","PC4_29",
             "PC4_30","PC4_31","PC4_32","PC4_33","PC4_34","PC4_35","PC4_36",
             "FNC2_01","FNC2_02","FNC2_03","FNC2_04","FNC2_05","FNC2_06","FNC2_07","FNC2_08","FNC2_09","FNC2_10","FNC2_11","FNC2_12","FNC2_13","FNC2_14",
             "FNC2_15","FNC2_16","FNC2_17","FNC2_18","FNC2_19","FNC2_20","FNC2_21","FNC2_22","FNC2_23","FNC2_24","FNC2_25","FNC2_26","FNC2_27","FNC2_28",
             "FNC2_29","FNC2_30","FNC2_31","FNC2_32","FNC2_33","FNC2_34","FNC2_35","FNC2_36",
             "FNC3_01","FNC3_02","FNC3_03","FNC3_04","FNC3_05","FNC3_06","FNC3_07","FNC3_08","FNC3_09","FNC3_10","FNC3_11","FNC3_12","FNC3_13","FNC3_14",
             "FNC3_15","FNC3_16","FNC3_17","FNC3_18","FNC3_19","FNC3_20","FNC3_21","FNC3_22","FNC3_23","FNC3_24","FNC3_25","FNC3_26","FNC3_27","FNC3_28",
             "FNC3_29","FNC3_30","FNC3_31","FNC3_32","FNC3_33","FNC3_34","FNC3_35","FNC3_36",
             "FNC4_01","FNC4_02","FNC4_03","FNC4_04","FNC4_05","FNC4_06","FNC4_07","FNC4_08","FNC4_09","FNC4_10","FNC4_11","FNC4_12","FNC4_13","FNC4_14",
             "FNC4_15","FNC4_16","FNC4_17","FNC4_18","FNC4_19","FNC4_20","FNC4_21","FNC4_22","FNC4_23","FNC4_24","FNC4_25","FNC4_26","FNC4_27","FNC4_28",
             "FNC4_29","FNC4_30","FNC4_31","FNC4_32","FNC4_33","FNC4_34","FNC4_35","FNC4_36")
df[, columns] <- lapply(columns, function(x) as.numeric(as.character(df[[x]])))


iters <- c(1:10)

for (i in iters) {
  dfa <- filter(df, `simulation replicate`==i)
  
  output<-data.frame(matrix(nrow = 9, ncol = 41))
  colnames(	output) <- c("gene flow intensity per generation", "metric", "value01", "value02", "value03", "value04", "value05", "value06", "value07", "value08", "value09", "value10", "value11", "value12", "value13", "value14", "value15", "value16", "value17", "value18", "value19", "value20", "value21", "value22", "value23", "value24", "value25", "value26", "value27", "value28", "value29", "value30", "value31", "value32", "value33", "value34", "value35", "value36", "pop_no", "qpAdm setup", "simulation replicate")
  
  c<-1
  output[c,1] <- dfa[1,3]
  output[c,2] <- "false discovery rate, 2-way"  
  output[c,3] <- (nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01>0 & dfa$FNC2_01>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01>0 & dfa$FNC2_01>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01>0 & dfa$FNC2_01==0, ])/nrow(dfa))
  output[c,4] <- (nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02>0 & dfa$FNC2_02>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02>0 & dfa$FNC2_02>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02>0 & dfa$FNC2_02==0, ])/nrow(dfa))
  output[c,5] <- (nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03>0 & dfa$FNC2_03>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03>0 & dfa$FNC2_03>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03>0 & dfa$FNC2_03==0, ])/nrow(dfa))
  output[c,6] <- (nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04>0 & dfa$FNC2_04>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04>0 & dfa$FNC2_04>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04>0 & dfa$FNC2_04==0, ])/nrow(dfa))
  output[c,7] <- (nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05>0 & dfa$FNC2_05>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05>0 & dfa$FNC2_05>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05>0 & dfa$FNC2_05==0, ])/nrow(dfa))
  output[c,8] <- (nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06>0 & dfa$FNC2_06>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06>0 & dfa$FNC2_06>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06>0 & dfa$FNC2_06==0, ])/nrow(dfa))
  output[c,9] <- (nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07>0 & dfa$FNC2_07>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07>0 & dfa$FNC2_07>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07>0 & dfa$FNC2_07==0, ])/nrow(dfa))
  output[c,10 ] <- (nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08>0 & dfa$FNC2_08>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08>0 & dfa$FNC2_08>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08>0 & dfa$FNC2_08==0, ])/nrow(dfa))
  output[c,11 ] <- (nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09>0 & dfa$FNC2_09>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09>0 & dfa$FNC2_09>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09>0 & dfa$FNC2_09==0, ])/nrow(dfa))
  output[c,12 ] <- (nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10>0 & dfa$FNC2_10>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10>0 & dfa$FNC2_10>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10>0 & dfa$FNC2_10==0, ])/nrow(dfa))
  output[c,13 ] <- (nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11>0 & dfa$FNC2_11>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11>0 & dfa$FNC2_11>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11>0 & dfa$FNC2_11==0, ])/nrow(dfa))
  output[c,14 ] <- (nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12>0 & dfa$FNC2_12>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12>0 & dfa$FNC2_12>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12>0 & dfa$FNC2_12==0, ])/nrow(dfa))
  output[c,15 ] <- (nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13>0 & dfa$FNC2_13>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13>0 & dfa$FNC2_13>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13>0 & dfa$FNC2_13==0, ])/nrow(dfa))
  output[c,16 ] <- (nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14>0 & dfa$FNC2_14>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14>0 & dfa$FNC2_14>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14>0 & dfa$FNC2_14==0, ])/nrow(dfa))
  output[c,17 ] <- (nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15>0 & dfa$FNC2_15>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15>0 & dfa$FNC2_15>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15>0 & dfa$FNC2_15==0, ])/nrow(dfa))
  output[c,18 ] <- (nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16>0 & dfa$FNC2_16>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16>0 & dfa$FNC2_16>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16>0 & dfa$FNC2_16==0, ])/nrow(dfa))
  output[c,19 ] <- (nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17>0 & dfa$FNC2_17>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17>0 & dfa$FNC2_17>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17>0 & dfa$FNC2_17==0, ])/nrow(dfa))
  output[c,20 ] <- (nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18>0 & dfa$FNC2_18>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18>0 & dfa$FNC2_18>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18>0 & dfa$FNC2_18==0, ])/nrow(dfa))
  output[c,21 ] <- (nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19>0 & dfa$FNC2_19>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19>0 & dfa$FNC2_19>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19>0 & dfa$FNC2_19==0, ])/nrow(dfa))
  output[c,22 ] <- (nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20>0 & dfa$FNC2_20>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20>0 & dfa$FNC2_20>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20>0 & dfa$FNC2_20==0, ])/nrow(dfa))
  output[c,23 ] <- (nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21>0 & dfa$FNC2_21>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21>0 & dfa$FNC2_21>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21>0 & dfa$FNC2_21==0, ])/nrow(dfa))
  output[c,24 ] <- (nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22>0 & dfa$FNC2_22>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22>0 & dfa$FNC2_22>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22>0 & dfa$FNC2_22==0, ])/nrow(dfa))
  output[c,25 ] <- (nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23>0 & dfa$FNC2_23>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23>0 & dfa$FNC2_23>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23>0 & dfa$FNC2_23==0, ])/nrow(dfa))
  output[c,26 ] <- (nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24>0 & dfa$FNC2_24>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24>0 & dfa$FNC2_24>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24>0 & dfa$FNC2_24==0, ])/nrow(dfa))
  output[c,27 ] <- (nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25>0 & dfa$FNC2_25>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25>0 & dfa$FNC2_25>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25>0 & dfa$FNC2_25==0, ])/nrow(dfa))
  output[c,28 ] <- (nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26>0 & dfa$FNC2_26>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26>0 & dfa$FNC2_26>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26>0 & dfa$FNC2_26==0, ])/nrow(dfa))
  output[c,29 ] <- (nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27>0 & dfa$FNC2_27>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27>0 & dfa$FNC2_27>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27>0 & dfa$FNC2_27==0, ])/nrow(dfa))
  output[c,30 ] <- (nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28>0 & dfa$FNC2_28>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28>0 & dfa$FNC2_28>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28>0 & dfa$FNC2_28==0, ])/nrow(dfa))
  output[c,31 ] <- (nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29>0 & dfa$FNC2_29>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29>0 & dfa$FNC2_29>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29>0 & dfa$FNC2_29==0, ])/nrow(dfa))
  output[c,32 ] <- (nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30>0 & dfa$FNC2_30>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30>0 & dfa$FNC2_30>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30>0 & dfa$FNC2_30==0, ])/nrow(dfa))
  output[c,33 ] <- (nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31>0 & dfa$FNC2_31>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31>0 & dfa$FNC2_31>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31>0 & dfa$FNC2_31==0, ])/nrow(dfa))
  output[c,34 ] <- (nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32>0 & dfa$FNC2_32>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32>0 & dfa$FNC2_32>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32>0 & dfa$FNC2_32==0, ])/nrow(dfa))
  output[c,35 ] <- (nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33>0 & dfa$FNC2_33>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33>0 & dfa$FNC2_33>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33>0 & dfa$FNC2_33==0, ])/nrow(dfa))
  output[c,36 ] <- (nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34>0 & dfa$FNC2_34>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34>0 & dfa$FNC2_34>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34>0 & dfa$FNC2_34==0, ])/nrow(dfa))
  output[c,37 ] <- (nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35>0 & dfa$FNC2_35>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35>0 & dfa$FNC2_35>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35>0 & dfa$FNC2_35==0, ])/nrow(dfa))
  output[c,38 ] <- (nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36>0 & dfa$FNC2_36>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36>0 & dfa$FNC2_36>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36>0 & dfa$FNC2_36==0, ])/nrow(dfa))
  output[c,39] <- dfa[1,1]
  output[c,40] <- dfa[1,2]
  output[c,41] <- dfa[1,4]
  output[c+1,1] <- dfa[1,3]
  output[c+1,2] <- "false discovery rate, 3-way"  
  output[c+1,3] <- (nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01>0 & dfa$FNC3_01>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01>0 & dfa$FNC3_01>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 &dfa$PC3_01>0 & dfa$FNC3_01==0, ])/nrow(dfa))
  output[c+1,4] <- (nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02>0 & dfa$FNC3_02>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02>0 & dfa$FNC3_02>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 &dfa$PC3_02>0 & dfa$FNC3_02==0, ])/nrow(dfa))
  output[c+1,5] <- (nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03>0 & dfa$FNC3_03>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03>0 & dfa$FNC3_03>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 &dfa$PC3_03>0 & dfa$FNC3_03==0, ])/nrow(dfa))
  output[c+1,6] <- (nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04>0 & dfa$FNC3_04>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04>0 & dfa$FNC3_04>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 &dfa$PC3_04>0 & dfa$FNC3_04==0, ])/nrow(dfa))
  output[c+1,7] <- (nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05>0 & dfa$FNC3_05>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05>0 & dfa$FNC3_05>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 &dfa$PC3_05>0 & dfa$FNC3_05==0, ])/nrow(dfa))
  output[c+1,8] <- (nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06>0 & dfa$FNC3_06>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06>0 & dfa$FNC3_06>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 &dfa$PC3_06>0 & dfa$FNC3_06==0, ])/nrow(dfa))
  output[c+1,9] <- (nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07>0 & dfa$FNC3_07>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07>0 & dfa$FNC3_07>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 &dfa$PC3_07>0 & dfa$FNC3_07==0, ])/nrow(dfa))
  output[c+1,10 ] <- (nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08>0 & dfa$FNC3_08>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08>0 & dfa$FNC3_08>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 &dfa$PC3_08>0 & dfa$FNC3_08==0, ])/nrow(dfa))
  output[c+1,11 ] <- (nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09>0 & dfa$FNC3_09>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09>0 & dfa$FNC3_09>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 &dfa$PC3_09>0 & dfa$FNC3_09==0, ])/nrow(dfa))
  output[c+1,12 ] <- (nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10>0 & dfa$FNC3_10>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10>0 & dfa$FNC3_10>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 &dfa$PC3_10>0 & dfa$FNC3_10==0, ])/nrow(dfa))
  output[c+1,13 ] <- (nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11>0 & dfa$FNC3_11>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11>0 & dfa$FNC3_11>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 &dfa$PC3_11>0 & dfa$FNC3_11==0, ])/nrow(dfa))
  output[c+1,14 ] <- (nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12>0 & dfa$FNC3_12>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12>0 & dfa$FNC3_12>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 &dfa$PC3_12>0 & dfa$FNC3_12==0, ])/nrow(dfa))
  output[c+1,15 ] <- (nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13>0 & dfa$FNC3_13>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13>0 & dfa$FNC3_13>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 &dfa$PC3_13>0 & dfa$FNC3_13==0, ])/nrow(dfa))
  output[c+1,16 ] <- (nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14>0 & dfa$FNC3_14>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14>0 & dfa$FNC3_14>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 &dfa$PC3_14>0 & dfa$FNC3_14==0, ])/nrow(dfa))
  output[c+1,17 ] <- (nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15>0 & dfa$FNC3_15>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15>0 & dfa$FNC3_15>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 &dfa$PC3_15>0 & dfa$FNC3_15==0, ])/nrow(dfa))
  output[c+1,18 ] <- (nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16>0 & dfa$FNC3_16>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16>0 & dfa$FNC3_16>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 &dfa$PC3_16>0 & dfa$FNC3_16==0, ])/nrow(dfa))
  output[c+1,19 ] <- (nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17>0 & dfa$FNC3_17>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17>0 & dfa$FNC3_17>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 &dfa$PC3_17>0 & dfa$FNC3_17==0, ])/nrow(dfa))
  output[c+1,20 ] <- (nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18>0 & dfa$FNC3_18>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18>0 & dfa$FNC3_18>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 &dfa$PC3_18>0 & dfa$FNC3_18==0, ])/nrow(dfa))
  output[c+1,21 ] <- (nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19>0 & dfa$FNC3_19>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19>0 & dfa$FNC3_19>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 &dfa$PC3_19>0 & dfa$FNC3_19==0, ])/nrow(dfa))
  output[c+1,22 ] <- (nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20>0 & dfa$FNC3_20>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20>0 & dfa$FNC3_20>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 &dfa$PC3_20>0 & dfa$FNC3_20==0, ])/nrow(dfa))
  output[c+1,23 ] <- (nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21>0 & dfa$FNC3_21>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21>0 & dfa$FNC3_21>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 &dfa$PC3_21>0 & dfa$FNC3_21==0, ])/nrow(dfa))
  output[c+1,24 ] <- (nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22>0 & dfa$FNC3_22>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22>0 & dfa$FNC3_22>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 &dfa$PC3_22>0 & dfa$FNC3_22==0, ])/nrow(dfa))
  output[c+1,25 ] <- (nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23>0 & dfa$FNC3_23>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23>0 & dfa$FNC3_23>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 &dfa$PC3_23>0 & dfa$FNC3_23==0, ])/nrow(dfa))
  output[c+1,26 ] <- (nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24>0 & dfa$FNC3_24>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24>0 & dfa$FNC3_24>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 &dfa$PC3_24>0 & dfa$FNC3_24==0, ])/nrow(dfa))
  output[c+1,27 ] <- (nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25>0 & dfa$FNC3_25>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25>0 & dfa$FNC3_25>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 &dfa$PC3_25>0 & dfa$FNC3_25==0, ])/nrow(dfa))
  output[c+1,28 ] <- (nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26>0 & dfa$FNC3_26>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26>0 & dfa$FNC3_26>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 &dfa$PC3_26>0 & dfa$FNC3_26==0, ])/nrow(dfa))
  output[c+1,29 ] <- (nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27>0 & dfa$FNC3_27>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27>0 & dfa$FNC3_27>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 &dfa$PC3_27>0 & dfa$FNC3_27==0, ])/nrow(dfa))
  output[c+1,30 ] <- (nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28>0 & dfa$FNC3_28>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28>0 & dfa$FNC3_28>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 &dfa$PC3_28>0 & dfa$FNC3_28==0, ])/nrow(dfa))
  output[c+1,31 ] <- (nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29>0 & dfa$FNC3_29>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29>0 & dfa$FNC3_29>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 &dfa$PC3_29>0 & dfa$FNC3_29==0, ])/nrow(dfa))
  output[c+1,32 ] <- (nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30>0 & dfa$FNC3_30>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30>0 & dfa$FNC3_30>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 &dfa$PC3_30>0 & dfa$FNC3_30==0, ])/nrow(dfa))
  output[c+1,33 ] <- (nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31>0 & dfa$FNC3_31>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31>0 & dfa$FNC3_31>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 &dfa$PC3_31>0 & dfa$FNC3_31==0, ])/nrow(dfa))
  output[c+1,34 ] <- (nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32>0 & dfa$FNC3_32>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32>0 & dfa$FNC3_32>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 &dfa$PC3_32>0 & dfa$FNC3_32==0, ])/nrow(dfa))
  output[c+1,35 ] <- (nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33>0 & dfa$FNC3_33>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33>0 & dfa$FNC3_33>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 &dfa$PC3_33>0 & dfa$FNC3_33==0, ])/nrow(dfa))
  output[c+1,36 ] <- (nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34>0 & dfa$FNC3_34>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34>0 & dfa$FNC3_34>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 &dfa$PC3_34>0 & dfa$FNC3_34==0, ])/nrow(dfa))
  output[c+1,37 ] <- (nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35>0 & dfa$FNC3_35>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35>0 & dfa$FNC3_35>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 &dfa$PC3_35>0 & dfa$FNC3_35==0, ])/nrow(dfa))
  output[c+1,38 ] <- (nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36>0 & dfa$FNC3_36>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36>0 & dfa$FNC3_36>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 &dfa$PC3_36>0 & dfa$FNC3_36==0, ])/nrow(dfa))
  output[c+1,39] <- dfa[1,1]
  output[c+1,40] <- dfa[1,2]
  output[c+1,41] <- dfa[1,4]
  output[c+2,1] <- dfa[1,3]
  output[c+2,2] <- "false discovery rate, 4-way"  
  output[c+2,3] <- (nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01==0 & dfa$PC4_01>0 & dfa$FNC4_01>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01==0 & dfa$PC4_01>0 & dfa$FNC4_01>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01==0 & dfa$PC4_01>0 & dfa$FNC4_01==0, ])/nrow(dfa))
  output[c+2,4] <- (nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02==0 & dfa$PC4_02>0 & dfa$FNC4_02>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02==0 & dfa$PC4_02>0 & dfa$FNC4_02>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02==0 & dfa$PC4_02>0 & dfa$FNC4_02==0, ])/nrow(dfa))
  output[c+2,5] <- (nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03==0 & dfa$PC4_03>0 & dfa$FNC4_03>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03==0 & dfa$PC4_03>0 & dfa$FNC4_03>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03==0 & dfa$PC4_03>0 & dfa$FNC4_03==0, ])/nrow(dfa))
  output[c+2,6] <- (nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04==0 & dfa$PC4_04>0 & dfa$FNC4_04>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04==0 & dfa$PC4_04>0 & dfa$FNC4_04>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04==0 & dfa$PC4_04>0 & dfa$FNC4_04==0, ])/nrow(dfa))
  output[c+2,7] <- (nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05==0 & dfa$PC4_05>0 & dfa$FNC4_05>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05==0 & dfa$PC4_05>0 & dfa$FNC4_05>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05==0 & dfa$PC4_05>0 & dfa$FNC4_05==0, ])/nrow(dfa))
  output[c+2,8] <- (nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06==0 & dfa$PC4_06>0 & dfa$FNC4_06>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06==0 & dfa$PC4_06>0 & dfa$FNC4_06>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06==0 & dfa$PC4_06>0 & dfa$FNC4_06==0, ])/nrow(dfa))
  output[c+2,9] <- (nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07==0 & dfa$PC4_07>0 & dfa$FNC4_07>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07==0 & dfa$PC4_07>0 & dfa$FNC4_07>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07==0 & dfa$PC4_07>0 & dfa$FNC4_07==0, ])/nrow(dfa))
  output[c+2,10] <- (nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08==0 & dfa$PC4_08>0 & dfa$FNC4_08>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08==0 & dfa$PC4_08>0 & dfa$FNC4_08>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08==0 & dfa$PC4_08>0 & dfa$FNC4_08==0, ])/nrow(dfa))
  output[c+2,11] <- (nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09==0 & dfa$PC4_09>0 & dfa$FNC4_09>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09==0 & dfa$PC4_09>0 & dfa$FNC4_09>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09==0 & dfa$PC4_09>0 & dfa$FNC4_09==0, ])/nrow(dfa))
  output[c+2,12] <- (nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10==0 & dfa$PC4_10>0 & dfa$FNC4_10>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10==0 & dfa$PC4_10>0 & dfa$FNC4_10>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10==0 & dfa$PC4_10>0 & dfa$FNC4_10==0, ])/nrow(dfa))
  output[c+2,13] <- (nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11==0 & dfa$PC4_11>0 & dfa$FNC4_11>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11==0 & dfa$PC4_11>0 & dfa$FNC4_11>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11==0 & dfa$PC4_11>0 & dfa$FNC4_11==0, ])/nrow(dfa))
  output[c+2,14] <- (nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12==0 & dfa$PC4_12>0 & dfa$FNC4_12>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12==0 & dfa$PC4_12>0 & dfa$FNC4_12>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12==0 & dfa$PC4_12>0 & dfa$FNC4_12==0, ])/nrow(dfa))
  output[c+2,15] <- (nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13==0 & dfa$PC4_13>0 & dfa$FNC4_13>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13==0 & dfa$PC4_13>0 & dfa$FNC4_13>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13==0 & dfa$PC4_13>0 & dfa$FNC4_13==0, ])/nrow(dfa))
  output[c+2,16] <- (nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14==0 & dfa$PC4_14>0 & dfa$FNC4_14>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14==0 & dfa$PC4_14>0 & dfa$FNC4_14>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14==0 & dfa$PC4_14>0 & dfa$FNC4_14==0, ])/nrow(dfa))
  output[c+2,17] <- (nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15==0 & dfa$PC4_15>0 & dfa$FNC4_15>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15==0 & dfa$PC4_15>0 & dfa$FNC4_15>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15==0 & dfa$PC4_15>0 & dfa$FNC4_15==0, ])/nrow(dfa))
  output[c+2,18] <- (nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16==0 & dfa$PC4_16>0 & dfa$FNC4_16>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16==0 & dfa$PC4_16>0 & dfa$FNC4_16>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16==0 & dfa$PC4_16>0 & dfa$FNC4_16==0, ])/nrow(dfa))
  output[c+2,19] <- (nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17==0 & dfa$PC4_17>0 & dfa$FNC4_17>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17==0 & dfa$PC4_17>0 & dfa$FNC4_17>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17==0 & dfa$PC4_17>0 & dfa$FNC4_17==0, ])/nrow(dfa))
  output[c+2,20] <- (nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18==0 & dfa$PC4_18>0 & dfa$FNC4_18>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18==0 & dfa$PC4_18>0 & dfa$FNC4_18>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18==0 & dfa$PC4_18>0 & dfa$FNC4_18==0, ])/nrow(dfa))
  output[c+2,21] <- (nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19==0 & dfa$PC4_19>0 & dfa$FNC4_19>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19==0 & dfa$PC4_19>0 & dfa$FNC4_19>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19==0 & dfa$PC4_19>0 & dfa$FNC4_19==0, ])/nrow(dfa))
  output[c+2,22] <- (nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20==0 & dfa$PC4_20>0 & dfa$FNC4_20>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20==0 & dfa$PC4_20>0 & dfa$FNC4_20>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20==0 & dfa$PC4_20>0 & dfa$FNC4_20==0, ])/nrow(dfa))
  output[c+2,23] <- (nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21==0 & dfa$PC4_21>0 & dfa$FNC4_21>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21==0 & dfa$PC4_21>0 & dfa$FNC4_21>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21==0 & dfa$PC4_21>0 & dfa$FNC4_21==0, ])/nrow(dfa))
  output[c+2,24] <- (nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22==0 & dfa$PC4_22>0 & dfa$FNC4_22>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22==0 & dfa$PC4_22>0 & dfa$FNC4_22>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22==0 & dfa$PC4_22>0 & dfa$FNC4_22==0, ])/nrow(dfa))
  output[c+2,25] <- (nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23==0 & dfa$PC4_23>0 & dfa$FNC4_23>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23==0 & dfa$PC4_23>0 & dfa$FNC4_23>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23==0 & dfa$PC4_23>0 & dfa$FNC4_23==0, ])/nrow(dfa))
  output[c+2,26] <- (nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24==0 & dfa$PC4_24>0 & dfa$FNC4_24>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24==0 & dfa$PC4_24>0 & dfa$FNC4_24>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24==0 & dfa$PC4_24>0 & dfa$FNC4_24==0, ])/nrow(dfa))
  output[c+2,27] <- (nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25==0 & dfa$PC4_25>0 & dfa$FNC4_25>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25==0 & dfa$PC4_25>0 & dfa$FNC4_25>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25==0 & dfa$PC4_25>0 & dfa$FNC4_25==0, ])/nrow(dfa))
  output[c+2,28] <- (nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26==0 & dfa$PC4_26>0 & dfa$FNC4_26>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26==0 & dfa$PC4_26>0 & dfa$FNC4_26>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26==0 & dfa$PC4_26>0 & dfa$FNC4_26==0, ])/nrow(dfa))
  output[c+2,29] <- (nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27==0 & dfa$PC4_27>0 & dfa$FNC4_27>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27==0 & dfa$PC4_27>0 & dfa$FNC4_27>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27==0 & dfa$PC4_27>0 & dfa$FNC4_27==0, ])/nrow(dfa))
  output[c+2,30] <- (nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28==0 & dfa$PC4_28>0 & dfa$FNC4_28>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28==0 & dfa$PC4_28>0 & dfa$FNC4_28>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28==0 & dfa$PC4_28>0 & dfa$FNC4_28==0, ])/nrow(dfa))
  output[c+2,31] <- (nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29==0 & dfa$PC4_29>0 & dfa$FNC4_29>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29==0 & dfa$PC4_29>0 & dfa$FNC4_29>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29==0 & dfa$PC4_29>0 & dfa$FNC4_29==0, ])/nrow(dfa))
  output[c+2,32] <- (nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30==0 & dfa$PC4_30>0 & dfa$FNC4_30>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30==0 & dfa$PC4_30>0 & dfa$FNC4_30>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30==0 & dfa$PC4_30>0 & dfa$FNC4_30==0, ])/nrow(dfa))
  output[c+2,33] <- (nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31==0 & dfa$PC4_31>0 & dfa$FNC4_31>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31==0 & dfa$PC4_31>0 & dfa$FNC4_31>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31==0 & dfa$PC4_31>0 & dfa$FNC4_31==0, ])/nrow(dfa))
  output[c+2,34] <- (nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32==0 & dfa$PC4_32>0 & dfa$FNC4_32>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32==0 & dfa$PC4_32>0 & dfa$FNC4_32>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32==0 & dfa$PC4_32>0 & dfa$FNC4_32==0, ])/nrow(dfa))
  output[c+2,35] <- (nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33==0 & dfa$PC4_33>0 & dfa$FNC4_33>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33==0 & dfa$PC4_33>0 & dfa$FNC4_33>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33==0 & dfa$PC4_33>0 & dfa$FNC4_33==0, ])/nrow(dfa))
  output[c+2,36] <- (nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34==0 & dfa$PC4_34>0 & dfa$FNC4_34>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34==0 & dfa$PC4_34>0 & dfa$FNC4_34>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34==0 & dfa$PC4_34>0 & dfa$FNC4_34==0, ])/nrow(dfa))
  output[c+2,37] <- (nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35==0 & dfa$PC4_35>0 & dfa$FNC4_35>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35==0 & dfa$PC4_35>0 & dfa$FNC4_35>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35==0 & dfa$PC4_35>0 & dfa$FNC4_35==0, ])/nrow(dfa))
  output[c+2,38] <- (nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36==0 & dfa$PC4_36>0 & dfa$FNC4_36>0, ])/nrow(dfa))/(nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36==0 & dfa$PC4_36>0 & dfa$FNC4_36>0, ])/nrow(dfa) + nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36==0 & dfa$PC4_36>0 & dfa$FNC4_36==0, ])/nrow(dfa))
  output[c+2,39] <- dfa[1,1]
  output[c+2,40] <- dfa[1,2]
  output[c+2,41] <- dfa[1,4]
  c<-c+3
  output[c,1] <- dfa[1,3]
  output[c,2] <- "best available 2-way models NOT rejected, all simpler models rejected"
  output[c,3] <- nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01>0 & dfa$FNC2_01==0, ])/nrow(dfa)
  output[c,4] <- nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02>0 & dfa$FNC2_02==0, ])/nrow(dfa)
  output[c,5] <- nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03>0 & dfa$FNC2_03==0, ])/nrow(dfa)
  output[c,6] <- nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04>0 & dfa$FNC2_04==0, ])/nrow(dfa)
  output[c,7] <- nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05>0 & dfa$FNC2_05==0, ])/nrow(dfa)
  output[c,8] <- nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06>0 & dfa$FNC2_06==0, ])/nrow(dfa)
  output[c,9] <- nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07>0 & dfa$FNC2_07==0, ])/nrow(dfa)
  output[c,10 ] <- nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08>0 & dfa$FNC2_08==0, ])/nrow(dfa)
  output[c,11 ] <- nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09>0 & dfa$FNC2_09==0, ])/nrow(dfa)
  output[c,12 ] <- nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10>0 & dfa$FNC2_10==0, ])/nrow(dfa)
  output[c,13 ] <- nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11>0 & dfa$FNC2_11==0, ])/nrow(dfa)
  output[c,14 ] <- nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12>0 & dfa$FNC2_12==0, ])/nrow(dfa)
  output[c,15 ] <- nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13>0 & dfa$FNC2_13==0, ])/nrow(dfa)
  output[c,16 ] <- nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14>0 & dfa$FNC2_14==0, ])/nrow(dfa)
  output[c,17 ] <- nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15>0 & dfa$FNC2_15==0, ])/nrow(dfa)
  output[c,18 ] <- nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16>0 & dfa$FNC2_16==0, ])/nrow(dfa)
  output[c,19 ] <- nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17>0 & dfa$FNC2_17==0, ])/nrow(dfa)
  output[c,20 ] <- nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18>0 & dfa$FNC2_18==0, ])/nrow(dfa)
  output[c,21 ] <- nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19>0 & dfa$FNC2_19==0, ])/nrow(dfa)
  output[c,22 ] <- nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20>0 & dfa$FNC2_20==0, ])/nrow(dfa)
  output[c,23 ] <- nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21>0 & dfa$FNC2_21==0, ])/nrow(dfa)
  output[c,24 ] <- nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22>0 & dfa$FNC2_22==0, ])/nrow(dfa)
  output[c,25 ] <- nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23>0 & dfa$FNC2_23==0, ])/nrow(dfa)
  output[c,26 ] <- nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24>0 & dfa$FNC2_24==0, ])/nrow(dfa)
  output[c,27 ] <- nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25>0 & dfa$FNC2_25==0, ])/nrow(dfa)
  output[c,28 ] <- nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26>0 & dfa$FNC2_26==0, ])/nrow(dfa)
  output[c,29 ] <- nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27>0 & dfa$FNC2_27==0, ])/nrow(dfa)
  output[c,30 ] <- nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28>0 & dfa$FNC2_28==0, ])/nrow(dfa)
  output[c,31 ] <- nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29>0 & dfa$FNC2_29==0, ])/nrow(dfa)
  output[c,32 ] <- nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30>0 & dfa$FNC2_30==0, ])/nrow(dfa)
  output[c,33 ] <- nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31>0 & dfa$FNC2_31==0, ])/nrow(dfa)
  output[c,34 ] <- nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32>0 & dfa$FNC2_32==0, ])/nrow(dfa)
  output[c,35 ] <- nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33>0 & dfa$FNC2_33==0, ])/nrow(dfa)
  output[c,36 ] <- nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34>0 & dfa$FNC2_34==0, ])/nrow(dfa)
  output[c,37 ] <- nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35>0 & dfa$FNC2_35==0, ])/nrow(dfa)
  output[c,38 ] <- nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36>0 & dfa$FNC2_36==0, ])/nrow(dfa)
  output[c,39] <- dfa[1,1]
  output[c,40] <- dfa[1,2]
  output[c,41] <- dfa[1,4]
  output[c+1,1] <- dfa[1,3]
  output[c+1,2] <- "best available 3-way models NOT rejected, all simpler models rejected"
  output[c+1,3] <- nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01>0 & dfa$FNC3_01==0, ])/nrow(dfa)
  output[c+1,4] <- nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02>0 & dfa$FNC3_02==0, ])/nrow(dfa)
  output[c+1,5] <- nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03>0 & dfa$FNC3_03==0, ])/nrow(dfa)
  output[c+1,6] <- nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04>0 & dfa$FNC3_04==0, ])/nrow(dfa)
  output[c+1,7] <- nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05>0 & dfa$FNC3_05==0, ])/nrow(dfa)
  output[c+1,8] <- nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06>0 & dfa$FNC3_06==0, ])/nrow(dfa)
  output[c+1,9] <- nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07>0 & dfa$FNC3_07==0, ])/nrow(dfa)
  output[c+1,10 ] <- nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08>0 & dfa$FNC3_08==0, ])/nrow(dfa)
  output[c+1,11 ] <- nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09>0 & dfa$FNC3_09==0, ])/nrow(dfa)
  output[c+1,12 ] <- nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10>0 & dfa$FNC3_10==0, ])/nrow(dfa)
  output[c+1,13 ] <- nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11>0 & dfa$FNC3_11==0, ])/nrow(dfa)
  output[c+1,14 ] <- nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12>0 & dfa$FNC3_12==0, ])/nrow(dfa)
  output[c+1,15 ] <- nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13>0 & dfa$FNC3_13==0, ])/nrow(dfa)
  output[c+1,16 ] <- nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14>0 & dfa$FNC3_14==0, ])/nrow(dfa)
  output[c+1,17 ] <- nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15>0 & dfa$FNC3_15==0, ])/nrow(dfa)
  output[c+1,18 ] <- nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16>0 & dfa$FNC3_16==0, ])/nrow(dfa)
  output[c+1,19 ] <- nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17>0 & dfa$FNC3_17==0, ])/nrow(dfa)
  output[c+1,20 ] <- nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18>0 & dfa$FNC3_18==0, ])/nrow(dfa)
  output[c+1,21 ] <- nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19>0 & dfa$FNC3_19==0, ])/nrow(dfa)
  output[c+1,22 ] <- nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20>0 & dfa$FNC3_20==0, ])/nrow(dfa)
  output[c+1,23 ] <- nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21>0 & dfa$FNC3_21==0, ])/nrow(dfa)
  output[c+1,24 ] <- nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22>0 & dfa$FNC3_22==0, ])/nrow(dfa)
  output[c+1,25 ] <- nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23>0 & dfa$FNC3_23==0, ])/nrow(dfa)
  output[c+1,26 ] <- nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24>0 & dfa$FNC3_24==0, ])/nrow(dfa)
  output[c+1,27 ] <- nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25>0 & dfa$FNC3_25==0, ])/nrow(dfa)
  output[c+1,28 ] <- nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26>0 & dfa$FNC3_26==0, ])/nrow(dfa)
  output[c+1,29 ] <- nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27>0 & dfa$FNC3_27==0, ])/nrow(dfa)
  output[c+1,30 ] <- nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28>0 & dfa$FNC3_28==0, ])/nrow(dfa)
  output[c+1,31 ] <- nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29>0 & dfa$FNC3_29==0, ])/nrow(dfa)
  output[c+1,32 ] <- nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30>0 & dfa$FNC3_30==0, ])/nrow(dfa)
  output[c+1,33 ] <- nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31>0 & dfa$FNC3_31==0, ])/nrow(dfa)
  output[c+1,34 ] <- nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32>0 & dfa$FNC3_32==0, ])/nrow(dfa)
  output[c+1,35 ] <- nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33>0 & dfa$FNC3_33==0, ])/nrow(dfa)
  output[c+1,36 ] <- nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34>0 & dfa$FNC3_34==0, ])/nrow(dfa)
  output[c+1,37 ] <- nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35>0 & dfa$FNC3_35==0, ])/nrow(dfa)
  output[c+1,38 ] <- nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36>0 & dfa$FNC3_36==0, ])/nrow(dfa)
  output[c+1,39] <- dfa[1,1]
  output[c+1,40] <- dfa[1,2]
  output[c+1,41] <- dfa[1,4]
  output[c+2,1] <- dfa[1,3]
  output[c+2,2] <- "best available 2-way models rejected, all simpler models rejected"
  output[c+2,3] <- nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01>0 & dfa$FNC2_01>0, ])/nrow(dfa)
  output[c+2,4] <- nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02>0 & dfa$FNC2_02>0, ])/nrow(dfa)
  output[c+2,5] <- nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03>0 & dfa$FNC2_03>0, ])/nrow(dfa)
  output[c+2,6] <- nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04>0 & dfa$FNC2_04>0, ])/nrow(dfa)
  output[c+2,7] <- nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05>0 & dfa$FNC2_05>0, ])/nrow(dfa)
  output[c+2,8] <- nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06>0 & dfa$FNC2_06>0, ])/nrow(dfa)
  output[c+2,9] <- nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07>0 & dfa$FNC2_07>0, ])/nrow(dfa)
  output[c+2,10 ] <- nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08>0 & dfa$FNC2_08>0, ])/nrow(dfa)
  output[c+2,11 ] <- nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09>0 & dfa$FNC2_09>0, ])/nrow(dfa)
  output[c+2,12 ] <- nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10>0 & dfa$FNC2_10>0, ])/nrow(dfa)
  output[c+2,13 ] <- nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11>0 & dfa$FNC2_11>0, ])/nrow(dfa)
  output[c+2,14 ] <- nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12>0 & dfa$FNC2_12>0, ])/nrow(dfa)
  output[c+2,15 ] <- nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13>0 & dfa$FNC2_13>0, ])/nrow(dfa)
  output[c+2,16 ] <- nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14>0 & dfa$FNC2_14>0, ])/nrow(dfa)
  output[c+2,17 ] <- nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15>0 & dfa$FNC2_15>0, ])/nrow(dfa)
  output[c+2,18 ] <- nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16>0 & dfa$FNC2_16>0, ])/nrow(dfa)
  output[c+2,19 ] <- nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17>0 & dfa$FNC2_17>0, ])/nrow(dfa)
  output[c+2,20 ] <- nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18>0 & dfa$FNC2_18>0, ])/nrow(dfa)
  output[c+2,21 ] <- nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19>0 & dfa$FNC2_19>0, ])/nrow(dfa)
  output[c+2,22 ] <- nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20>0 & dfa$FNC2_20>0, ])/nrow(dfa)
  output[c+2,23 ] <- nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21>0 & dfa$FNC2_21>0, ])/nrow(dfa)
  output[c+2,24 ] <- nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22>0 & dfa$FNC2_22>0, ])/nrow(dfa)
  output[c+2,25 ] <- nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23>0 & dfa$FNC2_23>0, ])/nrow(dfa)
  output[c+2,26 ] <- nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24>0 & dfa$FNC2_24>0, ])/nrow(dfa)
  output[c+2,27 ] <- nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25>0 & dfa$FNC2_25>0, ])/nrow(dfa)
  output[c+2,28 ] <- nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26>0 & dfa$FNC2_26>0, ])/nrow(dfa)
  output[c+2,29 ] <- nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27>0 & dfa$FNC2_27>0, ])/nrow(dfa)
  output[c+2,30 ] <- nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28>0 & dfa$FNC2_28>0, ])/nrow(dfa)
  output[c+2,31 ] <- nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29>0 & dfa$FNC2_29>0, ])/nrow(dfa)
  output[c+2,32 ] <- nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30>0 & dfa$FNC2_30>0, ])/nrow(dfa)
  output[c+2,33 ] <- nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31>0 & dfa$FNC2_31>0, ])/nrow(dfa)
  output[c+2,34 ] <- nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32>0 & dfa$FNC2_32>0, ])/nrow(dfa)
  output[c+2,35 ] <- nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33>0 & dfa$FNC2_33>0, ])/nrow(dfa)
  output[c+2,36 ] <- nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34>0 & dfa$FNC2_34>0, ])/nrow(dfa)
  output[c+2,37 ] <- nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35>0 & dfa$FNC2_35>0, ])/nrow(dfa)
  output[c+2,38 ] <- nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36>0 & dfa$FNC2_36>0, ])/nrow(dfa)
  output[c+2,39] <- dfa[1,1]
  output[c+2,40] <- dfa[1,2]
  output[c+2,41] <- dfa[1,4]
  output[c+3,1] <- dfa[1,3]
  output[c+3,2] <- "best available 3-way models rejected, all simpler models rejected"
  output[c+3,3] <- nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01>0 & dfa$FNC3_01>0, ])/nrow(dfa)
  output[c+3,4] <- nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02>0 & dfa$FNC3_02>0, ])/nrow(dfa)
  output[c+3,5] <- nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03>0 & dfa$FNC3_03>0, ])/nrow(dfa)
  output[c+3,6] <- nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04>0 & dfa$FNC3_04>0, ])/nrow(dfa)
  output[c+3,7] <- nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05>0 & dfa$FNC3_05>0, ])/nrow(dfa)
  output[c+3,8] <- nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06>0 & dfa$FNC3_06>0, ])/nrow(dfa)
  output[c+3,9] <- nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07>0 & dfa$FNC3_07>0, ])/nrow(dfa)
  output[c+3,10 ] <- nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08>0 & dfa$FNC3_08>0, ])/nrow(dfa)
  output[c+3,11 ] <- nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09>0 & dfa$FNC3_09>0, ])/nrow(dfa)
  output[c+3,12 ] <- nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10>0 & dfa$FNC3_10>0, ])/nrow(dfa)
  output[c+3,13 ] <- nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11>0 & dfa$FNC3_11>0, ])/nrow(dfa)
  output[c+3,14 ] <- nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12>0 & dfa$FNC3_12>0, ])/nrow(dfa)
  output[c+3,15 ] <- nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13>0 & dfa$FNC3_13>0, ])/nrow(dfa)
  output[c+3,16 ] <- nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14>0 & dfa$FNC3_14>0, ])/nrow(dfa)
  output[c+3,17 ] <- nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15>0 & dfa$FNC3_15>0, ])/nrow(dfa)
  output[c+3,18 ] <- nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16>0 & dfa$FNC3_16>0, ])/nrow(dfa)
  output[c+3,19 ] <- nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17>0 & dfa$FNC3_17>0, ])/nrow(dfa)
  output[c+3,20 ] <- nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18>0 & dfa$FNC3_18>0, ])/nrow(dfa)
  output[c+3,21 ] <- nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19>0 & dfa$FNC3_19>0, ])/nrow(dfa)
  output[c+3,22 ] <- nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20>0 & dfa$FNC3_20>0, ])/nrow(dfa)
  output[c+3,23 ] <- nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21>0 & dfa$FNC3_21>0, ])/nrow(dfa)
  output[c+3,24 ] <- nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22>0 & dfa$FNC3_22>0, ])/nrow(dfa)
  output[c+3,25 ] <- nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23>0 & dfa$FNC3_23>0, ])/nrow(dfa)
  output[c+3,26 ] <- nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24>0 & dfa$FNC3_24>0, ])/nrow(dfa)
  output[c+3,27 ] <- nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25>0 & dfa$FNC3_25>0, ])/nrow(dfa)
  output[c+3,28 ] <- nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26>0 & dfa$FNC3_26>0, ])/nrow(dfa)
  output[c+3,29 ] <- nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27>0 & dfa$FNC3_27>0, ])/nrow(dfa)
  output[c+3,30 ] <- nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28>0 & dfa$FNC3_28>0, ])/nrow(dfa)
  output[c+3,31 ] <- nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29>0 & dfa$FNC3_29>0, ])/nrow(dfa)
  output[c+3,32 ] <- nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30>0 & dfa$FNC3_30>0, ])/nrow(dfa)
  output[c+3,33 ] <- nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31>0 & dfa$FNC3_31>0, ])/nrow(dfa)
  output[c+3,34 ] <- nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32>0 & dfa$FNC3_32>0, ])/nrow(dfa)
  output[c+3,35 ] <- nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33>0 & dfa$FNC3_33>0, ])/nrow(dfa)
  output[c+3,36 ] <- nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34>0 & dfa$FNC3_34>0, ])/nrow(dfa)
  output[c+3,37 ] <- nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35>0 & dfa$FNC3_35>0, ])/nrow(dfa)
  output[c+3,38 ] <- nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36>0 & dfa$FNC3_36>0, ])/nrow(dfa)
  output[c+3,39] <- dfa[1,1]
  output[c+3,40] <- dfa[1,2]
  output[c+3,41] <- dfa[1,4]
  output[c+4,1] <- dfa[1,3]
  output[c+4,2] <- "best available 4-way models NOT rejected, all simpler models rejected"
  output[c+4,3] <- nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01==0 & dfa$PC4_01>0 & dfa$FNC4_01==0, ])/nrow(dfa)
  output[c+4,4] <- nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02==0 & dfa$PC4_02>0 & dfa$FNC4_02==0, ])/nrow(dfa)
  output[c+4,5] <- nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03==0 & dfa$PC4_03>0 & dfa$FNC4_03==0, ])/nrow(dfa)
  output[c+4,6] <- nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04==0 & dfa$PC4_04>0 & dfa$FNC4_04==0, ])/nrow(dfa)
  output[c+4,7] <- nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05==0 & dfa$PC4_05>0 & dfa$FNC4_05==0, ])/nrow(dfa)
  output[c+4,8] <- nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06==0 & dfa$PC4_06>0 & dfa$FNC4_06==0, ])/nrow(dfa)
  output[c+4,9] <- nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07==0 & dfa$PC4_07>0 & dfa$FNC4_07==0, ])/nrow(dfa)
  output[c+4,10 ] <- nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08==0 & dfa$PC4_08>0 & dfa$FNC4_08==0, ])/nrow(dfa)
  output[c+4,11 ] <- nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09==0 & dfa$PC4_09>0 & dfa$FNC4_09==0, ])/nrow(dfa)
  output[c+4,12 ] <- nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10==0 & dfa$PC4_10>0 & dfa$FNC4_10==0, ])/nrow(dfa)
  output[c+4,13 ] <- nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11==0 & dfa$PC4_11>0 & dfa$FNC4_11==0, ])/nrow(dfa)
  output[c+4,14 ] <- nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12==0 & dfa$PC4_12>0 & dfa$FNC4_12==0, ])/nrow(dfa)
  output[c+4,15 ] <- nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13==0 & dfa$PC4_13>0 & dfa$FNC4_13==0, ])/nrow(dfa)
  output[c+4,16 ] <- nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14==0 & dfa$PC4_14>0 & dfa$FNC4_14==0, ])/nrow(dfa)
  output[c+4,17 ] <- nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15==0 & dfa$PC4_15>0 & dfa$FNC4_15==0, ])/nrow(dfa)
  output[c+4,18 ] <- nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16==0 & dfa$PC4_16>0 & dfa$FNC4_16==0, ])/nrow(dfa)
  output[c+4,19 ] <- nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17==0 & dfa$PC4_17>0 & dfa$FNC4_17==0, ])/nrow(dfa)
  output[c+4,20 ] <- nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18==0 & dfa$PC4_18>0 & dfa$FNC4_18==0, ])/nrow(dfa)
  output[c+4,21 ] <- nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19==0 & dfa$PC4_19>0 & dfa$FNC4_19==0, ])/nrow(dfa)
  output[c+4,22 ] <- nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20==0 & dfa$PC4_20>0 & dfa$FNC4_20==0, ])/nrow(dfa)
  output[c+4,23 ] <- nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21==0 & dfa$PC4_21>0 & dfa$FNC4_21==0, ])/nrow(dfa)
  output[c+4,24 ] <- nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22==0 & dfa$PC4_22>0 & dfa$FNC4_22==0, ])/nrow(dfa)
  output[c+4,25 ] <- nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23==0 & dfa$PC4_23>0 & dfa$FNC4_23==0, ])/nrow(dfa)
  output[c+4,26 ] <- nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24==0 & dfa$PC4_24>0 & dfa$FNC4_24==0, ])/nrow(dfa)
  output[c+4,27 ] <- nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25==0 & dfa$PC4_25>0 & dfa$FNC4_25==0, ])/nrow(dfa)
  output[c+4,28 ] <- nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26==0 & dfa$PC4_26>0 & dfa$FNC4_26==0, ])/nrow(dfa)
  output[c+4,29 ] <- nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27==0 & dfa$PC4_27>0 & dfa$FNC4_27==0, ])/nrow(dfa)
  output[c+4,30 ] <- nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28==0 & dfa$PC4_28>0 & dfa$FNC4_28==0, ])/nrow(dfa)
  output[c+4,31 ] <- nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29==0 & dfa$PC4_29>0 & dfa$FNC4_29==0, ])/nrow(dfa)
  output[c+4,32 ] <- nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30==0 & dfa$PC4_30>0 & dfa$FNC4_30==0, ])/nrow(dfa)
  output[c+4,33 ] <- nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31==0 & dfa$PC4_31>0 & dfa$FNC4_31==0, ])/nrow(dfa)
  output[c+4,34 ] <- nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32==0 & dfa$PC4_32>0 & dfa$FNC4_32==0, ])/nrow(dfa)
  output[c+4,35 ] <- nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33==0 & dfa$PC4_33>0 & dfa$FNC4_33==0, ])/nrow(dfa)
  output[c+4,36 ] <- nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34==0 & dfa$PC4_34>0 & dfa$FNC4_34==0, ])/nrow(dfa)
  output[c+4,37 ] <- nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35==0 & dfa$PC4_35>0 & dfa$FNC4_35==0, ])/nrow(dfa)
  output[c+4,38 ] <- nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36==0 & dfa$PC4_36>0 & dfa$FNC4_36==0, ])/nrow(dfa)
  output[c+4,39] <- dfa[1,1]
  output[c+4,40] <- dfa[1,2]
  output[c+4,41] <- dfa[1,4]
  output[c+5,1] <- dfa[1,3]
  output[c+5,2] <- "best available 4-way models rejected, all simpler models rejected"
  output[c+5,3] <- nrow(dfa[dfa$PC1_01==0 & dfa$PC2_01==0 & dfa$PC3_01==0 & dfa$PC4_01>0 & dfa$FNC4_01>0, ])/nrow(dfa)
  output[c+5,4] <- nrow(dfa[dfa$PC1_02==0 & dfa$PC2_02==0 & dfa$PC3_02==0 & dfa$PC4_02>0 & dfa$FNC4_02>0, ])/nrow(dfa)
  output[c+5,5] <- nrow(dfa[dfa$PC1_03==0 & dfa$PC2_03==0 & dfa$PC3_03==0 & dfa$PC4_03>0 & dfa$FNC4_03>0, ])/nrow(dfa)
  output[c+5,6] <- nrow(dfa[dfa$PC1_04==0 & dfa$PC2_04==0 & dfa$PC3_04==0 & dfa$PC4_04>0 & dfa$FNC4_04>0, ])/nrow(dfa)
  output[c+5,7] <- nrow(dfa[dfa$PC1_05==0 & dfa$PC2_05==0 & dfa$PC3_05==0 & dfa$PC4_05>0 & dfa$FNC4_05>0, ])/nrow(dfa)
  output[c+5,8] <- nrow(dfa[dfa$PC1_06==0 & dfa$PC2_06==0 & dfa$PC3_06==0 & dfa$PC4_06>0 & dfa$FNC4_06>0, ])/nrow(dfa)
  output[c+5,9] <- nrow(dfa[dfa$PC1_07==0 & dfa$PC2_07==0 & dfa$PC3_07==0 & dfa$PC4_07>0 & dfa$FNC4_07>0, ])/nrow(dfa)
  output[c+5,10 ] <- nrow(dfa[dfa$PC1_08==0 & dfa$PC2_08==0 & dfa$PC3_08==0 & dfa$PC4_08>0 & dfa$FNC4_08>0, ])/nrow(dfa)
  output[c+5,11 ] <- nrow(dfa[dfa$PC1_09==0 & dfa$PC2_09==0 & dfa$PC3_09==0 & dfa$PC4_09>0 & dfa$FNC4_09>0, ])/nrow(dfa)
  output[c+5,12 ] <- nrow(dfa[dfa$PC1_10==0 & dfa$PC2_10==0 & dfa$PC3_10==0 & dfa$PC4_10>0 & dfa$FNC4_10>0, ])/nrow(dfa)
  output[c+5,13 ] <- nrow(dfa[dfa$PC1_11==0 & dfa$PC2_11==0 & dfa$PC3_11==0 & dfa$PC4_11>0 & dfa$FNC4_11>0, ])/nrow(dfa)
  output[c+5,14 ] <- nrow(dfa[dfa$PC1_12==0 & dfa$PC2_12==0 & dfa$PC3_12==0 & dfa$PC4_12>0 & dfa$FNC4_12>0, ])/nrow(dfa)
  output[c+5,15 ] <- nrow(dfa[dfa$PC1_13==0 & dfa$PC2_13==0 & dfa$PC3_13==0 & dfa$PC4_13>0 & dfa$FNC4_13>0, ])/nrow(dfa)
  output[c+5,16 ] <- nrow(dfa[dfa$PC1_14==0 & dfa$PC2_14==0 & dfa$PC3_14==0 & dfa$PC4_14>0 & dfa$FNC4_14>0, ])/nrow(dfa)
  output[c+5,17 ] <- nrow(dfa[dfa$PC1_15==0 & dfa$PC2_15==0 & dfa$PC3_15==0 & dfa$PC4_15>0 & dfa$FNC4_15>0, ])/nrow(dfa)
  output[c+5,18 ] <- nrow(dfa[dfa$PC1_16==0 & dfa$PC2_16==0 & dfa$PC3_16==0 & dfa$PC4_16>0 & dfa$FNC4_16>0, ])/nrow(dfa)
  output[c+5,19 ] <- nrow(dfa[dfa$PC1_17==0 & dfa$PC2_17==0 & dfa$PC3_17==0 & dfa$PC4_17>0 & dfa$FNC4_17>0, ])/nrow(dfa)
  output[c+5,20 ] <- nrow(dfa[dfa$PC1_18==0 & dfa$PC2_18==0 & dfa$PC3_18==0 & dfa$PC4_18>0 & dfa$FNC4_18>0, ])/nrow(dfa)
  output[c+5,21 ] <- nrow(dfa[dfa$PC1_19==0 & dfa$PC2_19==0 & dfa$PC3_19==0 & dfa$PC4_19>0 & dfa$FNC4_19>0, ])/nrow(dfa)
  output[c+5,22 ] <- nrow(dfa[dfa$PC1_20==0 & dfa$PC2_20==0 & dfa$PC3_20==0 & dfa$PC4_20>0 & dfa$FNC4_20>0, ])/nrow(dfa)
  output[c+5,23 ] <- nrow(dfa[dfa$PC1_21==0 & dfa$PC2_21==0 & dfa$PC3_21==0 & dfa$PC4_21>0 & dfa$FNC4_21>0, ])/nrow(dfa)
  output[c+5,24 ] <- nrow(dfa[dfa$PC1_22==0 & dfa$PC2_22==0 & dfa$PC3_22==0 & dfa$PC4_22>0 & dfa$FNC4_22>0, ])/nrow(dfa)
  output[c+5,25 ] <- nrow(dfa[dfa$PC1_23==0 & dfa$PC2_23==0 & dfa$PC3_23==0 & dfa$PC4_23>0 & dfa$FNC4_23>0, ])/nrow(dfa)
  output[c+5,26 ] <- nrow(dfa[dfa$PC1_24==0 & dfa$PC2_24==0 & dfa$PC3_24==0 & dfa$PC4_24>0 & dfa$FNC4_24>0, ])/nrow(dfa)
  output[c+5,27 ] <- nrow(dfa[dfa$PC1_25==0 & dfa$PC2_25==0 & dfa$PC3_25==0 & dfa$PC4_25>0 & dfa$FNC4_25>0, ])/nrow(dfa)
  output[c+5,28 ] <- nrow(dfa[dfa$PC1_26==0 & dfa$PC2_26==0 & dfa$PC3_26==0 & dfa$PC4_26>0 & dfa$FNC4_26>0, ])/nrow(dfa)
  output[c+5,29 ] <- nrow(dfa[dfa$PC1_27==0 & dfa$PC2_27==0 & dfa$PC3_27==0 & dfa$PC4_27>0 & dfa$FNC4_27>0, ])/nrow(dfa)
  output[c+5,30 ] <- nrow(dfa[dfa$PC1_28==0 & dfa$PC2_28==0 & dfa$PC3_28==0 & dfa$PC4_28>0 & dfa$FNC4_28>0, ])/nrow(dfa)
  output[c+5,31 ] <- nrow(dfa[dfa$PC1_29==0 & dfa$PC2_29==0 & dfa$PC3_29==0 & dfa$PC4_29>0 & dfa$FNC4_29>0, ])/nrow(dfa)
  output[c+5,32 ] <- nrow(dfa[dfa$PC1_30==0 & dfa$PC2_30==0 & dfa$PC3_30==0 & dfa$PC4_30>0 & dfa$FNC4_30>0, ])/nrow(dfa)
  output[c+5,33 ] <- nrow(dfa[dfa$PC1_31==0 & dfa$PC2_31==0 & dfa$PC3_31==0 & dfa$PC4_31>0 & dfa$FNC4_31>0, ])/nrow(dfa)
  output[c+5,34 ] <- nrow(dfa[dfa$PC1_32==0 & dfa$PC2_32==0 & dfa$PC3_32==0 & dfa$PC4_32>0 & dfa$FNC4_32>0, ])/nrow(dfa)
  output[c+5,35 ] <- nrow(dfa[dfa$PC1_33==0 & dfa$PC2_33==0 & dfa$PC3_33==0 & dfa$PC4_33>0 & dfa$FNC4_33>0, ])/nrow(dfa)
  output[c+5,36 ] <- nrow(dfa[dfa$PC1_34==0 & dfa$PC2_34==0 & dfa$PC3_34==0 & dfa$PC4_34>0 & dfa$FNC4_34>0, ])/nrow(dfa)
  output[c+5,37 ] <- nrow(dfa[dfa$PC1_35==0 & dfa$PC2_35==0 & dfa$PC3_35==0 & dfa$PC4_35>0 & dfa$FNC4_35>0, ])/nrow(dfa)
  output[c+5,38 ] <- nrow(dfa[dfa$PC1_36==0 & dfa$PC2_36==0 & dfa$PC3_36==0 & dfa$PC4_36>0 & dfa$FNC4_36>0, ])/nrow(dfa)
  output[c+5,39] <- dfa[1,1]
  output[c+5,40] <- dfa[1,2]
  output[c+5,41] <- dfa[1,4]
  c<-c+6
  output <- output %>% relocate(`simulation replicate`, .after=`gene flow intensity per generation`)
  output <- output %>% relocate(`qpAdm setup`, .before=`gene flow intensity per generation`)
  output <- output %>% relocate(pop_no, .before=`qpAdm setup`)
  write.table(output, file=paste0(f,"_experiment_summaries_dist_to_corner_maxdist_diffsqrt8_iter",i,".txt"), row.names=F, sep='\t')
  rm(dfa)
  rm(output)
}
rm(df)
}


