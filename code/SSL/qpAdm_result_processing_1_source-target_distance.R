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

setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm')

###initial processing of 2-way tables
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_2way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_2way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_standard_2way_p2")

for (f in files) {
df2 <- read.table(paste0(f,".txt"), sep="\t", header=T)
df2 <- df2[,-4]

colnames(df2) <- c("simulation replicate", "gene flow intensity per generation", "qpAdm setup", "target", "s1", "s2", "weight.s1", "se.s1", "weight.s2", "se.s2", "p_2way", "p_1way_t.s1", "p_1way_t.s2", "p_1way_s1.s2", "right")

columns2 <-c("weight.s1", "se.s1", "weight.s2", "se.s2", "p_2way", "p_1way_t.s1", "p_1way_t.s2", "p_1way_s1.s2")
df2[, columns2] <- lapply(columns2, function(x) as.numeric(as.character(df2[[x]])))

df2 <- df2 %>% relocate(`qpAdm setup`, .before=`simulation replicate`)  %>% relocate(`gene flow intensity per generation`, .before=`simulation replicate`)


###calculation of distances between demes
# load the graph and convert it to igraph format
g = read.table('filledcircle_graph.txt', header = T) %>% 
  transmute(from = pop1, to = pop2) %>% 
  igraph::graph_from_data_frame()

test <- as.data.frame(degree(g))
colnames(test) <- "degree"
median(test$degree)
mean(test$degree)
min(test$degree)


# 2-way models
targets = as.data.frame(t(as.data.frame(strsplit(df2$target, split="_"))))
s1 = as.data.frame(t(as.data.frame(strsplit(df2$s1, split="_"))))
s2 = as.data.frame(t(as.data.frame(strsplit(df2$s2, split="_"))))
models <- as.data.frame(cbind(targets$V1, s1$V1, s2$V1))
colnames(models) <- c("targets", "s1", "s2")

res1 = future_apply(models, 1, function(x) {
  out<-shortest_paths(g, from=x[1], to=x[2], mode="all", output="vpath", algorithm="unweighted")$vpath
})
res2 = future_apply(models, 1, function(x) {
  out<-shortest_paths(g, from=x[1], to=x[3], mode="all", output="vpath", algorithm="unweighted")$vpath
})


combs <- c(1:length(models$targets))
outtab2 = sapply(combs, function(comb){
  t_s1_dist = length(res1[[comb]][[1]])-1
  t_s2_dist = length(res2[[comb]][[1]])-1
  c("distance_t.s1" = t_s1_dist, "distance_t.s2" = t_s2_dist, 
    "max. distance" = max(t_s1_dist, t_s2_dist), "avg. distance" = (t_s1_dist + t_s2_dist)/2)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)

df2 <- cbind(df2, outtab2)

write.table(df2, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}




###initial processing of 3-way tables
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_3way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_3way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_standard_3way_p2")

for (f in files) {
df3 <- read.table(paste0(f,".txt"), sep="\t", header=T)
df3 <- df3[,-4]

colnames(df3) <- c("simulation replicate", "gene flow intensity per generation", "qpAdm setup", "target", "s1", "s2", "s3", "weight.s1", "se.s1", "weight.s2", "se.s2", "weight.s3", "se.s3", "p_3way", "p_1way_t.s1", "p_1way_t.s2", "p_1way_t.s3", "p_1way_s1.s2", "p_1way_s1.s3", "p_1way_s2.s3", "p_2way_t.s1.s2", "p_2way_t.s1.s3", "p_2way_t.s2.s3", "p_2way_s1.s2.s3", "right")

columns3 <-c("weight.s1", "se.s1", "weight.s2", "se.s2", "weight.s3", "se.s3", "p_3way", "p_1way_t.s1", "p_1way_t.s2", "p_1way_t.s3", "p_1way_s1.s2", "p_1way_s1.s3", "p_1way_s2.s3", "p_2way_t.s1.s2", "p_2way_t.s1.s3", "p_2way_t.s2.s3", "p_2way_s1.s2.s3")
df3[, columns3] <- lapply(columns3, function(x) as.numeric(as.character(df3[[x]])))

df3 <- df3 %>% relocate(`qpAdm setup`, .before=`simulation replicate`)  %>% relocate(`gene flow intensity per generation`, .before=`simulation replicate`)


###calculation of distances between demes
# load the graph and convert it to igraph format
g = read.table('filledcircle_graph.txt', header = T) %>% 
  transmute(from = pop1, to = pop2) %>% 
  igraph::graph_from_data_frame()

# 3-way models
targets = as.data.frame(t(as.data.frame(strsplit(df3$target, split="_"))))
s1 = as.data.frame(t(as.data.frame(strsplit(df3$s1, split="_"))))
s2 = as.data.frame(t(as.data.frame(strsplit(df3$s2, split="_"))))
s3 = as.data.frame(t(as.data.frame(strsplit(df3$s3, split="_"))))
models <- as.data.frame(cbind(targets$V1, s1$V1, s2$V1, s3$V1))
colnames(models) <- c("targets", "s1", "s2", "s3")

res1 = future_apply(models, 1, function(x) {
  out<-shortest_paths(g, from=x[1], to=x[2], mode="all", output="vpath", algorithm="unweighted")$vpath
})
res2 = future_apply(models, 1, function(x) {
  out<-shortest_paths(g, from=x[1], to=x[3], mode="all", output="vpath", algorithm="unweighted")$vpath
})
res3 = future_apply(models, 1, function(x) {
  out<-shortest_paths(g, from=x[1], to=x[4], mode="all", output="vpath", algorithm="unweighted")$vpath
})


combs <- c(1:length(models$targets))
outtab3 = sapply(combs, function(comb){
  t_s1_dist = length(res1[[comb]][[1]])-1
  t_s2_dist = length(res2[[comb]][[1]])-1
  t_s3_dist = length(res3[[comb]][[1]])-1
  c("distance_t.s1" = t_s1_dist, "distance_t.s2" = t_s2_dist, "distance_t.s3" = t_s3_dist, 
    "max. distance" = max(t_s1_dist, t_s2_dist, t_s3_dist), "avg. distance" = (t_s1_dist + t_s2_dist + t_s3_dist)/3)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)

df3 <- cbind(df3, outtab3)

write.table(df3, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}



###initial processing of 4-way tables
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_4way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_4way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_standard_4way_p2")

for (f in files) {
  df4 <- read.table(paste0(f,".txt"), sep="\t", header=T)
  df4 <- df4[,-4]
  
  colnames(df4) <- c("simulation replicate","gene flow intensity per generation","qpAdm setup","target","s1","s2","s3","s4","weight.s1","se.s1","weight.s2","se.s2","weight.s3","se.s3","weight.s4","se.s4","p_4way","p_1way_t.s1","p_1way_t.s2","p_1way_t.s3","p_1way_t.s4","p_1way_s1.s2","p_1way_s1.s3","p_1way_s1.s4","p_1way_s2.s3","p_1way_s2.s4","p_1way_s3.s4",
                     "p_2way_t.s1.s2","p_2way_t.s1.s3","p_2way_t.s1.s4","p_2way_t.s2.s3","p_2way_t.s2.s4","p_2way_t.s3.s4","p_2way_s1.s2.s3","p_2way_s1.s2.s4","p_2way_s1.s3.s4","p_2way_s2.s3.s4","p_3way_t.s1.s2.s3","p_3way_t.s1.s2.s4","p_3way_t.s1.s3.s4","p_3way_t.s2.s3.s4","p_3way_s1.s2.s3.s4","right")
  
  columns4 <-c("weight.s1","se.s1","weight.s2","se.s2","weight.s3","se.s3","weight.s4","se.s4","p_4way","p_1way_t.s1","p_1way_t.s2","p_1way_t.s3","p_1way_t.s4","p_1way_s1.s2","p_1way_s1.s3","p_1way_s1.s4","p_1way_s2.s3","p_1way_s2.s4","p_1way_s3.s4",
               "p_2way_t.s1.s2","p_2way_t.s1.s3","p_2way_t.s1.s4","p_2way_t.s2.s3","p_2way_t.s2.s4","p_2way_t.s3.s4","p_2way_s1.s2.s3","p_2way_s1.s2.s4","p_2way_s1.s3.s4","p_2way_s2.s3.s4","p_3way_t.s1.s2.s3","p_3way_t.s1.s2.s4","p_3way_t.s1.s3.s4","p_3way_t.s2.s3.s4","p_3way_s1.s2.s3.s4")
  df4[, columns4] <- lapply(columns4, function(x) as.numeric(as.character(df4[[x]])))
  
  df4 <- df4 %>% relocate(`qpAdm setup`, .before=`simulation replicate`)  %>% relocate(`gene flow intensity per generation`, .before=`simulation replicate`)
  
  
  ###calculation of distances between demes
  # load the graph and convert it to igraph format
  g = read.table('filledcircle_graph.txt', header = T) %>% 
    transmute(from = pop1, to = pop2) %>% 
    igraph::graph_from_data_frame()
  
  # 4-way models
  targets = as.data.frame(t(as.data.frame(strsplit(df4$target, split="_"))))
  s1 = as.data.frame(t(as.data.frame(strsplit(df4$s1, split="_"))))
  s2 = as.data.frame(t(as.data.frame(strsplit(df4$s2, split="_"))))
  s3 = as.data.frame(t(as.data.frame(strsplit(df4$s3, split="_"))))
  s4 = as.data.frame(t(as.data.frame(strsplit(df4$s4, split="_"))))
  models <- as.data.frame(cbind(targets$V1, s1$V1, s2$V1, s3$V1, s4$V1))
  colnames(models) <- c("targets", "s1", "s2", "s3", "s4")
  
  res1 = future_apply(models, 1, function(x) {
    out<-shortest_paths(g, from=x[1], to=x[2], mode="all", output="vpath", algorithm="unweighted")$vpath
  })
  res2 = future_apply(models, 1, function(x) {
    out<-shortest_paths(g, from=x[1], to=x[3], mode="all", output="vpath", algorithm="unweighted")$vpath
  })
  res3 = future_apply(models, 1, function(x) {
    out<-shortest_paths(g, from=x[1], to=x[4], mode="all", output="vpath", algorithm="unweighted")$vpath
  })
  res4 = future_apply(models, 1, function(x) {
    out<-shortest_paths(g, from=x[1], to=x[5], mode="all", output="vpath", algorithm="unweighted")$vpath
  })
  
  combs <- c(1:length(models$targets))
  outtab3 = sapply(combs, function(comb){
    t_s1_dist = length(res1[[comb]][[1]])-1
    t_s2_dist = length(res2[[comb]][[1]])-1
    t_s3_dist = length(res3[[comb]][[1]])-1
    t_s4_dist = length(res4[[comb]][[1]])-1
    c("distance_t.s1" = t_s1_dist, "distance_t.s2" = t_s2_dist, "distance_t.s3" = t_s3_dist, "distance_t.s4" = t_s4_dist,
      "max. distance" = max(t_s1_dist, t_s2_dist, t_s3_dist, t_s4_dist), "avg. distance" = (t_s1_dist + t_s2_dist + t_s3_dist + t_s4_dist)/4)
  }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
  
  df4 <- cbind(df4, outtab3)
  
  write.table(df4, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}

