library(readr)
library(plyr)
library(dplyr)
library(egg)
library(tidyr)
library(purrr)
library(reshape2)
library(pcds)
library(tidyverse)
library(igraph)
library(future)
library(future.apply)
library(magrittr)
future::plan(multisession)

setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm')
euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))

###re-processing of 2-way tables, 13-pop.
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_2way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_2way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_standard_2way_p2")

for (f in files) {
df2 <- read_delim(paste0(f,".txt"), delim="\t", trim_ws=T)

combs <- c(1:length(df2$target))
outtab = sapply(combs, function(x){
  if (is.na(df2$min_S1_T_S2_angle[x])) {
    c("dist_to_corner_maxdist" = NA, "dist_to_corner_avgdist" = NA)
  } else {
    c("dist_to_corner_maxdist" = euclidean_dist(c(df2$`max. distance`[x], df2$min_S1_T_S2_angle[x]/20), c(1,9)),
      "dist_to_corner_avgdist" = euclidean_dist(c(df2$`avg. distance`[x], df2$min_S1_T_S2_angle[x]/20), c(1,9)))
  } 
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("dist_to_corner_maxdist", "dist_to_corner_avgdist")
df2 <- cbind(df2, outtab)

write.table(df2, file=paste0(f,"_v2.txt"), sep='\t', row.names=F)
}



###re-processing of 2-way tables, 18-pop.
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_standard_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_high_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low-high2_2way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_standard_2way_p2")

for (f in files) {
  df2 <- read_delim(paste0(f,".txt"), delim="\t", trim_ws=T)
  
  combs <- c(1:length(df2$target))
  outtab = sapply(combs, function(x){
    if (is.na(df2$min_S1_T_S2_angle[x])) {
      c("dist_to_corner_maxdist" = NA, "dist_to_corner_avgdist" = NA)
    } else {
      c("dist_to_corner_maxdist" = euclidean_dist(c(df2$`max. distance`[x], df2$min_S1_T_S2_angle[x]/20), c(1,9)),
        "dist_to_corner_avgdist" = euclidean_dist(c(df2$`avg. distance`[x], df2$min_S1_T_S2_angle[x]/20), c(1,9)))
    } 
  }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
  colnames(outtab) <- c("dist_to_corner_maxdist", "dist_to_corner_avgdist")
  df2 <- cbind(df2, outtab)

  write.table(df2, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}



###re-processing of 3-way tables, 13 pop.
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_3way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_3way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_standard_3way_p2")

for (f in files) {
df3 <- read_delim(paste0(f,".txt"), delim="\t", trim_ws=T)

combs <- c(1:length(df3$target))
outtab = sapply(combs, function(x){
  if (is.na(df3$min_S1_T_S2_angle[x])) {
    c("dist_to_corner_maxdist" = NA, "dist_to_corner_avgdist" = NA)
  } else {
    c("dist_to_corner_maxdist" = euclidean_dist(c(df3$`max. distance`[x], df3$min_S1_T_S2_angle[x]/13.33333333), c(1,9)),
      "dist_to_corner_avgdist" = euclidean_dist(c(df3$`avg. distance`[x], df3$min_S1_T_S2_angle[x]/13.33333333), c(1,9)))
  } 
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("dist_to_corner_maxdist", "dist_to_corner_avgdist")
df3 <- cbind(df3, outtab)

write.table(df3, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}



###re-processing of 3-way tables, 18 pop.
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_standard_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_high_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low-high2_3way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_standard_3way_p2")

for (f in files) {
  df3 <- read_delim(paste0(f,".txt"), delim="\t", trim_ws=T)
  
  combs <- c(1:length(df3$target))
  outtab = sapply(combs, function(x){
    if (is.na(df3$min_S1_T_S2_angle[x])) {
      c("dist_to_corner_maxdist" = NA, "dist_to_corner_avgdist" = NA)
    } else {
      c("dist_to_corner_maxdist" = euclidean_dist(c(df3$`max. distance`[x], df3$min_S1_T_S2_angle[x]/13.33333333), c(1,9)),
        "dist_to_corner_avgdist" = euclidean_dist(c(df3$`avg. distance`[x], df3$min_S1_T_S2_angle[x]/13.33333333), c(1,9)))
    } 
  }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
  colnames(outtab) <- c("dist_to_corner_maxdist", "dist_to_corner_avgdist")
  df3 <- cbind(df3, outtab)

  write.table(df3, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}




###re-processing of 4-way tables, 13 pop.
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_13pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_low-high2_4way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_13pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_13pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_low-high2_4way_p3",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_13pop_standard_4way_p2")

for (f in files) {
  df4 <- read_delim(paste0(f,".txt"), delim="\t", trim_ws=T)
  
  combs <- c(1:length(df4$target))
  outtab = sapply(combs, function(x){
    if (is.na(df4$min_S1_T_S2_angle[x])) {
      c("dist_to_corner_maxdist" = NA, "dist_to_corner_avgdist" = NA)
    } else {
      c("dist_to_corner_maxdist" = euclidean_dist(c(df4$`max. distance`[x], df4$min_S1_T_S2_angle[x]/10), c(1,9)),
        "dist_to_corner_avgdist" = euclidean_dist(c(df4$`avg. distance`[x], df4$min_S1_T_S2_angle[x]/10), c(1,9)))
    } 
  }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
  colnames(outtab) <- c("dist_to_corner_maxdist", "dist_to_corner_avgdist")
  df4 <- cbind(df4, outtab)

  write.table(df4, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}




###re-processing of 4-way tables, 18 pop.
files <- c("./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_non-rotating_18pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_distal_rotating_18pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_non-rotating_18pop_standard_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_high_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_low-high2_4way_p2",
           "./input_files_basic_1-4way/filledcircle_qpAdm_proximal_rotating_18pop_standard_4way_p2")


for (f in files) {
  df4 <- read_delim(paste0(f,".txt"), delim="\t", trim_ws=T)
  
  combs <- c(1:length(df4$target))
  outtab = sapply(combs, function(x){
    if (is.na(df4$min_S1_T_S2_angle[x])) {
      c("dist_to_corner_maxdist" = NA, "dist_to_corner_avgdist" = NA)
    } else {
      c("dist_to_corner_maxdist" = euclidean_dist(c(df4$`max. distance`[x], df4$min_S1_T_S2_angle[x]/10), c(1,9)),
        "dist_to_corner_avgdist" = euclidean_dist(c(df4$`avg. distance`[x], df4$min_S1_T_S2_angle[x]/10), c(1,9)))
    } 
  }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
  colnames(outtab) <- c("dist_to_corner_maxdist", "dist_to_corner_avgdist")
  df4 <- cbind(df4, outtab)

  write.table(df4, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}

