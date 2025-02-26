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

df2$pop_no <- "13 demes"
df2 <- df2 %>% relocate(pop_no, .before=`qpAdm setup`)

###extract stable deme names and coordinates
dfca = as.data.frame(t(as.data.frame(strsplit(df2$target, split="_"))))[,1]
dfcb = as.data.frame(t(as.data.frame(strsplit(df2$s1, split="_"))))[,1]
dfcc = as.data.frame(t(as.data.frame(strsplit(df2$s2, split="_"))))[,1]
dfc <- as.data.frame(cbind(dfca, dfcb, dfcc))
colnames(dfc) <- c("target_no", "s1_no", "s2_no")

df <- read.table("filledcircle_coordinates.txt", sep="\t", header=T)
df$target_no <- df$deme_no
df$s1_no <- df$deme_no
df$s2_no <- df$deme_no

dfc <- inner_join(dfc, df[,c(4,2,3)], by="target_no")
dfc <- inner_join(dfc, df[,c(5,2,3)], by="s1_no")
dfc <- inner_join(dfc, df[,c(6,2,3)], by="s2_no")
colnames(dfc) <- c("target_no", "s1_no", "s2_no", "target_x", "target_y", "s1_x", "s1_y", "s2_x", "s2_y")

df2 <- cbind(df2, dfc[,c(4:9)])

###calculate_angles###
combs <- c(1:length(df2$target))
outtab = sapply(combs, function(x){
  A <- c(df2$s1_x[x], df2$s1_y[x])
  B <- c(df2$target_x[x], df2$target_y[x])
  C <- c(df2$s2_x[x], df2$s2_y[x])
  if (dist(rbind(A,B))==0 | dist(rbind(B,C))==0) {
    c("min_S1_T_S2_angle" = NA)
  } else {
    c("min_S1_T_S2_angle" = round(abs(angle.str2end(A,B,C,radian=F)$small.arc.angles[2]-angle.str2end(A,B,C,radian=F)$small.arc.angles[1]),3))
  } 
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("min_S1_T_S2_angle")
outtab[outtab == "NaN"] <- NA

df2 <- cbind(df2, outtab)
df2$`qpAdm setup` <- gsub("rotating", "rot.", df2$`qpAdm setup`)
df2$`qpAdm setup` <- gsub("proximal", "prox.", df2$`qpAdm setup`)
df2$`qpAdm setup` <- gsub("distal", "dist.", df2$`qpAdm setup`)

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
  
  df2$pop_no <- "18 demes"
  df2 <- df2 %>% relocate(pop_no, .before=`qpAdm setup`)
  
  ###extract stable deme names and coordinates
  dfca = as.data.frame(t(as.data.frame(strsplit(df2$target, split="_"))))[,1]
  dfcb = as.data.frame(t(as.data.frame(strsplit(df2$s1, split="_"))))[,1]
  dfcc = as.data.frame(t(as.data.frame(strsplit(df2$s2, split="_"))))[,1]
  dfc <- as.data.frame(cbind(dfca, dfcb, dfcc))
  colnames(dfc) <- c("target_no", "s1_no", "s2_no")
  
  df <- read.table("filledcircle_coordinates.txt", sep="\t", header=T)
  df$target_no <- df$deme_no
  df$s1_no <- df$deme_no
  df$s2_no <- df$deme_no
  
  dfc <- inner_join(dfc, df[,c(4,2,3)], by="target_no")
  dfc <- inner_join(dfc, df[,c(5,2,3)], by="s1_no")
  dfc <- inner_join(dfc, df[,c(6,2,3)], by="s2_no")
  colnames(dfc) <- c("target_no", "s1_no", "s2_no", "target_x", "target_y", "s1_x", "s1_y", "s2_x", "s2_y")
  
  df2 <- cbind(df2, dfc[,c(4:9)])
  
  ###calculate_angles###
  combs <- c(1:length(df2$target))
  outtab = sapply(combs, function(x){
    A <- c(df2$s1_x[x], df2$s1_y[x])
    B <- c(df2$target_x[x], df2$target_y[x])
    C <- c(df2$s2_x[x], df2$s2_y[x])
    if (dist(rbind(A,B))==0 | dist(rbind(B,C))==0) {
      c("min_S1_T_S2_angle" = NA)
    } else {
      c("min_S1_T_S2_angle" = round(abs(angle.str2end(A,B,C,radian=F)$small.arc.angles[2]-angle.str2end(A,B,C,radian=F)$small.arc.angles[1]),3))
    } 
  }) %>% as.data.frame() %>% set_rownames(NULL)
  colnames(outtab) <- c("min_S1_T_S2_angle")
  outtab[outtab == "NaN"] <- NA
  
  df2 <- cbind(df2, outtab)
  df2$`qpAdm setup` <- gsub("rotating", "rot.", df2$`qpAdm setup`)
  df2$`qpAdm setup` <- gsub("proximal", "prox.", df2$`qpAdm setup`)
  df2$`qpAdm setup` <- gsub("distal", "dist.", df2$`qpAdm setup`)
  
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

df3$pop_no <- "13 demes"
df3 <- df3 %>% relocate(pop_no, .before=`qpAdm setup`)

###extract stable deme names and coordinates
dfca = as.data.frame(t(as.data.frame(strsplit(df3$target, split="_"))))[,1]
dfcb = as.data.frame(t(as.data.frame(strsplit(df3$s1, split="_"))))[,1]
dfcc = as.data.frame(t(as.data.frame(strsplit(df3$s2, split="_"))))[,1]
dfcd = as.data.frame(t(as.data.frame(strsplit(df3$s3, split="_"))))[,1]
dfc <- as.data.frame(cbind(dfca, dfcb, dfcc, dfcd))
colnames(dfc) <- c("target_no", "s1_no", "s2_no", "s3_no")

df <- read.table("filledcircle_coordinates.txt", sep="\t", header=T)
df$target_no <- df$deme_no
df$s1_no <- df$deme_no
df$s2_no <- df$deme_no
df$s3_no <- df$deme_no

dfc <- inner_join(dfc, df[,c(4,2,3)], by="target_no")
dfc <- inner_join(dfc, df[,c(5,2,3)], by="s1_no")
dfc <- inner_join(dfc, df[,c(6,2,3)], by="s2_no")
dfc <- inner_join(dfc, df[,c(7,2,3)], by="s3_no")
colnames(dfc) <- c("target_no", "s1_no", "s2_no", "s3_no", "target_x", "target_y", "s1_x", "s1_y", "s2_x", "s2_y", "s3_x", "s3_y")

df3 <- cbind(df3, dfc[,5:12])

###calculate_angles###
combs <- c(1:length(df3$target))
outtab = sapply(combs, function(x){
  A <- c(df3$s1_x[x], df3$s1_y[x])
  B <- c(df3$target_x[x], df3$target_y[x])
  C <- c(df3$s2_x[x], df3$s2_y[x])
  D <- c(df3$s3_x[x], df3$s3_y[x])
  if (dist(rbind(A,B))==0 | dist(rbind(B,C))==0 | dist(rbind(B,D))==0) {
    c("min_S1_T_S2_angle" = NA)
  } else {
    c("min_S1_T_S2_angle" = round(min(abs(angle.str2end(A,B,C,radian=F)$small.arc.angles[2]-angle.str2end(A,B,C,radian=F)$small.arc.angles[1]),
                                      abs(angle.str2end(A,B,D,radian=F)$small.arc.angles[2]-angle.str2end(A,B,D,radian=F)$small.arc.angles[1]),
                                      abs(angle.str2end(C,B,D,radian=F)$small.arc.angles[2]-angle.str2end(C,B,D,radian=F)$small.arc.angles[1])),3))
  } 
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("min_S1_T_S2_angle")
outtab[outtab == "NaN"] <- NA

df3 <- cbind(df3, outtab)
df3$`qpAdm setup` <- gsub("rotating", "rot.", df3$`qpAdm setup`)
df3$`qpAdm setup` <- gsub("proximal", "prox.", df3$`qpAdm setup`)
df3$`qpAdm setup` <- gsub("distal", "dist.", df3$`qpAdm setup`)

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
  
  df3$pop_no <- "18 demes"
  df3 <- df3 %>% relocate(pop_no, .before=`qpAdm setup`)

  ###extract stable deme names and coordinates
  dfca = as.data.frame(t(as.data.frame(strsplit(df3$target, split="_"))))[,1]
  dfcb = as.data.frame(t(as.data.frame(strsplit(df3$s1, split="_"))))[,1]
  dfcc = as.data.frame(t(as.data.frame(strsplit(df3$s2, split="_"))))[,1]
  dfcd = as.data.frame(t(as.data.frame(strsplit(df3$s3, split="_"))))[,1]
  dfc <- as.data.frame(cbind(dfca, dfcb, dfcc, dfcd))
  colnames(dfc) <- c("target_no", "s1_no", "s2_no", "s3_no")
  
  df <- read.table("filledcircle_coordinates.txt", sep="\t", header=T)
  df$target_no <- df$deme_no
  df$s1_no <- df$deme_no
  df$s2_no <- df$deme_no
  df$s3_no <- df$deme_no
  
  dfc <- inner_join(dfc, df[,c(4,2,3)], by="target_no")
  dfc <- inner_join(dfc, df[,c(5,2,3)], by="s1_no")
  dfc <- inner_join(dfc, df[,c(6,2,3)], by="s2_no")
  dfc <- inner_join(dfc, df[,c(7,2,3)], by="s3_no")
  colnames(dfc) <- c("target_no", "s1_no", "s2_no", "s3_no", "target_x", "target_y", "s1_x", "s1_y", "s2_x", "s2_y", "s3_x", "s3_y")
  
  df3 <- cbind(df3, dfc[,5:12])
  
  ###calculate_angles###
  combs <- c(1:length(df3$target))
  outtab = sapply(combs, function(x){
    A <- c(df3$s1_x[x], df3$s1_y[x])
    B <- c(df3$target_x[x], df3$target_y[x])
    C <- c(df3$s2_x[x], df3$s2_y[x])
    D <- c(df3$s3_x[x], df3$s3_y[x])
    if (dist(rbind(A,B))==0 | dist(rbind(B,C))==0 | dist(rbind(B,D))==0) {
      c("min_S1_T_S2_angle" = NA)
    } else {
      c("min_S1_T_S2_angle" = round(min(abs(angle.str2end(A,B,C,radian=F)$small.arc.angles[2]-angle.str2end(A,B,C,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(A,B,D,radian=F)$small.arc.angles[2]-angle.str2end(A,B,D,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(C,B,D,radian=F)$small.arc.angles[2]-angle.str2end(C,B,D,radian=F)$small.arc.angles[1])),3))
    } 
  }) %>% as.data.frame() %>% set_rownames(NULL)
  colnames(outtab) <- c("min_S1_T_S2_angle")
  outtab[outtab == "NaN"] <- NA
  
  df3 <- cbind(df3, outtab)
  df3$`qpAdm setup` <- gsub("rotating", "rot.", df3$`qpAdm setup`)
  df3$`qpAdm setup` <- gsub("proximal", "prox.", df3$`qpAdm setup`)
  df3$`qpAdm setup` <- gsub("distal", "dist.", df3$`qpAdm setup`)
  
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
  
  df4$pop_no <- "13 demes"
  df4 <- df4 %>% relocate(pop_no, .before=`qpAdm setup`)
  
  ###extract stable deme names and coordinates
  dfca = as.data.frame(t(as.data.frame(strsplit(df4$target, split="_"))))[,1]
  dfcb = as.data.frame(t(as.data.frame(strsplit(df4$s1, split="_"))))[,1]
  dfcc = as.data.frame(t(as.data.frame(strsplit(df4$s2, split="_"))))[,1]
  dfcd = as.data.frame(t(as.data.frame(strsplit(df4$s3, split="_"))))[,1]
  dfce = as.data.frame(t(as.data.frame(strsplit(df4$s4, split="_"))))[,1]
  dfc <- as.data.frame(cbind(dfca, dfcb, dfcc, dfcd, dfce))
  colnames(dfc) <- c("target_no", "s1_no", "s2_no", "s3_no", "s4_no")
  
  df <- read.table("filledcircle_coordinates.txt", sep="\t", header=T)
  df$target_no <- df$deme_no
  df$s1_no <- df$deme_no
  df$s2_no <- df$deme_no
  df$s3_no <- df$deme_no
  df$s4_no <- df$deme_no
  
  dfc <- inner_join(dfc, df[,c(4,2,3)], by="target_no")
  dfc <- inner_join(dfc, df[,c(5,2,3)], by="s1_no")
  dfc <- inner_join(dfc, df[,c(6,2,3)], by="s2_no")
  dfc <- inner_join(dfc, df[,c(7,2,3)], by="s3_no")
  dfc <- inner_join(dfc, df[,c(8,2,3)], by="s4_no")
  colnames(dfc) <- c("target_no", "s1_no", "s2_no", "s3_no", "s4_no", "target_x", "target_y", "s1_x", "s1_y", "s2_x", "s2_y", "s3_x", "s3_y", "s4_x", "s4_y")
  
  df4 <- cbind(df4, dfc[,6:15])
  
  ###calculate_angles###
  combs <- c(1:length(df4$target))
  outtab = sapply(combs, function(x){
    A <- c(df4$s1_x[x], df4$s1_y[x])
    B <- c(df4$target_x[x], df4$target_y[x])
    C <- c(df4$s2_x[x], df4$s2_y[x])
    D <- c(df4$s3_x[x], df4$s3_y[x])
    E <- c(df4$s4_x[x], df4$s4_y[x])
    if (dist(rbind(A,B))==0 | dist(rbind(B,C))==0 | dist(rbind(B,D))==0 | dist(rbind(B,E))==0) {
      c("min_S1_T_S2_angle" = NA)
    } else {
      c("min_S1_T_S2_angle" = round(min(abs(angle.str2end(A,B,C,radian=F)$small.arc.angles[2]-angle.str2end(A,B,C,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(A,B,D,radian=F)$small.arc.angles[2]-angle.str2end(A,B,D,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(A,B,E,radian=F)$small.arc.angles[2]-angle.str2end(A,B,E,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(C,B,D,radian=F)$small.arc.angles[2]-angle.str2end(C,B,D,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(C,B,E,radian=F)$small.arc.angles[2]-angle.str2end(C,B,E,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(D,B,E,radian=F)$small.arc.angles[2]-angle.str2end(D,B,E,radian=F)$small.arc.angles[1])),3))
    } 
  }) %>% as.data.frame() %>% set_rownames(NULL)
  colnames(outtab) <- c("min_S1_T_S2_angle")
  outtab[outtab == "NaN"] <- NA
  
  df4 <- cbind(df4, outtab)
  df4$`qpAdm setup` <- gsub("rotating", "rot.", df4$`qpAdm setup`)
  df4$`qpAdm setup` <- gsub("proximal", "prox.", df4$`qpAdm setup`)
  df4$`qpAdm setup` <- gsub("distal", "dist.", df4$`qpAdm setup`)
  
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
  
  df4$pop_no <- "18 demes"
  df4 <- df4 %>% relocate(pop_no, .before=`qpAdm setup`)
  
  ###extract stable deme names and coordinates
  dfca = as.data.frame(t(as.data.frame(strsplit(df4$target, split="_"))))[,1]
  dfcb = as.data.frame(t(as.data.frame(strsplit(df4$s1, split="_"))))[,1]
  dfcc = as.data.frame(t(as.data.frame(strsplit(df4$s2, split="_"))))[,1]
  dfcd = as.data.frame(t(as.data.frame(strsplit(df4$s3, split="_"))))[,1]
  dfce = as.data.frame(t(as.data.frame(strsplit(df4$s4, split="_"))))[,1]
  dfc <- as.data.frame(cbind(dfca, dfcb, dfcc, dfcd, dfce))
  colnames(dfc) <- c("target_no", "s1_no", "s2_no", "s3_no", "s4_no")
  
  df <- read.table("filledcircle_coordinates.txt", sep="\t", header=T)
  df$target_no <- df$deme_no
  df$s1_no <- df$deme_no
  df$s2_no <- df$deme_no
  df$s3_no <- df$deme_no
  df$s4_no <- df$deme_no
  
  dfc <- inner_join(dfc, df[,c(4,2,3)], by="target_no")
  dfc <- inner_join(dfc, df[,c(5,2,3)], by="s1_no")
  dfc <- inner_join(dfc, df[,c(6,2,3)], by="s2_no")
  dfc <- inner_join(dfc, df[,c(7,2,3)], by="s3_no")
  dfc <- inner_join(dfc, df[,c(8,2,3)], by="s4_no")
  colnames(dfc) <- c("target_no", "s1_no", "s2_no", "s3_no", "s4_no", "target_x", "target_y", "s1_x", "s1_y", "s2_x", "s2_y", "s3_x", "s3_y", "s4_x", "s4_y")
  
  df4 <- cbind(df4, dfc[,6:15])
  
  ###calculate_angles###
  combs <- c(1:length(df4$target))
  outtab = sapply(combs, function(x){
    A <- c(df4$s1_x[x], df4$s1_y[x])
    B <- c(df4$target_x[x], df4$target_y[x])
    C <- c(df4$s2_x[x], df4$s2_y[x])
    D <- c(df4$s3_x[x], df4$s3_y[x])
    E <- c(df4$s4_x[x], df4$s4_y[x])
    if (dist(rbind(A,B))==0 | dist(rbind(B,C))==0 | dist(rbind(B,D))==0 | dist(rbind(B,E))==0) {
      c("min_S1_T_S2_angle" = NA)
    } else {
      c("min_S1_T_S2_angle" = round(min(abs(angle.str2end(A,B,C,radian=F)$small.arc.angles[2]-angle.str2end(A,B,C,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(A,B,D,radian=F)$small.arc.angles[2]-angle.str2end(A,B,D,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(A,B,E,radian=F)$small.arc.angles[2]-angle.str2end(A,B,E,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(C,B,D,radian=F)$small.arc.angles[2]-angle.str2end(C,B,D,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(C,B,E,radian=F)$small.arc.angles[2]-angle.str2end(C,B,E,radian=F)$small.arc.angles[1]),
                                        abs(angle.str2end(D,B,E,radian=F)$small.arc.angles[2]-angle.str2end(D,B,E,radian=F)$small.arc.angles[1])),3))
    } 
  }) %>% as.data.frame() %>% set_rownames(NULL)
  colnames(outtab) <- c("min_S1_T_S2_angle")
  outtab[outtab == "NaN"] <- NA
  
  df4 <- cbind(df4, outtab)
  df4$`qpAdm setup` <- gsub("rotating", "rot.", df4$`qpAdm setup`)
  df4$`qpAdm setup` <- gsub("proximal", "prox.", df4$`qpAdm setup`)
  df4$`qpAdm setup` <- gsub("distal", "dist.", df4$`qpAdm setup`)
  
  write.table(df4, file=paste0(f,"_v2.txt"), row.names=F, sep='\t')
}

