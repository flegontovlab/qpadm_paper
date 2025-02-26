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

###modify feasibility criteria###
df2$positive01 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                           df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                           df2$p_2way >=0.001 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
df2$positive02 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                           df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                           df2$p_2way >=0.01 & df2$p_1way_t.s1 <0.01 & df2$p_1way_t.s2 <0.01 & df2$p_1way_s1.s2 <0.01, "fitting", "non-fitting")
df2$positive03 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                           df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                           df2$p_2way >=0.05 & df2$p_1way_t.s1 <0.05 & df2$p_1way_t.s2 <0.05 & df2$p_1way_s1.s2 <0.05, "fitting", "non-fitting")
df2$positive04 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                           df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                           df2$p_2way >=0.1 & df2$p_1way_t.s1 <0.1 & df2$p_1way_t.s2 <0.1 & df2$p_1way_s1.s2 <0.1, "fitting", "non-fitting")
df2$positive05 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                           df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                           df2$p_2way >=0.5 & df2$p_1way_t.s1 <0.5 & df2$p_1way_t.s2 <0.5 & df2$p_1way_s1.s2 <0.5, "fitting", "non-fitting")
df2$positive06 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                           df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 &
                           df2$p_2way >=0.05 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
df2$positive07 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                           df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 &
                           df2$p_2way >=0.5 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")

df2$positive08 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.001 
                       & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
df2$positive09 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.01 
                       & df2$p_1way_t.s1 <0.01 & df2$p_1way_t.s2 <0.01 & df2$p_1way_s1.s2 <0.01, "fitting", "non-fitting")
df2$positive10 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.05 
                       & df2$p_1way_t.s1 <0.05 & df2$p_1way_t.s2 <0.05 & df2$p_1way_s1.s2 <0.05, "fitting", "non-fitting")
df2$positive11 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.1 
                       & df2$p_1way_t.s1 <0.1 & df2$p_1way_t.s2 <0.1 & df2$p_1way_s1.s2 <0.1, "fitting", "non-fitting")
df2$positive12 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.5 
                       & df2$p_1way_t.s1 <0.5 & df2$p_1way_t.s2 <0.5 & df2$p_1way_s1.s2 <0.5, "fitting", "non-fitting")
df2$positive13 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.05 
                       & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
df2$positive14 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.5 
                       & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")

df2$positive15 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.001 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
df2$positive16 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.01 & df2$p_1way_t.s1 <0.01 & df2$p_1way_t.s2 <0.01 & df2$p_1way_s1.s2 <0.01, "fitting", "non-fitting")
df2$positive17 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.05 & df2$p_1way_t.s1 <0.05 & df2$p_1way_t.s2 <0.05 & df2$p_1way_s1.s2 <0.05, "fitting", "non-fitting")
df2$positive18 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.1 & df2$p_1way_t.s1 <0.1 & df2$p_1way_t.s2 <0.1 & df2$p_1way_s1.s2 <0.1, "fitting", "non-fitting")
df2$positive19 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.5 & df2$p_1way_t.s1 <0.5 & df2$p_1way_t.s2 <0.5 & df2$p_1way_s1.s2 <0.5, "fitting", "non-fitting")
df2$positive20 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.05 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
df2$positive21 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.5 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")

df2$positive22 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.001, "fitting", "non-fitting")
df2$positive23 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.01, "fitting", "non-fitting")
df2$positive24 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.05, "fitting", "non-fitting")
df2$positive25 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.1, "fitting", "non-fitting")
df2$positive26 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.5, "fitting", "non-fitting")

df2$positive27 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.001, "fitting", "non-fitting")
df2$positive28 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.01, "fitting", "non-fitting")
df2$positive29 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.05, "fitting", "non-fitting")
df2$positive30 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.1, "fitting", "non-fitting")
df2$positive31 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.5, "fitting", "non-fitting")

df2$positive32 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.001, "fitting", "non-fitting")
df2$positive33 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.01, "fitting", "non-fitting")
df2$positive34 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.05, "fitting", "non-fitting")
df2$positive35 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.1, "fitting", "non-fitting")
df2$positive36 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.5, "fitting", "non-fitting")

write.table(df2, file=paste0(f,"_correct.txt"), sep='\t', row.names=F)
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
  
  ###modify feasibility criteria###
  df2$positive01 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                             df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                             df2$p_2way >=0.001 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  df2$positive02 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                             df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                             df2$p_2way >=0.01 & df2$p_1way_t.s1 <0.01 & df2$p_1way_t.s2 <0.01 & df2$p_1way_s1.s2 <0.01, "fitting", "non-fitting")
  df2$positive03 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                             df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                             df2$p_2way >=0.05 & df2$p_1way_t.s1 <0.05 & df2$p_1way_t.s2 <0.05 & df2$p_1way_s1.s2 <0.05, "fitting", "non-fitting")
  df2$positive04 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                             df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                             df2$p_2way >=0.1 & df2$p_1way_t.s1 <0.1 & df2$p_1way_t.s2 <0.1 & df2$p_1way_s1.s2 <0.1, "fitting", "non-fitting")
  df2$positive05 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                             df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & 
                             df2$p_2way >=0.5 & df2$p_1way_t.s1 <0.5 & df2$p_1way_t.s2 <0.5 & df2$p_1way_s1.s2 <0.5, "fitting", "non-fitting")
  df2$positive06 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                             df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 &
                             df2$p_2way >=0.05 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  df2$positive07 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & 
                             df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 &
                             df2$p_2way >=0.5 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  
  df2$positive08 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.001 
                           & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  df2$positive09 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.01 
                           & df2$p_1way_t.s1 <0.01 & df2$p_1way_t.s2 <0.01 & df2$p_1way_s1.s2 <0.01, "fitting", "non-fitting")
  df2$positive10 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.05 
                           & df2$p_1way_t.s1 <0.05 & df2$p_1way_t.s2 <0.05 & df2$p_1way_s1.s2 <0.05, "fitting", "non-fitting")
  df2$positive11 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.1 
                           & df2$p_1way_t.s1 <0.1 & df2$p_1way_t.s2 <0.1 & df2$p_1way_s1.s2 <0.1, "fitting", "non-fitting")
  df2$positive12 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.5 
                           & df2$p_1way_t.s1 <0.5 & df2$p_1way_t.s2 <0.5 & df2$p_1way_s1.s2 <0.5, "fitting", "non-fitting")
  df2$positive13 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.05 
                           & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  df2$positive14 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.5 
                           & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  
  df2$positive15 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.001 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  df2$positive16 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.01 & df2$p_1way_t.s1 <0.01 & df2$p_1way_t.s2 <0.01 & df2$p_1way_s1.s2 <0.01, "fitting", "non-fitting")
  df2$positive17 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.05 & df2$p_1way_t.s1 <0.05 & df2$p_1way_t.s2 <0.05 & df2$p_1way_s1.s2 <0.05, "fitting", "non-fitting")
  df2$positive18 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.1 & df2$p_1way_t.s1 <0.1 & df2$p_1way_t.s2 <0.1 & df2$p_1way_s1.s2 <0.1, "fitting", "non-fitting")
  df2$positive19 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.5 & df2$p_1way_t.s1 <0.5 & df2$p_1way_t.s2 <0.5 & df2$p_1way_s1.s2 <0.5, "fitting", "non-fitting")
  df2$positive20 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.05 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  df2$positive21 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.5 & df2$p_1way_t.s1 <0.001 & df2$p_1way_t.s2 <0.001 & df2$p_1way_s1.s2 <0.001, "fitting", "non-fitting")
  
  df2$positive22 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.001, "fitting", "non-fitting")
  df2$positive23 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.01, "fitting", "non-fitting")
  df2$positive24 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.05, "fitting", "non-fitting")
  df2$positive25 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.1, "fitting", "non-fitting")
  df2$positive26 <- ifelse(df2$weight.s1-2*df2$se.s1 >0 & df2$weight.s1+2*df2$se.s1 >0 & df2$weight.s2-2*df2$se.s2 >0 & df2$weight.s2+2*df2$se.s2 >0 & df2$weight.s1-2*df2$se.s1 <1 & df2$weight.s1+2*df2$se.s1 <1 & df2$weight.s2-2*df2$se.s2 <1 & df2$weight.s2+2*df2$se.s2 <1 & df2$p_2way >=0.5, "fitting", "non-fitting")
  
  df2$positive27 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.001, "fitting", "non-fitting")
  df2$positive28 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.01, "fitting", "non-fitting")
  df2$positive29 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.05, "fitting", "non-fitting")
  df2$positive30 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.1, "fitting", "non-fitting")
  df2$positive31 <- ifelse(df2$weight.s1<1 & df2$weight.s2<1 & df2$weight.s1>0 & df2$weight.s2>0 & df2$p_2way >=0.5, "fitting", "non-fitting")
  
  df2$positive32 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.001, "fitting", "non-fitting")
  df2$positive33 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.01, "fitting", "non-fitting")
  df2$positive34 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.05, "fitting", "non-fitting")
  df2$positive35 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.1, "fitting", "non-fitting")
  df2$positive36 <- ifelse(df2$weight.s1<1.3 & df2$weight.s2<1.3 & df2$weight.s1>-0.3 & df2$weight.s2>-0.3 & df2$p_2way >=0.5, "fitting", "non-fitting")
  
  write.table(df2, file=paste0(f,"_correct.txt"), sep='\t', row.names=F)
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

###add feasibility criteria###
df3$positive01 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                         df3$p_3way >=0.001 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 &  
                         df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
df3$positive02 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                         df3$p_3way >=0.01 & df3$p_1way_t.s1 <0.01 & df3$p_1way_t.s2 <0.01 & df3$p_1way_t.s3 <0.01 & df3$p_1way_s1.s2 <0.01 & df3$p_1way_s1.s3 <0.01 & df3$p_1way_s2.s3 <0.01 &  
                         df3$p_2way_t.s1.s2 <0.01 & df3$p_2way_t.s1.s3 <0.01 & df3$p_2way_t.s2.s3 <0.01 & df3$p_2way_s1.s2.s3 <0.01, "fitting", "non-fitting")
df3$positive03 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                         df3$p_3way >=0.05 & df3$p_1way_t.s1 <0.05 & df3$p_1way_t.s2 <0.05 & df3$p_1way_t.s3 <0.05 & df3$p_1way_s1.s2 <0.05 & df3$p_1way_s1.s3 <0.05 & df3$p_1way_s2.s3 <0.05 &  
                         df3$p_2way_t.s1.s2 <0.05 & df3$p_2way_t.s1.s3 <0.05 & df3$p_2way_t.s2.s3 <0.05 & df3$p_2way_s1.s2.s3 <0.05, "fitting", "non-fitting")
df3$positive04 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                         df3$p_3way >=0.1 & df3$p_1way_t.s1 <0.1 & df3$p_1way_t.s2 <0.1 & df3$p_1way_t.s3 <0.1 & df3$p_1way_s1.s2 <0.1 & df3$p_1way_s1.s3 <0.1 & df3$p_1way_s2.s3 <0.1 &  
                         df3$p_2way_t.s1.s2 <0.1 & df3$p_2way_t.s1.s3 <0.1 & df3$p_2way_t.s2.s3 <0.1 & df3$p_2way_s1.s2.s3 <0.1, "fitting", "non-fitting")
df3$positive05 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                         df3$p_3way >=0.5 & df3$p_1way_t.s1 <0.5 & df3$p_1way_t.s2 <0.5 & df3$p_1way_t.s3 <0.5 & df3$p_1way_s1.s2 <0.5 & df3$p_1way_s1.s3 <0.5 & df3$p_1way_s2.s3 <0.5 &  
                         df3$p_2way_t.s1.s2 <0.5 & df3$p_2way_t.s1.s3 <0.5 & df3$p_2way_t.s2.s3 <0.5 & df3$p_2way_s1.s2.s3 <0.5, "fitting", "non-fitting")
df3$positive06 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                         df3$p_3way >=0.05 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 &  
                         df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
df3$positive07 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                         df3$p_3way >=0.5 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 &  
                         df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")

df3$positive08 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.001 
                       & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                       & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
df3$positive09 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.01 
                         & df3$p_1way_t.s1 <0.01 & df3$p_1way_t.s2 <0.01 & df3$p_1way_t.s3 <0.01 & df3$p_1way_s1.s2 <0.01 & df3$p_1way_s1.s3 <0.01 & df3$p_1way_s2.s3 <0.01 
                         & df3$p_2way_t.s1.s2 <0.01 & df3$p_2way_t.s1.s3 <0.01 & df3$p_2way_t.s2.s3 <0.01 & df3$p_2way_s1.s2.s3 <0.01, "fitting", "non-fitting")
df3$positive10 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.05 
                         & df3$p_1way_t.s1 <0.05 & df3$p_1way_t.s2 <0.05 & df3$p_1way_t.s3 <0.05 & df3$p_1way_s1.s2 <0.05 & df3$p_1way_s1.s3 <0.05 & df3$p_1way_s2.s3 <0.05 
                         & df3$p_2way_t.s1.s2 <0.05 & df3$p_2way_t.s1.s3 <0.05 & df3$p_2way_t.s2.s3 <0.05 & df3$p_2way_s1.s2.s3 <0.05, "fitting", "non-fitting")
df3$positive11 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.1 
                       & df3$p_1way_t.s1 <0.1 & df3$p_1way_t.s2 <0.1 & df3$p_1way_t.s3 <0.1 & df3$p_1way_s1.s2 <0.1 & df3$p_1way_s1.s3 <0.1 & df3$p_1way_s2.s3 <0.1 
                       & df3$p_2way_t.s1.s2 <0.1 & df3$p_2way_t.s1.s3 <0.1 & df3$p_2way_t.s2.s3 <0.1 & df3$p_2way_s1.s2.s3 <0.1, "fitting", "non-fitting")
df3$positive12 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.5 
                       & df3$p_1way_t.s1 <0.5 & df3$p_1way_t.s2 <0.5 & df3$p_1way_t.s3 <0.5 & df3$p_1way_s1.s2 <0.5 & df3$p_1way_s1.s3 <0.5 & df3$p_1way_s2.s3 <0.5 
                       & df3$p_2way_t.s1.s2 <0.5 & df3$p_2way_t.s1.s3 <0.5 & df3$p_2way_t.s2.s3 <0.5 & df3$p_2way_s1.s2.s3 <0.5, "fitting", "non-fitting")
df3$positive13 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.05 
                       & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                       & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
df3$positive14 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 &df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 &  df3$p_3way >=0.5 
                       & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                       & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")

df3$positive15 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.001 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
df3$positive16 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.01 & df3$p_1way_t.s1 <0.01 & df3$p_1way_t.s2 <0.01 & df3$p_1way_t.s3 <0.01 & df3$p_1way_s1.s2 <0.01 & df3$p_1way_s1.s3 <0.01 & df3$p_1way_s2.s3 <0.01 
                           & df3$p_2way_t.s1.s2 <0.01 & df3$p_2way_t.s1.s3 <0.01 & df3$p_2way_t.s2.s3 <0.01 & df3$p_2way_s1.s2.s3 <0.01, "fitting", "non-fitting")
df3$positive17 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.05 & df3$p_1way_t.s1 <0.05 & df3$p_1way_t.s2 <0.05 & df3$p_1way_t.s3 <0.05 & df3$p_1way_s1.s2 <0.05 & df3$p_1way_s1.s3 <0.05 & df3$p_1way_s2.s3 <0.05 
                           & df3$p_2way_t.s1.s2 <0.05 & df3$p_2way_t.s1.s3 <0.05 & df3$p_2way_t.s2.s3 <0.05 & df3$p_2way_s1.s2.s3 <0.05, "fitting", "non-fitting")
df3$positive18 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.1 & df3$p_1way_t.s1 <0.1 & df3$p_1way_t.s2 <0.1 & df3$p_1way_t.s3 <0.1 & df3$p_1way_s1.s2 <0.1 & df3$p_1way_s1.s3 <0.1 & df3$p_1way_s2.s3 <0.1 
                           & df3$p_2way_t.s1.s2 <0.1 & df3$p_2way_t.s1.s3 <0.1 & df3$p_2way_t.s2.s3 <0.1 & df3$p_2way_s1.s2.s3 <0.1, "fitting", "non-fitting")
df3$positive19 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.5 & df3$p_1way_t.s1 <0.5 & df3$p_1way_t.s2 <0.5 & df3$p_1way_t.s3 <0.5 & df3$p_1way_s1.s2 <0.5 & df3$p_1way_s1.s3 <0.5 & df3$p_1way_s2.s3 <0.5 
                           & df3$p_2way_t.s1.s2 <0.5 & df3$p_2way_t.s1.s3 <0.5 & df3$p_2way_t.s2.s3 <0.5 & df3$p_2way_s1.s2.s3 <0.5, "fitting", "non-fitting")
df3$positive20 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.05 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
df3$positive21 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.5 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")

df3$positive22 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.001, "fitting", "non-fitting")
df3$positive23 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.01, "fitting", "non-fitting")
df3$positive24 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.05, "fitting", "non-fitting")
df3$positive25 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.1, "fitting", "non-fitting")
df3$positive26 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                         df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.5, "fitting", "non-fitting")

df3$positive27 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.001, "fitting", "non-fitting")
df3$positive28 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.01, "fitting", "non-fitting")
df3$positive29 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.05, "fitting", "non-fitting")
df3$positive30 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.1, "fitting", "non-fitting")
df3$positive31 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.5, "fitting", "non-fitting")

df3$positive32 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.001, "fitting", "non-fitting")
df3$positive33 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.01, "fitting", "non-fitting")
df3$positive34 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.05, "fitting", "non-fitting")
df3$positive35 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.1, "fitting", "non-fitting")
df3$positive36 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.5, "fitting", "non-fitting")

write.table(df3, file=paste0(f,"_correct.txt"), row.names=F, sep='\t')
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
  
  ###add feasibility criteria###
  df3$positive01 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                             df3$p_3way >=0.001 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 &  
                             df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  df3$positive02 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                             df3$p_3way >=0.01 & df3$p_1way_t.s1 <0.01 & df3$p_1way_t.s2 <0.01 & df3$p_1way_t.s3 <0.01 & df3$p_1way_s1.s2 <0.01 & df3$p_1way_s1.s3 <0.01 & df3$p_1way_s2.s3 <0.01 &  
                             df3$p_2way_t.s1.s2 <0.01 & df3$p_2way_t.s1.s3 <0.01 & df3$p_2way_t.s2.s3 <0.01 & df3$p_2way_s1.s2.s3 <0.01, "fitting", "non-fitting")
  df3$positive03 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                             df3$p_3way >=0.05 & df3$p_1way_t.s1 <0.05 & df3$p_1way_t.s2 <0.05 & df3$p_1way_t.s3 <0.05 & df3$p_1way_s1.s2 <0.05 & df3$p_1way_s1.s3 <0.05 & df3$p_1way_s2.s3 <0.05 &  
                             df3$p_2way_t.s1.s2 <0.05 & df3$p_2way_t.s1.s3 <0.05 & df3$p_2way_t.s2.s3 <0.05 & df3$p_2way_s1.s2.s3 <0.05, "fitting", "non-fitting")
  df3$positive04 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                             df3$p_3way >=0.1 & df3$p_1way_t.s1 <0.1 & df3$p_1way_t.s2 <0.1 & df3$p_1way_t.s3 <0.1 & df3$p_1way_s1.s2 <0.1 & df3$p_1way_s1.s3 <0.1 & df3$p_1way_s2.s3 <0.1 &  
                             df3$p_2way_t.s1.s2 <0.1 & df3$p_2way_t.s1.s3 <0.1 & df3$p_2way_t.s2.s3 <0.1 & df3$p_2way_s1.s2.s3 <0.1, "fitting", "non-fitting")
  df3$positive05 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                             df3$p_3way >=0.5 & df3$p_1way_t.s1 <0.5 & df3$p_1way_t.s2 <0.5 & df3$p_1way_t.s3 <0.5 & df3$p_1way_s1.s2 <0.5 & df3$p_1way_s1.s3 <0.5 & df3$p_1way_s2.s3 <0.5 &  
                             df3$p_2way_t.s1.s2 <0.5 & df3$p_2way_t.s1.s3 <0.5 & df3$p_2way_t.s2.s3 <0.5 & df3$p_2way_s1.s2.s3 <0.5, "fitting", "non-fitting")
  df3$positive06 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                             df3$p_3way >=0.05 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 &  
                             df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  df3$positive07 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 &
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & 
                             df3$p_3way >=0.5 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 &  
                             df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  
  df3$positive08 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.001 
                           & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  df3$positive09 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.01 
                           & df3$p_1way_t.s1 <0.01 & df3$p_1way_t.s2 <0.01 & df3$p_1way_t.s3 <0.01 & df3$p_1way_s1.s2 <0.01 & df3$p_1way_s1.s3 <0.01 & df3$p_1way_s2.s3 <0.01 
                           & df3$p_2way_t.s1.s2 <0.01 & df3$p_2way_t.s1.s3 <0.01 & df3$p_2way_t.s2.s3 <0.01 & df3$p_2way_s1.s2.s3 <0.01, "fitting", "non-fitting")
  df3$positive10 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.05 
                           & df3$p_1way_t.s1 <0.05 & df3$p_1way_t.s2 <0.05 & df3$p_1way_t.s3 <0.05 & df3$p_1way_s1.s2 <0.05 & df3$p_1way_s1.s3 <0.05 & df3$p_1way_s2.s3 <0.05 
                           & df3$p_2way_t.s1.s2 <0.05 & df3$p_2way_t.s1.s3 <0.05 & df3$p_2way_t.s2.s3 <0.05 & df3$p_2way_s1.s2.s3 <0.05, "fitting", "non-fitting")
  df3$positive11 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.1 
                           & df3$p_1way_t.s1 <0.1 & df3$p_1way_t.s2 <0.1 & df3$p_1way_t.s3 <0.1 & df3$p_1way_s1.s2 <0.1 & df3$p_1way_s1.s3 <0.1 & df3$p_1way_s2.s3 <0.1 
                           & df3$p_2way_t.s1.s2 <0.1 & df3$p_2way_t.s1.s3 <0.1 & df3$p_2way_t.s2.s3 <0.1 & df3$p_2way_s1.s2.s3 <0.1, "fitting", "non-fitting")
  df3$positive12 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.5 
                           & df3$p_1way_t.s1 <0.5 & df3$p_1way_t.s2 <0.5 & df3$p_1way_t.s3 <0.5 & df3$p_1way_s1.s2 <0.5 & df3$p_1way_s1.s3 <0.5 & df3$p_1way_s2.s3 <0.5 
                           & df3$p_2way_t.s1.s2 <0.5 & df3$p_2way_t.s1.s3 <0.5 & df3$p_2way_t.s2.s3 <0.5 & df3$p_2way_s1.s2.s3 <0.5, "fitting", "non-fitting")
  df3$positive13 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.05 
                           & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  df3$positive14 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 &df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 &  df3$p_3way >=0.5 
                           & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  
  df3$positive15 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.001 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  df3$positive16 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.01 & df3$p_1way_t.s1 <0.01 & df3$p_1way_t.s2 <0.01 & df3$p_1way_t.s3 <0.01 & df3$p_1way_s1.s2 <0.01 & df3$p_1way_s1.s3 <0.01 & df3$p_1way_s2.s3 <0.01 
                           & df3$p_2way_t.s1.s2 <0.01 & df3$p_2way_t.s1.s3 <0.01 & df3$p_2way_t.s2.s3 <0.01 & df3$p_2way_s1.s2.s3 <0.01, "fitting", "non-fitting")
  df3$positive17 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.05 & df3$p_1way_t.s1 <0.05 & df3$p_1way_t.s2 <0.05 & df3$p_1way_t.s3 <0.05 & df3$p_1way_s1.s2 <0.05 & df3$p_1way_s1.s3 <0.05 & df3$p_1way_s2.s3 <0.05 
                           & df3$p_2way_t.s1.s2 <0.05 & df3$p_2way_t.s1.s3 <0.05 & df3$p_2way_t.s2.s3 <0.05 & df3$p_2way_s1.s2.s3 <0.05, "fitting", "non-fitting")
  df3$positive18 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.1 & df3$p_1way_t.s1 <0.1 & df3$p_1way_t.s2 <0.1 & df3$p_1way_t.s3 <0.1 & df3$p_1way_s1.s2 <0.1 & df3$p_1way_s1.s3 <0.1 & df3$p_1way_s2.s3 <0.1 
                           & df3$p_2way_t.s1.s2 <0.1 & df3$p_2way_t.s1.s3 <0.1 & df3$p_2way_t.s2.s3 <0.1 & df3$p_2way_s1.s2.s3 <0.1, "fitting", "non-fitting")
  df3$positive19 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.5 & df3$p_1way_t.s1 <0.5 & df3$p_1way_t.s2 <0.5 & df3$p_1way_t.s3 <0.5 & df3$p_1way_s1.s2 <0.5 & df3$p_1way_s1.s3 <0.5 & df3$p_1way_s2.s3 <0.5 
                           & df3$p_2way_t.s1.s2 <0.5 & df3$p_2way_t.s1.s3 <0.5 & df3$p_2way_t.s2.s3 <0.5 & df3$p_2way_s1.s2.s3 <0.5, "fitting", "non-fitting")
  df3$positive20 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.05 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  df3$positive21 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.5 & df3$p_1way_t.s1 <0.001 & df3$p_1way_t.s2 <0.001 & df3$p_1way_t.s3 <0.001 & df3$p_1way_s1.s2 <0.001 & df3$p_1way_s1.s3 <0.001 & df3$p_1way_s2.s3 <0.001 
                           & df3$p_2way_t.s1.s2 <0.001 & df3$p_2way_t.s1.s3 <0.001 & df3$p_2way_t.s2.s3 <0.001 & df3$p_2way_s1.s2.s3 <0.001, "fitting", "non-fitting")
  
  df3$positive22 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.001, "fitting", "non-fitting")
  df3$positive23 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.01, "fitting", "non-fitting")
  df3$positive24 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.05, "fitting", "non-fitting")
  df3$positive25 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.1, "fitting", "non-fitting")
  df3$positive26 <- ifelse(df3$weight.s1-2*df3$se.s1 >0 & df3$weight.s1+2*df3$se.s1 >0 & df3$weight.s2-2*df3$se.s2 >0 & df3$weight.s2+2*df3$se.s2 >0 & df3$weight.s3-2*df3$se.s3 >0 & df3$weight.s3+2*df3$se.s3 >0 & 
                             df3$weight.s1-2*df3$se.s1 <1 & df3$weight.s1+2*df3$se.s1 <1 & df3$weight.s2-2*df3$se.s2 <1 & df3$weight.s2+2*df3$se.s2 <1 & df3$weight.s3-2*df3$se.s3 <1 & df3$weight.s3+2*df3$se.s3 <1 & df3$p_3way >=0.5, "fitting", "non-fitting")
  
  df3$positive27 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.001, "fitting", "non-fitting")
  df3$positive28 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.01, "fitting", "non-fitting")
  df3$positive29 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.05, "fitting", "non-fitting")
  df3$positive30 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.1, "fitting", "non-fitting")
  df3$positive31 <- ifelse(df3$weight.s1<1 & df3$weight.s2<1 & df3$weight.s3<1 & df3$weight.s1>0 & df3$weight.s2>0 & df3$weight.s3>0 & df3$p_3way >=0.5, "fitting", "non-fitting")
  
  df3$positive32 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.001, "fitting", "non-fitting")
  df3$positive33 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.01, "fitting", "non-fitting")
  df3$positive34 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.05, "fitting", "non-fitting")
  df3$positive35 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.1, "fitting", "non-fitting")
  df3$positive36 <- ifelse(df3$weight.s1<1.3 & df3$weight.s2<1.3 & df3$weight.s3<1.3 & df3$weight.s1>-0.3 & df3$weight.s2>-0.3 & df3$weight.s3>-0.3 & df3$p_3way >=0.5, "fitting", "non-fitting")
  
  write.table(df3, file=paste0(f,"_correct.txt"), row.names=F, sep='\t')
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

  ###add feasibility criteria###
  df4$positive01 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.001 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & 
                           df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & 
                           df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive02 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.01 & df4$p_1way_t.s1<0.01 & df4$p_1way_t.s2<0.01 & df4$p_1way_t.s3<0.01 & df4$p_1way_t.s4<0.01 & df4$p_1way_s1.s2<0.01 & df4$p_1way_s1.s3<0.01 & df4$p_1way_s1.s4<0.01 & df4$p_1way_s2.s3<0.01 & df4$p_1way_s2.s4<0.01 & 
                           df4$p_1way_s3.s4<0.01 & df4$p_2way_t.s1.s2<0.01 & df4$p_2way_t.s1.s3<0.01 & df4$p_2way_t.s1.s4<0.01 & df4$p_2way_t.s2.s3<0.01 & df4$p_2way_t.s2.s4<0.01 & df4$p_2way_t.s3.s4<0.01 & df4$p_2way_s1.s2.s3<0.01 & df4$p_2way_s1.s2.s4<0.01 & 
                           df4$p_2way_s1.s3.s4<0.01 & df4$p_2way_s2.s3.s4<0.01 & df4$p_3way_t.s1.s2.s3<0.01 & df4$p_3way_t.s1.s2.s4<0.01 & df4$p_3way_t.s1.s3.s4<0.01 & df4$p_3way_t.s2.s3.s4<0.01 & df4$p_3way_s1.s2.s3.s4<0.01, "fitting", "non-fitting")
  df4$positive03 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.05 & df4$p_1way_t.s1<0.05 & df4$p_1way_t.s2<0.05 & df4$p_1way_t.s3<0.05 & df4$p_1way_t.s4<0.05 & df4$p_1way_s1.s2<0.05 & df4$p_1way_s1.s3<0.05 & df4$p_1way_s1.s4<0.05 & df4$p_1way_s2.s3<0.05 & df4$p_1way_s2.s4<0.05 & 
                           df4$p_1way_s3.s4<0.05 & df4$p_2way_t.s1.s2<0.05 & df4$p_2way_t.s1.s3<0.05 & df4$p_2way_t.s1.s4<0.05 & df4$p_2way_t.s2.s3<0.05 & df4$p_2way_t.s2.s4<0.05 & df4$p_2way_t.s3.s4<0.05 & df4$p_2way_s1.s2.s3<0.05 & df4$p_2way_s1.s2.s4<0.05 & 
                           df4$p_2way_s1.s3.s4<0.05 & df4$p_2way_s2.s3.s4<0.05 & df4$p_3way_t.s1.s2.s3<0.05 & df4$p_3way_t.s1.s2.s4<0.05 & df4$p_3way_t.s1.s3.s4<0.05 & df4$p_3way_t.s2.s3.s4<0.05 & df4$p_3way_s1.s2.s3.s4<0.05, "fitting", "non-fitting")
  df4$positive04 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.1 & df4$p_1way_t.s1<0.1 & df4$p_1way_t.s2<0.1 & df4$p_1way_t.s3<0.1 & df4$p_1way_t.s4<0.1 & df4$p_1way_s1.s2<0.1 & df4$p_1way_s1.s3<0.1 & df4$p_1way_s1.s4<0.1 & df4$p_1way_s2.s3<0.1 & df4$p_1way_s2.s4<0.1 & 
                           df4$p_1way_s3.s4<0.1 & df4$p_2way_t.s1.s2<0.1 & df4$p_2way_t.s1.s3<0.1 & df4$p_2way_t.s1.s4<0.1 & df4$p_2way_t.s2.s3<0.1 & df4$p_2way_t.s2.s4<0.1 & df4$p_2way_t.s3.s4<0.1 & df4$p_2way_s1.s2.s3<0.1 & df4$p_2way_s1.s2.s4<0.1 & 
                           df4$p_2way_s1.s3.s4<0.1 & df4$p_2way_s2.s3.s4<0.1 & df4$p_3way_t.s1.s2.s3<0.1 & df4$p_3way_t.s1.s2.s4<0.1 & df4$p_3way_t.s1.s3.s4<0.1 & df4$p_3way_t.s2.s3.s4<0.1 & df4$p_3way_s1.s2.s3.s4<0.1, "fitting", "non-fitting")
  df4$positive05 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.5 & df4$p_1way_t.s1<0.5 & df4$p_1way_t.s2<0.5 & df4$p_1way_t.s3<0.5 & df4$p_1way_t.s4<0.5 & df4$p_1way_s1.s2<0.5 & df4$p_1way_s1.s3<0.5 & df4$p_1way_s1.s4<0.5 & df4$p_1way_s2.s3<0.5 & df4$p_1way_s2.s4<0.5 & 
                           df4$p_1way_s3.s4<0.5 & df4$p_2way_t.s1.s2<0.5 & df4$p_2way_t.s1.s3<0.5 & df4$p_2way_t.s1.s4<0.5 & df4$p_2way_t.s2.s3<0.5 & df4$p_2way_t.s2.s4<0.5 & df4$p_2way_t.s3.s4<0.5 & df4$p_2way_s1.s2.s3<0.5 & df4$p_2way_s1.s2.s4<0.5 & 
                           df4$p_2way_s1.s3.s4<0.5 & df4$p_2way_s2.s3.s4<0.5 & df4$p_3way_t.s1.s2.s3<0.5 & df4$p_3way_t.s1.s2.s4<0.5 & df4$p_3way_t.s1.s3.s4<0.5 & df4$p_3way_t.s2.s3.s4<0.5 & df4$p_3way_s1.s2.s3.s4<0.5, "fitting", "non-fitting")
  df4$positive06 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.05 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & 
                           df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & 
                           df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive07 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.5 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & 
                           df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & 
                           df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  
  df4$positive08 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.001 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive09 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.01 & df4$p_1way_t.s1<0.01 & df4$p_1way_t.s2<0.01 & df4$p_1way_t.s3<0.01 & df4$p_1way_t.s4<0.01 & df4$p_1way_s1.s2<0.01 & df4$p_1way_s1.s3<0.01 & df4$p_1way_s1.s4<0.01 & df4$p_1way_s2.s3<0.01 & df4$p_1way_s2.s4<0.01 & df4$p_1way_s3.s4<0.01 & df4$p_2way_t.s1.s2<0.01 & df4$p_2way_t.s1.s3<0.01 & df4$p_2way_t.s1.s4<0.01 & df4$p_2way_t.s2.s3<0.01 & df4$p_2way_t.s2.s4<0.01 & df4$p_2way_t.s3.s4<0.01 & df4$p_2way_s1.s2.s3<0.01 & df4$p_2way_s1.s2.s4<0.01 & df4$p_2way_s1.s3.s4<0.01 & df4$p_2way_s2.s3.s4<0.01 & df4$p_3way_t.s1.s2.s3<0.01 & df4$p_3way_t.s1.s2.s4<0.01 & df4$p_3way_t.s1.s3.s4<0.01 & df4$p_3way_t.s2.s3.s4<0.01 & df4$p_3way_s1.s2.s3.s4<0.01, "fitting", "non-fitting")
  df4$positive10 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.05 & df4$p_1way_t.s1<0.05 & df4$p_1way_t.s2<0.05 & df4$p_1way_t.s3<0.05 & df4$p_1way_t.s4<0.05 & df4$p_1way_s1.s2<0.05 & df4$p_1way_s1.s3<0.05 & df4$p_1way_s1.s4<0.05 & df4$p_1way_s2.s3<0.05 & df4$p_1way_s2.s4<0.05 & df4$p_1way_s3.s4<0.05 & df4$p_2way_t.s1.s2<0.05 & df4$p_2way_t.s1.s3<0.05 & df4$p_2way_t.s1.s4<0.05 & df4$p_2way_t.s2.s3<0.05 & df4$p_2way_t.s2.s4<0.05 & df4$p_2way_t.s3.s4<0.05 & df4$p_2way_s1.s2.s3<0.05 & df4$p_2way_s1.s2.s4<0.05 & df4$p_2way_s1.s3.s4<0.05 & df4$p_2way_s2.s3.s4<0.05 & df4$p_3way_t.s1.s2.s3<0.05 & df4$p_3way_t.s1.s2.s4<0.05 & df4$p_3way_t.s1.s3.s4<0.05 & df4$p_3way_t.s2.s3.s4<0.05 & df4$p_3way_s1.s2.s3.s4<0.05, "fitting", "non-fitting")
  df4$positive11 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.1 & df4$p_1way_t.s1<0.1 & df4$p_1way_t.s2<0.1 & df4$p_1way_t.s3<0.1 & df4$p_1way_t.s4<0.1 & df4$p_1way_s1.s2<0.1 & df4$p_1way_s1.s3<0.1 & df4$p_1way_s1.s4<0.1 & df4$p_1way_s2.s3<0.1 & df4$p_1way_s2.s4<0.1 & df4$p_1way_s3.s4<0.1 & df4$p_2way_t.s1.s2<0.1 & df4$p_2way_t.s1.s3<0.1 & df4$p_2way_t.s1.s4<0.1 & df4$p_2way_t.s2.s3<0.1 & df4$p_2way_t.s2.s4<0.1 & df4$p_2way_t.s3.s4<0.1 & df4$p_2way_s1.s2.s3<0.1 & df4$p_2way_s1.s2.s4<0.1 & df4$p_2way_s1.s3.s4<0.1 & df4$p_2way_s2.s3.s4<0.1 & df4$p_3way_t.s1.s2.s3<0.1 & df4$p_3way_t.s1.s2.s4<0.1 & df4$p_3way_t.s1.s3.s4<0.1 & df4$p_3way_t.s2.s3.s4<0.1 & df4$p_3way_s1.s2.s3.s4<0.1, "fitting", "non-fitting")
  df4$positive12 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.5 & df4$p_1way_t.s1<0.5 & df4$p_1way_t.s2<0.5 & df4$p_1way_t.s3<0.5 & df4$p_1way_t.s4<0.5 & df4$p_1way_s1.s2<0.5 & df4$p_1way_s1.s3<0.5 & df4$p_1way_s1.s4<0.5 & df4$p_1way_s2.s3<0.5 & df4$p_1way_s2.s4<0.5 & df4$p_1way_s3.s4<0.5 & df4$p_2way_t.s1.s2<0.5 & df4$p_2way_t.s1.s3<0.5 & df4$p_2way_t.s1.s4<0.5 & df4$p_2way_t.s2.s3<0.5 & df4$p_2way_t.s2.s4<0.5 & df4$p_2way_t.s3.s4<0.5 & df4$p_2way_s1.s2.s3<0.5 & df4$p_2way_s1.s2.s4<0.5 & df4$p_2way_s1.s3.s4<0.5 & df4$p_2way_s2.s3.s4<0.5 & df4$p_3way_t.s1.s2.s3<0.5 & df4$p_3way_t.s1.s2.s4<0.5 & df4$p_3way_t.s1.s3.s4<0.5 & df4$p_3way_t.s2.s3.s4<0.5 & df4$p_3way_s1.s2.s3.s4<0.5, "fitting", "non-fitting")
  df4$positive13 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.05 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive14 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.5 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  
  df4$positive15 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.001 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive16 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.01 & df4$p_1way_t.s1<0.01 & df4$p_1way_t.s2<0.01 & df4$p_1way_t.s3<0.01 & df4$p_1way_t.s4<0.01 & df4$p_1way_s1.s2<0.01 & df4$p_1way_s1.s3<0.01 & df4$p_1way_s1.s4<0.01 & df4$p_1way_s2.s3<0.01 & df4$p_1way_s2.s4<0.01 & df4$p_1way_s3.s4<0.01 & df4$p_2way_t.s1.s2<0.01 & df4$p_2way_t.s1.s3<0.01 & df4$p_2way_t.s1.s4<0.01 & df4$p_2way_t.s2.s3<0.01 & df4$p_2way_t.s2.s4<0.01 & df4$p_2way_t.s3.s4<0.01 & df4$p_2way_s1.s2.s3<0.01 & df4$p_2way_s1.s2.s4<0.01 & df4$p_2way_s1.s3.s4<0.01 & df4$p_2way_s2.s3.s4<0.01 & df4$p_3way_t.s1.s2.s3<0.01 & df4$p_3way_t.s1.s2.s4<0.01 & df4$p_3way_t.s1.s3.s4<0.01 & df4$p_3way_t.s2.s3.s4<0.01 & df4$p_3way_s1.s2.s3.s4<0.01, "fitting", "non-fitting")
  df4$positive17 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.05 & df4$p_1way_t.s1<0.05 & df4$p_1way_t.s2<0.05 & df4$p_1way_t.s3<0.05 & df4$p_1way_t.s4<0.05 & df4$p_1way_s1.s2<0.05 & df4$p_1way_s1.s3<0.05 & df4$p_1way_s1.s4<0.05 & df4$p_1way_s2.s3<0.05 & df4$p_1way_s2.s4<0.05 & df4$p_1way_s3.s4<0.05 & df4$p_2way_t.s1.s2<0.05 & df4$p_2way_t.s1.s3<0.05 & df4$p_2way_t.s1.s4<0.05 & df4$p_2way_t.s2.s3<0.05 & df4$p_2way_t.s2.s4<0.05 & df4$p_2way_t.s3.s4<0.05 & df4$p_2way_s1.s2.s3<0.05 & df4$p_2way_s1.s2.s4<0.05 & df4$p_2way_s1.s3.s4<0.05 & df4$p_2way_s2.s3.s4<0.05 & df4$p_3way_t.s1.s2.s3<0.05 & df4$p_3way_t.s1.s2.s4<0.05 & df4$p_3way_t.s1.s3.s4<0.05 & df4$p_3way_t.s2.s3.s4<0.05 & df4$p_3way_s1.s2.s3.s4<0.05, "fitting", "non-fitting")
  df4$positive18 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.1 & df4$p_1way_t.s1<0.1 & df4$p_1way_t.s2<0.1 & df4$p_1way_t.s3<0.1 & df4$p_1way_t.s4<0.1 & df4$p_1way_s1.s2<0.1 & df4$p_1way_s1.s3<0.1 & df4$p_1way_s1.s4<0.1 & df4$p_1way_s2.s3<0.1 & df4$p_1way_s2.s4<0.1 & df4$p_1way_s3.s4<0.1 & df4$p_2way_t.s1.s2<0.1 & df4$p_2way_t.s1.s3<0.1 & df4$p_2way_t.s1.s4<0.1 & df4$p_2way_t.s2.s3<0.1 & df4$p_2way_t.s2.s4<0.1 & df4$p_2way_t.s3.s4<0.1 & df4$p_2way_s1.s2.s3<0.1 & df4$p_2way_s1.s2.s4<0.1 & df4$p_2way_s1.s3.s4<0.1 & df4$p_2way_s2.s3.s4<0.1 & df4$p_3way_t.s1.s2.s3<0.1 & df4$p_3way_t.s1.s2.s4<0.1 & df4$p_3way_t.s1.s3.s4<0.1 & df4$p_3way_t.s2.s3.s4<0.1 & df4$p_3way_s1.s2.s3.s4<0.1, "fitting", "non-fitting")
  df4$positive19 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.5 & df4$p_1way_t.s1<0.5 & df4$p_1way_t.s2<0.5 & df4$p_1way_t.s3<0.5 & df4$p_1way_t.s4<0.5 & df4$p_1way_s1.s2<0.5 & df4$p_1way_s1.s3<0.5 & df4$p_1way_s1.s4<0.5 & df4$p_1way_s2.s3<0.5 & df4$p_1way_s2.s4<0.5 & df4$p_1way_s3.s4<0.5 & df4$p_2way_t.s1.s2<0.5 & df4$p_2way_t.s1.s3<0.5 & df4$p_2way_t.s1.s4<0.5 & df4$p_2way_t.s2.s3<0.5 & df4$p_2way_t.s2.s4<0.5 & df4$p_2way_t.s3.s4<0.5 & df4$p_2way_s1.s2.s3<0.5 & df4$p_2way_s1.s2.s4<0.5 & df4$p_2way_s1.s3.s4<0.5 & df4$p_2way_s2.s3.s4<0.5 & df4$p_3way_t.s1.s2.s3<0.5 & df4$p_3way_t.s1.s2.s4<0.5 & df4$p_3way_t.s1.s3.s4<0.5 & df4$p_3way_t.s2.s3.s4<0.5 & df4$p_3way_s1.s2.s3.s4<0.5, "fitting", "non-fitting")
  df4$positive20 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.05 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive21 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.5 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  
  df4$positive22 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.001, "fitting", "non-fitting")
  df4$positive23 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.01, "fitting", "non-fitting")
  df4$positive24 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.05, "fitting", "non-fitting")
  df4$positive25 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.1, "fitting", "non-fitting")
  df4$positive26 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                           df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                           df4$p_4way >=0.5, "fitting", "non-fitting")
 
  df4$positive27 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.001, "fitting", "non-fitting")
  df4$positive28 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.01, "fitting", "non-fitting")
  df4$positive29 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.05, "fitting", "non-fitting")
  df4$positive30 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.1, "fitting", "non-fitting")
  df4$positive31 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.5, "fitting", "non-fitting")
  
  df4$positive32 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.001, "fitting", "non-fitting")
  df4$positive33 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.01, "fitting", "non-fitting")
  df4$positive34 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.05, "fitting", "non-fitting")
  df4$positive35 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.1, "fitting", "non-fitting")
  df4$positive36 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.5, "fitting", "non-fitting")
  
  write.table(df4, file=paste0(f,"_correct.txt"), row.names=F, sep='\t')
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
  
  ###add feasibility criteria###
  df4$positive01 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.001 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & 
                             df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & 
                             df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive02 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.01 & df4$p_1way_t.s1<0.01 & df4$p_1way_t.s2<0.01 & df4$p_1way_t.s3<0.01 & df4$p_1way_t.s4<0.01 & df4$p_1way_s1.s2<0.01 & df4$p_1way_s1.s3<0.01 & df4$p_1way_s1.s4<0.01 & df4$p_1way_s2.s3<0.01 & df4$p_1way_s2.s4<0.01 & 
                             df4$p_1way_s3.s4<0.01 & df4$p_2way_t.s1.s2<0.01 & df4$p_2way_t.s1.s3<0.01 & df4$p_2way_t.s1.s4<0.01 & df4$p_2way_t.s2.s3<0.01 & df4$p_2way_t.s2.s4<0.01 & df4$p_2way_t.s3.s4<0.01 & df4$p_2way_s1.s2.s3<0.01 & df4$p_2way_s1.s2.s4<0.01 & 
                             df4$p_2way_s1.s3.s4<0.01 & df4$p_2way_s2.s3.s4<0.01 & df4$p_3way_t.s1.s2.s3<0.01 & df4$p_3way_t.s1.s2.s4<0.01 & df4$p_3way_t.s1.s3.s4<0.01 & df4$p_3way_t.s2.s3.s4<0.01 & df4$p_3way_s1.s2.s3.s4<0.01, "fitting", "non-fitting")
  df4$positive03 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.05 & df4$p_1way_t.s1<0.05 & df4$p_1way_t.s2<0.05 & df4$p_1way_t.s3<0.05 & df4$p_1way_t.s4<0.05 & df4$p_1way_s1.s2<0.05 & df4$p_1way_s1.s3<0.05 & df4$p_1way_s1.s4<0.05 & df4$p_1way_s2.s3<0.05 & df4$p_1way_s2.s4<0.05 & 
                             df4$p_1way_s3.s4<0.05 & df4$p_2way_t.s1.s2<0.05 & df4$p_2way_t.s1.s3<0.05 & df4$p_2way_t.s1.s4<0.05 & df4$p_2way_t.s2.s3<0.05 & df4$p_2way_t.s2.s4<0.05 & df4$p_2way_t.s3.s4<0.05 & df4$p_2way_s1.s2.s3<0.05 & df4$p_2way_s1.s2.s4<0.05 & 
                             df4$p_2way_s1.s3.s4<0.05 & df4$p_2way_s2.s3.s4<0.05 & df4$p_3way_t.s1.s2.s3<0.05 & df4$p_3way_t.s1.s2.s4<0.05 & df4$p_3way_t.s1.s3.s4<0.05 & df4$p_3way_t.s2.s3.s4<0.05 & df4$p_3way_s1.s2.s3.s4<0.05, "fitting", "non-fitting")
  df4$positive04 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.1 & df4$p_1way_t.s1<0.1 & df4$p_1way_t.s2<0.1 & df4$p_1way_t.s3<0.1 & df4$p_1way_t.s4<0.1 & df4$p_1way_s1.s2<0.1 & df4$p_1way_s1.s3<0.1 & df4$p_1way_s1.s4<0.1 & df4$p_1way_s2.s3<0.1 & df4$p_1way_s2.s4<0.1 & 
                             df4$p_1way_s3.s4<0.1 & df4$p_2way_t.s1.s2<0.1 & df4$p_2way_t.s1.s3<0.1 & df4$p_2way_t.s1.s4<0.1 & df4$p_2way_t.s2.s3<0.1 & df4$p_2way_t.s2.s4<0.1 & df4$p_2way_t.s3.s4<0.1 & df4$p_2way_s1.s2.s3<0.1 & df4$p_2way_s1.s2.s4<0.1 & 
                             df4$p_2way_s1.s3.s4<0.1 & df4$p_2way_s2.s3.s4<0.1 & df4$p_3way_t.s1.s2.s3<0.1 & df4$p_3way_t.s1.s2.s4<0.1 & df4$p_3way_t.s1.s3.s4<0.1 & df4$p_3way_t.s2.s3.s4<0.1 & df4$p_3way_s1.s2.s3.s4<0.1, "fitting", "non-fitting")
  df4$positive05 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.5 & df4$p_1way_t.s1<0.5 & df4$p_1way_t.s2<0.5 & df4$p_1way_t.s3<0.5 & df4$p_1way_t.s4<0.5 & df4$p_1way_s1.s2<0.5 & df4$p_1way_s1.s3<0.5 & df4$p_1way_s1.s4<0.5 & df4$p_1way_s2.s3<0.5 & df4$p_1way_s2.s4<0.5 & 
                             df4$p_1way_s3.s4<0.5 & df4$p_2way_t.s1.s2<0.5 & df4$p_2way_t.s1.s3<0.5 & df4$p_2way_t.s1.s4<0.5 & df4$p_2way_t.s2.s3<0.5 & df4$p_2way_t.s2.s4<0.5 & df4$p_2way_t.s3.s4<0.5 & df4$p_2way_s1.s2.s3<0.5 & df4$p_2way_s1.s2.s4<0.5 & 
                             df4$p_2way_s1.s3.s4<0.5 & df4$p_2way_s2.s3.s4<0.5 & df4$p_3way_t.s1.s2.s3<0.5 & df4$p_3way_t.s1.s2.s4<0.5 & df4$p_3way_t.s1.s3.s4<0.5 & df4$p_3way_t.s2.s3.s4<0.5 & df4$p_3way_s1.s2.s3.s4<0.5, "fitting", "non-fitting")
  df4$positive06 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.05 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & 
                             df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & 
                             df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive07 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.5 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & 
                             df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & 
                             df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  
  df4$positive08 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.001 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive09 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.01 & df4$p_1way_t.s1<0.01 & df4$p_1way_t.s2<0.01 & df4$p_1way_t.s3<0.01 & df4$p_1way_t.s4<0.01 & df4$p_1way_s1.s2<0.01 & df4$p_1way_s1.s3<0.01 & df4$p_1way_s1.s4<0.01 & df4$p_1way_s2.s3<0.01 & df4$p_1way_s2.s4<0.01 & df4$p_1way_s3.s4<0.01 & df4$p_2way_t.s1.s2<0.01 & df4$p_2way_t.s1.s3<0.01 & df4$p_2way_t.s1.s4<0.01 & df4$p_2way_t.s2.s3<0.01 & df4$p_2way_t.s2.s4<0.01 & df4$p_2way_t.s3.s4<0.01 & df4$p_2way_s1.s2.s3<0.01 & df4$p_2way_s1.s2.s4<0.01 & df4$p_2way_s1.s3.s4<0.01 & df4$p_2way_s2.s3.s4<0.01 & df4$p_3way_t.s1.s2.s3<0.01 & df4$p_3way_t.s1.s2.s4<0.01 & df4$p_3way_t.s1.s3.s4<0.01 & df4$p_3way_t.s2.s3.s4<0.01 & df4$p_3way_s1.s2.s3.s4<0.01, "fitting", "non-fitting")
  df4$positive10 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.05 & df4$p_1way_t.s1<0.05 & df4$p_1way_t.s2<0.05 & df4$p_1way_t.s3<0.05 & df4$p_1way_t.s4<0.05 & df4$p_1way_s1.s2<0.05 & df4$p_1way_s1.s3<0.05 & df4$p_1way_s1.s4<0.05 & df4$p_1way_s2.s3<0.05 & df4$p_1way_s2.s4<0.05 & df4$p_1way_s3.s4<0.05 & df4$p_2way_t.s1.s2<0.05 & df4$p_2way_t.s1.s3<0.05 & df4$p_2way_t.s1.s4<0.05 & df4$p_2way_t.s2.s3<0.05 & df4$p_2way_t.s2.s4<0.05 & df4$p_2way_t.s3.s4<0.05 & df4$p_2way_s1.s2.s3<0.05 & df4$p_2way_s1.s2.s4<0.05 & df4$p_2way_s1.s3.s4<0.05 & df4$p_2way_s2.s3.s4<0.05 & df4$p_3way_t.s1.s2.s3<0.05 & df4$p_3way_t.s1.s2.s4<0.05 & df4$p_3way_t.s1.s3.s4<0.05 & df4$p_3way_t.s2.s3.s4<0.05 & df4$p_3way_s1.s2.s3.s4<0.05, "fitting", "non-fitting")
  df4$positive11 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.1 & df4$p_1way_t.s1<0.1 & df4$p_1way_t.s2<0.1 & df4$p_1way_t.s3<0.1 & df4$p_1way_t.s4<0.1 & df4$p_1way_s1.s2<0.1 & df4$p_1way_s1.s3<0.1 & df4$p_1way_s1.s4<0.1 & df4$p_1way_s2.s3<0.1 & df4$p_1way_s2.s4<0.1 & df4$p_1way_s3.s4<0.1 & df4$p_2way_t.s1.s2<0.1 & df4$p_2way_t.s1.s3<0.1 & df4$p_2way_t.s1.s4<0.1 & df4$p_2way_t.s2.s3<0.1 & df4$p_2way_t.s2.s4<0.1 & df4$p_2way_t.s3.s4<0.1 & df4$p_2way_s1.s2.s3<0.1 & df4$p_2way_s1.s2.s4<0.1 & df4$p_2way_s1.s3.s4<0.1 & df4$p_2way_s2.s3.s4<0.1 & df4$p_3way_t.s1.s2.s3<0.1 & df4$p_3way_t.s1.s2.s4<0.1 & df4$p_3way_t.s1.s3.s4<0.1 & df4$p_3way_t.s2.s3.s4<0.1 & df4$p_3way_s1.s2.s3.s4<0.1, "fitting", "non-fitting")
  df4$positive12 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.5 & df4$p_1way_t.s1<0.5 & df4$p_1way_t.s2<0.5 & df4$p_1way_t.s3<0.5 & df4$p_1way_t.s4<0.5 & df4$p_1way_s1.s2<0.5 & df4$p_1way_s1.s3<0.5 & df4$p_1way_s1.s4<0.5 & df4$p_1way_s2.s3<0.5 & df4$p_1way_s2.s4<0.5 & df4$p_1way_s3.s4<0.5 & df4$p_2way_t.s1.s2<0.5 & df4$p_2way_t.s1.s3<0.5 & df4$p_2way_t.s1.s4<0.5 & df4$p_2way_t.s2.s3<0.5 & df4$p_2way_t.s2.s4<0.5 & df4$p_2way_t.s3.s4<0.5 & df4$p_2way_s1.s2.s3<0.5 & df4$p_2way_s1.s2.s4<0.5 & df4$p_2way_s1.s3.s4<0.5 & df4$p_2way_s2.s3.s4<0.5 & df4$p_3way_t.s1.s2.s3<0.5 & df4$p_3way_t.s1.s2.s4<0.5 & df4$p_3way_t.s1.s3.s4<0.5 & df4$p_3way_t.s2.s3.s4<0.5 & df4$p_3way_s1.s2.s3.s4<0.5, "fitting", "non-fitting")
  df4$positive13 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.05 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive14 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.5 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  
  df4$positive15 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.001 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive16 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.01 & df4$p_1way_t.s1<0.01 & df4$p_1way_t.s2<0.01 & df4$p_1way_t.s3<0.01 & df4$p_1way_t.s4<0.01 & df4$p_1way_s1.s2<0.01 & df4$p_1way_s1.s3<0.01 & df4$p_1way_s1.s4<0.01 & df4$p_1way_s2.s3<0.01 & df4$p_1way_s2.s4<0.01 & df4$p_1way_s3.s4<0.01 & df4$p_2way_t.s1.s2<0.01 & df4$p_2way_t.s1.s3<0.01 & df4$p_2way_t.s1.s4<0.01 & df4$p_2way_t.s2.s3<0.01 & df4$p_2way_t.s2.s4<0.01 & df4$p_2way_t.s3.s4<0.01 & df4$p_2way_s1.s2.s3<0.01 & df4$p_2way_s1.s2.s4<0.01 & df4$p_2way_s1.s3.s4<0.01 & df4$p_2way_s2.s3.s4<0.01 & df4$p_3way_t.s1.s2.s3<0.01 & df4$p_3way_t.s1.s2.s4<0.01 & df4$p_3way_t.s1.s3.s4<0.01 & df4$p_3way_t.s2.s3.s4<0.01 & df4$p_3way_s1.s2.s3.s4<0.01, "fitting", "non-fitting")
  df4$positive17 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.05 & df4$p_1way_t.s1<0.05 & df4$p_1way_t.s2<0.05 & df4$p_1way_t.s3<0.05 & df4$p_1way_t.s4<0.05 & df4$p_1way_s1.s2<0.05 & df4$p_1way_s1.s3<0.05 & df4$p_1way_s1.s4<0.05 & df4$p_1way_s2.s3<0.05 & df4$p_1way_s2.s4<0.05 & df4$p_1way_s3.s4<0.05 & df4$p_2way_t.s1.s2<0.05 & df4$p_2way_t.s1.s3<0.05 & df4$p_2way_t.s1.s4<0.05 & df4$p_2way_t.s2.s3<0.05 & df4$p_2way_t.s2.s4<0.05 & df4$p_2way_t.s3.s4<0.05 & df4$p_2way_s1.s2.s3<0.05 & df4$p_2way_s1.s2.s4<0.05 & df4$p_2way_s1.s3.s4<0.05 & df4$p_2way_s2.s3.s4<0.05 & df4$p_3way_t.s1.s2.s3<0.05 & df4$p_3way_t.s1.s2.s4<0.05 & df4$p_3way_t.s1.s3.s4<0.05 & df4$p_3way_t.s2.s3.s4<0.05 & df4$p_3way_s1.s2.s3.s4<0.05, "fitting", "non-fitting")
  df4$positive18 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.1 & df4$p_1way_t.s1<0.1 & df4$p_1way_t.s2<0.1 & df4$p_1way_t.s3<0.1 & df4$p_1way_t.s4<0.1 & df4$p_1way_s1.s2<0.1 & df4$p_1way_s1.s3<0.1 & df4$p_1way_s1.s4<0.1 & df4$p_1way_s2.s3<0.1 & df4$p_1way_s2.s4<0.1 & df4$p_1way_s3.s4<0.1 & df4$p_2way_t.s1.s2<0.1 & df4$p_2way_t.s1.s3<0.1 & df4$p_2way_t.s1.s4<0.1 & df4$p_2way_t.s2.s3<0.1 & df4$p_2way_t.s2.s4<0.1 & df4$p_2way_t.s3.s4<0.1 & df4$p_2way_s1.s2.s3<0.1 & df4$p_2way_s1.s2.s4<0.1 & df4$p_2way_s1.s3.s4<0.1 & df4$p_2way_s2.s3.s4<0.1 & df4$p_3way_t.s1.s2.s3<0.1 & df4$p_3way_t.s1.s2.s4<0.1 & df4$p_3way_t.s1.s3.s4<0.1 & df4$p_3way_t.s2.s3.s4<0.1 & df4$p_3way_s1.s2.s3.s4<0.1, "fitting", "non-fitting")
  df4$positive19 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.5 & df4$p_1way_t.s1<0.5 & df4$p_1way_t.s2<0.5 & df4$p_1way_t.s3<0.5 & df4$p_1way_t.s4<0.5 & df4$p_1way_s1.s2<0.5 & df4$p_1way_s1.s3<0.5 & df4$p_1way_s1.s4<0.5 & df4$p_1way_s2.s3<0.5 & df4$p_1way_s2.s4<0.5 & df4$p_1way_s3.s4<0.5 & df4$p_2way_t.s1.s2<0.5 & df4$p_2way_t.s1.s3<0.5 & df4$p_2way_t.s1.s4<0.5 & df4$p_2way_t.s2.s3<0.5 & df4$p_2way_t.s2.s4<0.5 & df4$p_2way_t.s3.s4<0.5 & df4$p_2way_s1.s2.s3<0.5 & df4$p_2way_s1.s2.s4<0.5 & df4$p_2way_s1.s3.s4<0.5 & df4$p_2way_s2.s3.s4<0.5 & df4$p_3way_t.s1.s2.s3<0.5 & df4$p_3way_t.s1.s2.s4<0.5 & df4$p_3way_t.s1.s3.s4<0.5 & df4$p_3way_t.s2.s3.s4<0.5 & df4$p_3way_s1.s2.s3.s4<0.5, "fitting", "non-fitting")
  df4$positive20 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.05 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  df4$positive21 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.5 & df4$p_1way_t.s1<0.001 & df4$p_1way_t.s2<0.001 & df4$p_1way_t.s3<0.001 & df4$p_1way_t.s4<0.001 & df4$p_1way_s1.s2<0.001 & df4$p_1way_s1.s3<0.001 & df4$p_1way_s1.s4<0.001 & df4$p_1way_s2.s3<0.001 & df4$p_1way_s2.s4<0.001 & df4$p_1way_s3.s4<0.001 & df4$p_2way_t.s1.s2<0.001 & df4$p_2way_t.s1.s3<0.001 & df4$p_2way_t.s1.s4<0.001 & df4$p_2way_t.s2.s3<0.001 & df4$p_2way_t.s2.s4<0.001 & df4$p_2way_t.s3.s4<0.001 & df4$p_2way_s1.s2.s3<0.001 & df4$p_2way_s1.s2.s4<0.001 & df4$p_2way_s1.s3.s4<0.001 & df4$p_2way_s2.s3.s4<0.001 & df4$p_3way_t.s1.s2.s3<0.001 & df4$p_3way_t.s1.s2.s4<0.001 & df4$p_3way_t.s1.s3.s4<0.001 & df4$p_3way_t.s2.s3.s4<0.001 & df4$p_3way_s1.s2.s3.s4<0.001, "fitting", "non-fitting")
  
  df4$positive22 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.001, "fitting", "non-fitting")
  df4$positive23 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.01, "fitting", "non-fitting")
  df4$positive24 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.05, "fitting", "non-fitting")
  df4$positive25 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.1, "fitting", "non-fitting")
  df4$positive26 <- ifelse(df4$weight.s1-2*df4$se.s1 >0 & df4$weight.s1+2*df4$se.s1 >0 & df4$weight.s2-2*df4$se.s2 >0 & df4$weight.s2+2*df4$se.s2 >0 & df4$weight.s3-2*df4$se.s3 >0 & df4$weight.s3+2*df4$se.s3 >0 & df4$weight.s4-2*df4$se.s4 >0 & df4$weight.s4+2*df4$se.s4 >0 & 
                             df4$weight.s1-2*df4$se.s1 <1 & df4$weight.s1+2*df4$se.s1 <1 & df4$weight.s2-2*df4$se.s2 <1 & df4$weight.s2+2*df4$se.s2 <1 & df4$weight.s3-2*df4$se.s3 <1 & df4$weight.s3+2*df4$se.s3 <1 & df4$weight.s4-2*df4$se.s4 <1 & df4$weight.s4+2*df4$se.s4 <1 & 
                             df4$p_4way >=0.5, "fitting", "non-fitting")
  
  df4$positive27 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.001, "fitting", "non-fitting")
  df4$positive28 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.01, "fitting", "non-fitting")
  df4$positive29 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.05, "fitting", "non-fitting")
  df4$positive30 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.1, "fitting", "non-fitting")
  df4$positive31 <- ifelse(df4$weight.s1<1 & df4$weight.s2<1 & df4$weight.s3<1 & df4$weight.s4<1 & df4$weight.s1>0 & df4$weight.s2>0 & df4$weight.s3>0 & df4$weight.s4>0 & df4$p_4way >=0.5, "fitting", "non-fitting")
  
  df4$positive32 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.001, "fitting", "non-fitting")
  df4$positive33 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.01, "fitting", "non-fitting")
  df4$positive34 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.05, "fitting", "non-fitting")
  df4$positive35 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.1, "fitting", "non-fitting")
  df4$positive36 <- ifelse(df4$weight.s1<1.3 & df4$weight.s2<1.3 & df4$weight.s3<1.3 & df4$weight.s4<1.3 & df4$weight.s1>-0.3 & df4$weight.s2>-0.3 & df4$weight.s3>-0.3 & df4$weight.s4>-0.3 & df4$p_4way >=0.5, "fitting", "non-fitting")
  
  write.table(df4, file=paste0(f,"_correct.txt"), row.names=F, sep='\t')
}


