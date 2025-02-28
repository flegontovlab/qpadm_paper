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
#df2 <- read_delim(paste0(f,"_2way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:20,63:65)]
#df3 <- read_delim(paste0(f,"_3way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:15,26,30,31,76:78)]
#df4 <- read_delim(paste0(f,"_4way_p3_correct.txt"), delim="\t", trim_ws=T)[,c(1:18,44,49,50,97:99)]

for (f in files) {
df2 <- read_delim(paste0(f,"_2way_p2_correct.txt"), delim="\t", trim_ws=T)[,c(1:20,63:65)]
df3 <- read_delim(paste0(f,"_3way_p2_correct.txt"), delim="\t", trim_ws=T)[,c(1:15,26,30,31,76:78)]
df4 <- read_delim(paste0(f,"_4way_p2_correct.txt"), delim="\t", trim_ws=T)[,c(1:18,44,49,50,97:99)]

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
df2 <- df2 %>% relocate(distance_s1.s2, .after=distance_t.s2)
rm(combs)
rm(s1)
rm(s2)
rm(res1)
rm(models)
rm(outtab2)

###generating a 1-way table
df5a <- df2[,c(1:4,5,6,13,16,17)]
df5b <- df2[,c(1:4,5,7,14,16,18)]
df5c <- df2[,c(1:4,6,7,15,16,19)] #s1*s2 pairs
colnames(df5a) <- c("pop_no", "qpAdm setup", "gene flow intensity per generation", "simulation replicate", "s1", "s2", "p_1way", "right", "distance_t.s1")
colnames(df5b) <- c("pop_no", "qpAdm setup", "gene flow intensity per generation", "simulation replicate", "s1", "s2", "p_1way", "right", "distance_t.s1")
colnames(df5c) <- c("pop_no", "qpAdm setup", "gene flow intensity per generation", "simulation replicate", "s1", "s2", "p_1way", "right", "distance_t.s1") #s1*s2 pairs

dff <- rbind(df5a, df5b, df5c)
dfff <- dff[with(dff, order(`gene flow intensity per generation`, `simulation replicate`, right, s1, s2)), ]
df1 <- dfff[!duplicated(dfff[, c("gene flow intensity per generation", "simulation replicate", "right", "s1", "s2")]), ]
colnames(df1)[colnames(df1) == "distance_t.s1"] <- "max. distance"
df1$`avg. distance` <- df1$`max. distance`
df1$dist_to_corner_maxdist <- df1$`max. distance`
df1$dist_to_corner_avgdist <- df1$`max. distance`
rm(df5a)
rm(df5b)
rm(df5c)
rm(dff)
rm(dfff)

###renaming and rearranging columns
colnames(df1)[colnames(df1) == "p_1way"] <- "p"
colnames(df2)[colnames(df2) == "p_2way"] <- "p"
colnames(df3)[colnames(df3) == "p_3way"] <- "p"
colnames(df4)[colnames(df4) == "p_4way"] <- "p"

df2$model <- "2-way"
df2 <- df2 %>% relocate(model, .after=target)
df3$model <- "3-way"
df3 <- df3 %>% relocate(model, .after=target)
df4$model <- "4-way"
df4 <- df4 %>% relocate(model, .after=target)

df1_1 <- df1
df1_2 <- df1
colnames(df1_1)[colnames(df1_1) == "s1"] <- "target"
colnames(df1_1)[colnames(df1_1) == "s2"] <- "s1"
colnames(df1_2)[colnames(df1_2) == "s2"] <- "target"
df1_2 <- df1_2 %>% relocate(target, .before=s1)
df1 <- rbind(df1_1, df1_2)
df1$model <- "1-way"
df1 <- df1 %>% relocate(model, .after=target)
rm(df1_1)
rm(df1_2)

df2 <- df2[,c(1:13,17,21:25)]


###calculating difference between estimated admixture fractions (EAF) and equal fractions (EF):
combs <- c(1:length(df2$target))
outtab2 = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(df2$weight.s1[x]-0.5), abs(df2$weight.s2[x]-0.5)))
  max_se = df2$se.s1[x]
  c("max_deviation_fromEF" = max_deviation_fromEF, "max_se" = max_se)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab2) <- c("max_deviation_fromEF", "max_se")
df2 <- cbind(df2, outtab2)

combs <- c(1:length(df3$target))
outtab2 = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(df3$weight.s1[x]-1/3), abs(df3$weight.s2[x]-1/3), abs(df3$weight.s3[x]-1/3)))
  max_se = max(c(df3$se.s1[x], df3$se.s2[x], df3$se.s3[x]), na.rm=T)
  c("max_deviation_fromEF" = max_deviation_fromEF, "max_se" = max_se)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab2) <- c("max_deviation_fromEF", "max_se")
df3 <- cbind(df3, outtab2)

combs <- c(1:length(df4$target))
outtab2 = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(df4$weight.s1[x]-0.25), abs(df4$weight.s2[x]-0.25), abs(df4$weight.s3[x]-0.25), abs(df4$weight.s4[x]-0.25)))
  max_se = max(c(df4$se.s1[x], df4$se.s2[x], df4$se.s3[x], df4$se.s4[x]), na.rm=T)
  c("max_deviation_fromEF" = max_deviation_fromEF, "max_se" = max_se)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab2) <- c("max_deviation_fromEF", "max_se")
df4 <- cbind(df4, outtab2)

###merging all the tables
dfa <- rbind.fill(df1, df2, df3, df4)
rm(combs)
rm(outtab2)
rm(df1)
rm(df2)
rm(df3)
rm(df4)

dfa <- dfa[with(dfa, order(model, pop_no, `qpAdm setup`, `gene flow intensity per generation`, `simulation replicate`, right, target, s1, s2, s3, s4)), ]
columns <-c("simulation replicate", "p", "max. distance", "avg. distance", "dist_to_corner_maxdist", "dist_to_corner_avgdist", "weight.s1", "se.s1", "weight.s2", "se.s2", "weight.s3", "se.s3", "weight.s4", "se.s4", "max_deviation_fromEF", "max_se")
dfa[, columns] <- lapply(columns, function(x) as.numeric(as.character(dfa[[x]])))

write.table(dfa, file=paste0(f,"_for2D.txt"), row.names=F, sep='\t')
}




###PLOTTING FIGURES###
setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm/2D_plots_p-values_adm_proportions/')

df <- read_bulk(directory = "./Isolation_by_distance/stepping-stone_simulations/qpAdm/2D_plots_p-values_adm_proportions/", name_contains = "_for2D.txt", fun = readr::read_delim, delim="\t")[,c(1:14,19:22,25)]

columns <-c("simulation.replicate", "p", "max..distance", "avg..distance", "dist_to_corner_maxdist", "dist_to_corner_avgdist", "max_deviation_fromEF", "max_se")
df[, columns] <- lapply(columns, function(x) as.numeric(as.character(df[[x]])))
df <- filter(df, model!="1-way")
#df_dedup <- df[!duplicated(df[, c("pop_no", "qpAdm.setup", "gene.flow.intensity.per.generation", "simulation.replicate", "target", "model", "s1", "s2", "s3", "s4", "right")]), ]

df21 <- filter(df, max..distance==1 & avg..distance==1 & model=="2-way" & min_S1_T_S2_angle==180)
df31 <- filter(df, max..distance==1 & avg..distance==1 & model=="3-way" & min_S1_T_S2_angle>116)
df41 <- filter(df, max..distance==1 & avg..distance==1 & model=="4-way" & min_S1_T_S2_angle>53)
df0 <- rbind(df21, df31, df41)


###correlations###
dfcorr_p1 <- filter(df, p>=1e-55)
dfcorr_p2 <- filter(df, p>=1e-5)
dfcorr_adm <- filter(df, max_deviation_fromEF<=0.5)
dfcorr_p1_adm <- filter(df, p>=1e-55, max_deviation_fromEF<=0.5)
dfcorr_p2_adm <- filter(df, p>=1e-5, max_deviation_fromEF<=0.5)

test21 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$max..distance, alternative="two.sided", method = "spearman"))) 
test22 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test23 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test24 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test25 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test21$filtering <- "max_deviation_fromEF<=0.5"
test22$filtering <- "max_deviation_fromEF<=0.5"
test23$filtering <- "max_deviation_fromEF<=0.5"
test24$filtering <- "max_deviation_fromEF<=0.5"
test25$filtering <- "max_deviation_fromEF<=0.5"
test21$var1 <- "max_deviation_fromEF"
test22$var1 <- "max_deviation_fromEF"
test23$var1 <- "max_deviation_fromEF"
test24$var1 <- "max_deviation_fromEF"
test25$var1 <- "max_deviation_fromEF"
test21$var2 <- "max..distance"
test22$var2 <- "avg..distance"
test23$var2 <- "min_S1_T_S2_angle"
test24$var2 <- "dist_to_corner_maxdist"
test25$var2 <- "dist_to_corner_avgdist"

test26 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$max..distance, alternative="two.sided", method = "spearman"))) 
test27 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test28 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test29 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test30 <- dfcorr_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test26$filtering <- "max_deviation_fromEF<=0.5"
test27$filtering <- "max_deviation_fromEF<=0.5"
test28$filtering <- "max_deviation_fromEF<=0.5"
test29$filtering <- "max_deviation_fromEF<=0.5"
test30$filtering <- "max_deviation_fromEF<=0.5"
test26$var1 <- "p"
test27$var1 <- "p"
test28$var1 <- "p"
test29$var1 <- "p"
test30$var1 <- "p"
test26$var2 <- "max..distance"
test27$var2 <- "avg..distance"
test28$var2 <- "min_S1_T_S2_angle"
test29$var2 <- "dist_to_corner_maxdist"
test30$var2 <- "dist_to_corner_avgdist"

test31 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$max..distance, alternative="two.sided", method = "spearman"))) 
test32 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test33 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test34 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test35 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test31$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test32$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test33$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test34$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test35$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test31$var1 <- "max_deviation_fromEF"
test32$var1 <- "max_deviation_fromEF"
test33$var1 <- "max_deviation_fromEF"
test34$var1 <- "max_deviation_fromEF"
test35$var1 <- "max_deviation_fromEF"
test31$var2 <- "max..distance"
test32$var2 <- "avg..distance"
test33$var2 <- "min_S1_T_S2_angle"
test34$var2 <- "dist_to_corner_maxdist"
test35$var2 <- "dist_to_corner_avgdist"

test36 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$max..distance, alternative="two.sided", method = "spearman"))) 
test37 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test38 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test39 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test40 <- dfcorr_p1_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test36$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test37$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test38$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test39$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test40$filtering <- "p>=1e-55, max_deviation_fromEF<=0.5"
test36$var1 <- "p"
test37$var1 <- "p"
test38$var1 <- "p"
test39$var1 <- "p"
test40$var1 <- "p"
test36$var2 <- "max..distance"
test37$var2 <- "avg..distance"
test38$var2 <- "min_S1_T_S2_angle"
test39$var2 <- "dist_to_corner_maxdist"
test40$var2 <- "dist_to_corner_avgdist"

test41 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$max..distance, alternative="two.sided", method = "spearman"))) 
test42 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test43 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test44 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test45 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test41$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test42$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test43$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test44$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test45$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test41$var1 <- "max_deviation_fromEF"
test42$var1 <- "max_deviation_fromEF"
test43$var1 <- "max_deviation_fromEF"
test44$var1 <- "max_deviation_fromEF"
test45$var1 <- "max_deviation_fromEF"
test41$var2 <- "max..distance"
test42$var2 <- "avg..distance"
test43$var2 <- "min_S1_T_S2_angle"
test44$var2 <- "dist_to_corner_maxdist"
test45$var2 <- "dist_to_corner_avgdist"

test46 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$max..distance, alternative="two.sided", method = "spearman"))) 
test47 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test48 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test49 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test50 <- dfcorr_p2_adm  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test46$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test47$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test48$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test49$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test50$filtering <- "p>=1e-5, max_deviation_fromEF<=0.5"
test46$var1 <- "p"
test47$var1 <- "p"
test48$var1 <- "p"
test49$var1 <- "p"
test50$var1 <- "p"
test46$var2 <- "max..distance"
test47$var2 <- "avg..distance"
test48$var2 <- "min_S1_T_S2_angle"
test49$var2 <- "dist_to_corner_maxdist"
test50$var2 <- "dist_to_corner_avgdist"

test01 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$max..distance, alternative="two.sided", method = "spearman"))) 
test02 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test03 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test04 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test05 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test01$filtering <- "p>=1e-55"
test02$filtering <- "p>=1e-55"
test03$filtering <- "p>=1e-55"
test04$filtering <- "p>=1e-55"
test05$filtering <- "p>=1e-55"
test01$var1 <- "max_deviation_fromEF"
test02$var1 <- "max_deviation_fromEF"
test03$var1 <- "max_deviation_fromEF"
test04$var1 <- "max_deviation_fromEF"
test05$var1 <- "max_deviation_fromEF"
test01$var2 <- "max..distance"
test02$var2 <- "avg..distance"
test03$var2 <- "min_S1_T_S2_angle"
test04$var2 <- "dist_to_corner_maxdist"
test05$var2 <- "dist_to_corner_avgdist"

test06 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$max..distance, alternative="two.sided", method = "spearman"))) 
test07 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test08 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test09 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test10 <- dfcorr_p1  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test06$filtering <- "p>=1e-55"
test07$filtering <- "p>=1e-55"
test08$filtering <- "p>=1e-55"
test09$filtering <- "p>=1e-55"
test10$filtering <- "p>=1e-55"
test06$var1 <- "p"
test07$var1 <- "p"
test08$var1 <- "p"
test09$var1 <- "p"
test10$var1 <- "p"
test06$var2 <- "max..distance"
test07$var2 <- "avg..distance"
test08$var2 <- "min_S1_T_S2_angle"
test09$var2 <- "dist_to_corner_maxdist"
test10$var2 <- "dist_to_corner_avgdist"

test11 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$max..distance, alternative="two.sided", method = "spearman"))) 
test12 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test13 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test14 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test15 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$max_deviation_fromEF, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test11$filtering <- "p>=1e-5"
test12$filtering <- "p>=1e-5"
test13$filtering <- "p>=1e-5"
test14$filtering <- "p>=1e-5"
test15$filtering <- "p>=1e-5"
test11$var1 <- "max_deviation_fromEF"
test12$var1 <- "max_deviation_fromEF"
test13$var1 <- "max_deviation_fromEF"
test14$var1 <- "max_deviation_fromEF"
test15$var1 <- "max_deviation_fromEF"
test11$var2 <- "max..distance"
test12$var2 <- "avg..distance"
test13$var2 <- "min_S1_T_S2_angle"
test14$var2 <- "dist_to_corner_maxdist"
test15$var2 <- "dist_to_corner_avgdist"

test16 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$max..distance, alternative="two.sided", method = "spearman"))) 
test17 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$avg..distance, alternative="two.sided", method = "spearman"))) 
test18 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$min_S1_T_S2_angle, alternative="two.sided", method = "spearman"))) 
test19 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_maxdist, alternative="two.sided", method = "spearman"))) 
test20 <- dfcorr_p2  %>% group_by(gene.flow.intensity.per.generation, model, qpAdm.setup, pop_no)  %>% do(tidy(cor.test(.$p, .$dist_to_corner_avgdist, alternative="two.sided", method = "spearman"))) 
test16$filtering <- "p>=1e-5"
test17$filtering <- "p>=1e-5"
test18$filtering <- "p>=1e-5"
test19$filtering <- "p>=1e-5"
test20$filtering <- "p>=1e-5"
test16$var1 <- "p"
test17$var1 <- "p"
test18$var1 <- "p"
test19$var1 <- "p"
test20$var1 <- "p"
test16$var2 <- "max..distance"
test17$var2 <- "avg..distance"
test18$var2 <- "min_S1_T_S2_angle"
test19$var2 <- "dist_to_corner_maxdist"
test20$var2 <- "dist_to_corner_avgdist"

corr<-rbind(test01,test02,test03,test04,test05,test06,test07,test08,test09,test10,test11,test12,test13,test14,test15,test16,test17,test18,test19,test20,test21,test22,test23,test24,test25,test26,test27,test28,test29,test30,test31,test32,test33,test34,test35,test36,test37,test38,test39,test40,test41,test42,test43,test44,test45,test46,test47,test48,test49,test50)
write.table(corr, file="Spearman_correlations_by_protocol.txt", sep='\t', row.names=F, quote=F)


###stratifying ideal symmetric models###
threshold05l <- df0 %>% group_by(gene.flow.intensity.per.generation, model) %>% filter(p<0.05) %>% summarise(count=n())
threshold01l <- df0 %>% group_by(gene.flow.intensity.per.generation, model) %>% filter(p<0.01) %>% summarise(count=n())
threshold05h <- df0 %>% group_by(gene.flow.intensity.per.generation, model) %>% filter(p>=0.05) %>% summarise(count=n())
threshold01h <- df0 %>% group_by(gene.flow.intensity.per.generation, model) %>% filter(p>=0.01) %>% summarise(count=n())

threshold2h <- df0 %>% group_by(gene.flow.intensity.per.generation, model) %>% filter(max_deviation_fromEF>1/2) %>% summarise(count=n())
threshold3h <- df0 %>% group_by(gene.flow.intensity.per.generation, model) %>% filter(max_deviation_fromEF>2/3) %>% summarise(count=n())
threshold4h <- df0 %>% group_by(gene.flow.intensity.per.generation, model) %>% filter(max_deviation_fromEF>3/4) %>% summarise(count=n())



###visualizing qpAdm models in the space formed by max.|estimated admixture fractions - equal fractions| and p-values###
###the results are stratified by model complexity and "sparsity" of stepping-stone landscapes
###("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2": these are approximate intervals for random sampling of gene-flow intensities from uniform or normal distributions)

###1. density of models populating the space###
p1density <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p)))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  geom_bin_2d(aes(fill=after_stat(count)), bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), breaks=c(-4, -2, 0, 2, 4), limits=c(-6,6))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(-50,-20,-10,-5), limits=c(-55,1))+
  scale_fill_gradientn(colours = pals::parula(20), trans="log", breaks = c(1,20,400,8100,162000), labels = c(1,20,400,8100,162000), na.value='grey90')+
  ggtitle("bins colored by: model count (logarithmic color scale)")+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,33), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p2density <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p)))+
    geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
    geom_bin_2d(aes(fill=after_stat(count)), bins = 35)+
    scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
    scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
    scale_fill_gradientn(colours = pals::parula(20), trans="log", breaks = c(1,20,400,8100,162000), labels = c(1,20,400,8100,162000), na.value='grey90')+
    geom_hline(yintercept=log10(0.001), linewidth=1, color="darkgreen", alpha=0.3)+
    geom_hline(yintercept=log10(0.5), linewidth=1, color="darkgreen", alpha=0.3)+
    ggtitle("bins colored by: model count (logarithmic color scale)")+
    #xlab("")+
    ylab("")+
    guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
    theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
          axis.text.x=element_text(size=28, face='bold'),
          axis.text.y=element_text(size=28, face='bold'),
          axis.title.x=element_blank(),
          axis.title.y=element_text(size=42, face='bold'),
          legend.position = "right",
          legend.key=element_rect(fill = 'white'),
          legend.title=element_text(size=32, family="mono", face = "bold"),
          legend.text=element_text(size=34),
          plot.margin = unit(c(0,0,0,45), "pt"),
          strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))



###2. the space colored by average source-target distance###
p1avgdist <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = avg..distance))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), breaks=c(-4, -2, 0, 2, 4), limits=c(-6,6))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(-50,-20,-10,-5), limits=c(-55,1))+
  scale_fill_gradientn(colours = pals::parula(20), breaks = c(0,1,2,3,4,5,6,7), labels = c(0,1,2,3,4,5,6,7), limits = c(0,7), na.value='grey90')+
  ggtitle('median average source-target distance')+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,33), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p2avgdist <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = avg..distance))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 35)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = pals::parula(20), breaks = c(0,1,2,3,4,5,6,7), labels = c(0,1,2,3,4,5,6,7), limits = c(0,7), na.value='grey90')+
  geom_hline(yintercept=log10(0.001), linewidth=1, color="darkgreen", alpha=0.3)+
  geom_hline(yintercept=log10(0.5), linewidth=1, color="darkgreen", alpha=0.3)+
  ggtitle('median average source-target distance')+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,45), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))


###3. the space colored by max. source-target distance###
p1maxdist <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = max..distance))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), breaks=c(-4, -2, 0, 2, 4), limits=c(-6,6))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(-50,-20,-10,-5), limits=c(-55,1))+
  scale_fill_gradientn(colours = pals::parula(20), breaks = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9), limits = c(0,9), na.value='grey90')+
  ggtitle('median max. source-target distance')+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,33), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p2maxdist <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = max..distance))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 35)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = pals::parula(20), breaks = c(0,1,2,3,4,5,6,7,8,9), labels = c(0,1,2,3,4,5,6,7,8,9), limits = c(0,9), na.value='grey90')+
  geom_hline(yintercept=log10(0.001), linewidth=1, color="darkgreen", alpha=0.3)+
  geom_hline(yintercept=log10(0.5), linewidth=1, color="darkgreen", alpha=0.3)+
  ggtitle('median max. source-target distance')+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,45), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))


###4. the space colored by min. source-target-source angle###
p1angle <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = min_S1_T_S2_angle))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), breaks=c(-4, -2, 0, 2, 4), limits=c(-6,6))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(-50,-20,-10,-5), limits=c(-55,1))+
  scale_fill_gradientn(colours = rev(pals::parula(20)), breaks = c(0,30,60,90,120,150,180), labels = c(0,30,60,90,120,150,180), limits = c(0,180), na.value='grey90')+
  ggtitle("median min. source-target-source angle")+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,33), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p2angle <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = min_S1_T_S2_angle))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 35)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = rev(pals::parula(20)), breaks = c(0,30,60,90,120,150,180), labels = c(0,30,60,90,120,150,180), limits = c(0,180), na.value='grey90')+
  geom_hline(yintercept=log10(0.001), linewidth=1, color="darkgreen", alpha=0.3)+
  geom_hline(yintercept=log10(0.5), linewidth=1, color="darkgreen", alpha=0.3)+
  ggtitle("median min. source-target-source angle")+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,45), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))


###5. the space colored by the distance to the ideal symmetric model based on average source-target distance###
p1dist_to_corner_avgdist <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = dist_to_corner_avgdist))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), breaks=c(-4, -2, 0, 2, 4), limits=c(-6,6))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(-50,-20,-10,-5), limits=c(-55,1))+
  scale_fill_gradientn(colours = pals::parula(20), breaks = c(0,2,4,6,8,10,12), labels = c(0,2,4,6,8,10,12), limits = c(0,12), na.value='grey90')+
  ggtitle('median Euclidean distance to the ideal symmetric model (in the space of avg. distance vs. min. angle)')+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,33), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p2dist_to_corner_avgdist <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = dist_to_corner_avgdist))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 35)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = pals::parula(20), breaks = c(0,2,4,6,8,10,12), labels = c(0,2,4,6,8,10,12), limits = c(0,12), na.value='grey90')+
  geom_hline(yintercept=log10(0.001), linewidth=1, color="darkgreen", alpha=0.3)+
  geom_hline(yintercept=log10(0.5), linewidth=1, color="darkgreen", alpha=0.3)+
  ggtitle('median Euclidean distance to the ideal symmetric model (in the space of avg. distance vs. min. angle)')+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,45), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))


###the space colored by the distance to the ideal symmetric model based on max. source-target distance###
p1dist_to_corner_maxdist <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = dist_to_corner_maxdist))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), breaks=c(-4, -2, 0, 2, 4), limits=c(-6,6))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(-50,-20,-10,-5), limits=c(-55,1))+
  scale_fill_gradientn(colours = pals::parula(20), breaks = c(0,2,4,6,8,10,12), labels = c(0,2,4,6,8,10,12), limits = c(0,12), na.value='grey90')+
  ggtitle('median Euclidean distance to the ideal symmetric model (in the space of max. distance vs. min. angle)')+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,33), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p2dist_to_corner_maxdist <- ggplot(df, aes(x = log10(max_deviation_fromEF), y = log10(p), z = dist_to_corner_maxdist))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="red", alpha=1)+
  stat_summary_2d(fun = median, bins = 35)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = pals::parula(20), breaks = c(0,2,4,6,8,10,12), labels = c(0,2,4,6,8,10,12), limits = c(0,12), na.value='grey90')+
  geom_hline(yintercept=log10(0.001), linewidth=1, color="darkgreen", alpha=0.3)+
  geom_hline(yintercept=log10(0.5), linewidth=1, color="darkgreen", alpha=0.3)+
  ggtitle('median Euclidean distance to the ideal symmetric model (in the space of max. distance vs. min. angle)')+
  #xlab("")+
  ylab("")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,45), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))


###6. visualizing ideal symmetric models
p1ideal <- ggplot(df0, aes(x = log10(max_deviation_fromEF), y = log10(p)))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="darkred", alpha=1)+
  geom_bin_2d(aes(fill=after_stat(count)), bins = 75)+
  scale_x_continuous(expand=c(0.01,0.01), breaks=c(-4, -2, 0, 2, 4), limits=c(-6,6))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(-50,-20,-10,-5), limits=c(-55,1))+
  scale_fill_gradientn(colours = pals::parula(20), na.value='grey90')+
  ggtitle('count of ideal symmetric models')+
  #xlab("")+
  ylab("log10(p-value)")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,33), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))

p2ideal <- ggplot(df0, aes(x = log10(max_deviation_fromEF), y = p))+
  geom_vline(xintercept=log10(0.5), linewidth=1, color="darkred", alpha=1)+
  geom_bin_2d(aes(fill=after_stat(count)), bins = 75)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(0,0.25,0.5,0.75,1))+
  scale_fill_gradientn(colours = pals::parula(20), na.value='grey90')+
  ggtitle('count of ideal symmetric models')+
  #xlab("")+
  ylab("p-value")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=18))+
  theme(plot.title = element_text(size=42, face='bold', color='darkmagenta'),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=42, face='bold'),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=34),
        plot.margin = unit(c(0,0,0,23), "pt"),
        strip.text=element_text(size=34, face='bold'))+
  facet_nested(. ~ model + factor(gene.flow.intensity.per.generation, levels=c("10-5 to 10-4", "10-4 to 10-3", "10-5 to 10-2", "10-3 to 10-2")))


###assembpling panels###
pall1 <- plot_grid(
    plot_grid(
      p1density+theme(legend.position="none"),
      p1avgdist+theme(legend.position="none"),
      p1angle+theme(legend.position="none"),
      p1dist_to_corner_avgdist+theme(legend.position="none"),
      p1ideal+theme(legend.position="none"),
      p2density+theme(legend.position="none"),
      p2avgdist+theme(legend.position="none"),
      p2angle+theme(legend.position="none"),
      p2dist_to_corner_avgdist+theme(legend.position="none"),
      p2ideal+theme(legend.position="none"),
      ncol = 1, align = "v"),
    plot_grid(
      get_legend(p1density),
      get_legend(p1avgdist),
      get_legend(p1angle),
      get_legend(p1dist_to_corner_avgdist),
      get_legend(p1ideal),
      get_legend(p2density),
      get_legend(p2avgdist),
      get_legend(p2angle),
      get_legend(p2dist_to_corner_avgdist),
      get_legend(p2ideal),
      ncol =1, align = "v"),
    rel_widths = c(24,2))

pall2 <- plot_grid(
  plot_grid(
    p1density+theme(legend.position="none"),
    p1maxdist+theme(legend.position="none"),
    p1angle+theme(legend.position="none"),
    p1dist_to_corner_maxdist+theme(legend.position="none"),
    p1ideal+theme(legend.position="none"),
    p2density+theme(legend.position="none"),
    p2maxdist+theme(legend.position="none"),
    p2angle+theme(legend.position="none"),
    p2dist_to_corner_maxdist+theme(legend.position="none"),
    p2ideal+theme(legend.position="none"),
    ncol = 1, align = "v"),
  plot_grid(
    get_legend(p1density),
    get_legend(p1maxdist),
    get_legend(p1angle),
    get_legend(p1dist_to_corner_maxdist),
    get_legend(p1ideal),
    get_legend(p2density),
    get_legend(p2maxdist),
    get_legend(p2angle),
    get_legend(p2dist_to_corner_maxdist),
    get_legend(p2ideal),
    ncol =1, align = "v"),
  rel_widths = c(24,2))
  

  ###this is Figure 5###
  pdf("2D_plot_p-values_adm_proportions_main_maxdist.pdf", width=46, height=56)
  print(annotate_figure(pall2, top=text_grob('qpAdm model metrics in the space of estimated admixture fractions and p-values: 3,200 independent randomized experiments per landscape type,\nall protocols combined, 24,568,800 2- to 4-way models in total', color="black", face="bold", size=44),
                               bottom=text_grob('log10( max. |estimated admixture fractions - equal fractions| )\n', color="black", face="bold", size=44)))
  dev.off()
  
  ###this is Figure S17###
  pdf("2D_plot_p-values_adm_proportions_main_avgdist.pdf", width=46, height=56)
  print(annotate_figure(pall1, top=text_grob('qpAdm model metrics in the space of estimated admixture fractions and p-values: 3,200 independent randomized experiments per landscape type,\nall protocols combined, 24,568,800 2- to 4-way models in total', color="black", face="bold", size=44),
                        bottom=text_grob('log10( max. |estimated admixture fractions - equal fractions| )\n', color="black", face="bold", size=44)))
  dev.off()
  
