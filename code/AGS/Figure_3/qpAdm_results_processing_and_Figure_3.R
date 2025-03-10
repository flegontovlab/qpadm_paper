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
library(ggsci)
library(pals)
library(rstatix)

future::plan(multisession)

setwd('./biases_project/qpAdm/simulations_Ulas/input_data_and_plots_for_revision')


###reading results for setups 1 and 2
df <- read_delim("qpadm_setups12_300-1000Mbp_800gen_rotatin_non-rotating.csv.gz", delim=",", trim_ws=T)[,c(1:16,19,20,23,24,27:31)]

###reading results for setup 3
dfa <- read_delim("qpadm_setup3_3000Mbp_3000gen_noisy_data_rotating.csv.gz", delim=",", trim_ws=T)
dfb <- read_delim("qpadm_setup3_3000Mbp_3000gen_rotating.csv.gz", delim=",", trim_ws=T)
dfc <- read_delim("qpadm_setup3_3000Mbp_3000gen_noisy_data_non-rotating.csv.gz", delim=",", trim_ws=T)
dfd <- read_delim("qpadm_setup3_3000Mbp_3000gen_non-rotating.csv.gz", delim=",", trim_ws=T)
dfa$max_depth_gen <- 3000
dfb$max_depth_gen <- 3000
dfc$max_depth_gen <- 3000
dfd$max_depth_gen <- 3000
dfa$genome_size <- "3,000 Mbp"
dfb$genome_size <- "3,000 Mbp"
dfc$genome_size <- "3,000 Mbp"
dfd$genome_size <- "3,000 Mbp"
dfa$data_quality <- "Noisy data"
dfb$data_quality <- "High-quality data"
dfc$data_quality <- "Noisy data"
dfd$data_quality <- "High-quality data"
dfa$qpAdm_protocol <- "rotating"
dfb$qpAdm_protocol <- "rotating"
dfc$qpAdm_protocol <- "non-rotating"
dfd$qpAdm_protocol <- "non-rotating"
df3 <- rbind(dfa, dfb, dfc, dfd)[,c(14:17,1:13)]
df3 <- filter(df3, target!="N")
df3 <- filter(df3, source2!="N")
rm(dfa)
rm(dfb)
rm(dfc)
rm(dfd)

###processing results for setup 3
df3$model_label <- paste0(df3$topology_no, "_", df3$target, "_", df3$source1, "_", df3$source2)
df3$topo_target <- paste0(df3$topology_no, "_", df3$target)
df3$topo_source1 <- paste0(df3$topology_no, "_", df3$source1)
df3$topo_source2 <- paste0(df3$topology_no, "_", df3$source2)
df3 <- df3 %>% relocate(model_label, .after=source2)
df3 <- df3 %>% relocate(topo_target, .after=model_label)
df3 <- df3 %>% relocate(topo_source1, .after=topo_target)
df3 <- df3 %>% relocate(topo_source2, .after=topo_source1)
df3$positive01 <- ifelse(df3$source1_weight-2*df3$source1_SE >0 & df3$source1_weight+2*df3$source1_SE >0 & df3$source2_weight-2*df3$source2_SE >0 & df3$source2_weight+2*df3$source2_SE >0 & df3$p_value >=0.01 
                         & df3$p_value_T_S1 <0.01 & df3$p_value_T_S2 <0.01 & df3$p_value_S1_S2 <0.01, "fitting", "non-fitting")

dates <- read.table("sampling_dates.txt", sep="\t", header=T)
df3 <- left_join(df3, dates[,c(1,4)], by="topo_target")
df3 <- left_join(df3, dates[,c(2,4)], by="topo_source1")
df3 <- left_join(df3, dates[,c(3,4)], by="topo_source2")
names(df3)[names(df3)=='sampling_time.x'] <- 'target_sampling_date'
names(df3)[names(df3)=='sampling_time.y'] <- 'source1_sampling_date'
names(df3)[names(df3)=='sampling_time'] <- 'source2_sampling_date'
df3 <- df3 %>% relocate(target_sampling_date, .after=topo_source2)
df3 <- df3 %>% relocate(source1_sampling_date, .after=target_sampling_date)
df3 <- df3 %>% relocate(source2_sampling_date, .after=source1_sampling_date)

###processing merged results###
df <- rbind(df, df3)
rm(df3)

df$temp_strat <- ifelse(df$source1_sampling_date<df$target_sampling_date | df$source2_sampling_date<df$target_sampling_date, "proximal", "distal")

write.table(filter(df, max_depth_gen==3000, positive01=="fitting"), "qpadm_setup3_30chr_fitting.txt", sep='\t', row.names=F)

combs <- c(1:length(df$target))
outtab2 = sapply(combs, function(x){
  max_deviation_from0.5 = max(c(abs(df$source1_weight[x]-0.5), abs(df$source2_weight[x]-0.5)))
  max_se = df$source1_SE[x]
  c("max_deviation_from0.5" = max_deviation_from0.5, "max_se" = max_se)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab2) <- c("max_deviation_from0.5", "max_se")
df <- cbind(df, outtab2)

#df$qpAdm_protocol <- gsub("rotating", "rot.", df$qpAdm_protocol)
df$max_depth_gen <- gsub(800, "800 gen.", df$max_depth_gen)
df$max_depth_gen <- gsub(3000, "3000 gen.", df$max_depth_gen)
df$genome_size <- gsub("300 Mbp", "300 Mbp, 800 gen.", df$genome_size)
df$genome_size <- gsub("1,000 Mbp", "1000 Mbp, 800 gen.", df$genome_size)
df$genome_size <- gsub("3,000 Mbp", "3000 Mbp, 3000 gen.", df$genome_size)

#dfd <- filter(df, temp_strat=="distal")
#dfdp <- df
#dfdp$temp_strat <- gsub("distal", "proximal", df$temp_strat)
#df <- rbind(dfd, dfdp)

conclusions <- read.table("model_classification_into_FP_TP.txt", sep="\t", header=T)
df <- left_join(df, conclusions, by="model_label")

###adding numadmix###
numadmix <- read.table("numadmix_summary.txt", sep="\t", header=T)
df <- left_join(df, numadmix[,c(1,4)], by="topo_target")
df <- left_join(df, numadmix[,c(2,4)], by="topo_source1")
df <- left_join(df, numadmix[,c(3,4)], by="topo_source2")
names(df)[names(df)=='nadmix.x'] <- 'target_nadmix'
names(df)[names(df)=='nadmix.y'] <- 'source1_nadmix'
names(df)[names(df)=='nadmix'] <- 'source2_nadmix'
df <- df %>% relocate(target_nadmix, .after=source2_sampling_date)
df <- df %>% relocate(source1_nadmix, .after=target_nadmix)
df <- df %>% relocate(source2_nadmix, .after=source1_nadmix)

combs <- c(1:length(df$target))
outtab2 = sapply(combs, function(x){
  max_source_nadmix = max(c(df$source1_nadmix[x], df$source2_nadmix[x]))
  avg_source_nadmix = (df$source1_nadmix[x] + df$source2_nadmix[x])/2
  c("max_source_nadmix" = max_source_nadmix, "avg_source_nadmix" = avg_source_nadmix)
}) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab2) <- c("max_source_nadmix", "avg_source_nadmix")
df <- cbind(df, outtab2)
df <- df %>% relocate(max_source_nadmix, .after=source2_nadmix)
df <- df %>% relocate(avg_source_nadmix, .after=max_source_nadmix)
df$target_avg_source_nadmix_diff <- df$target_nadmix - df$avg_source_nadmix
df <- df %>% relocate(target_avg_source_nadmix_diff, .after=avg_source_nadmix)


###calculating FDR###
FPc.dist <- df %>% group_by(qpAdm_protocol, data_quality, genome_size, simulation_replicate_no, .drop=F) %>% filter(conclusion=="FP", positive01=="fitting", temp_strat=="distal") %>% summarise(FP_count=n())
TPc.dist <- df %>% group_by(qpAdm_protocol, data_quality, genome_size, simulation_replicate_no, .drop=F) %>% filter(conclusion=="TP", positive01=="fitting", temp_strat=="distal") %>% summarise(TP_count=n())
FPc.dist$temp_strat <- "distal"
TPc.dist$temp_strat <- "distal"
FPc.dist <- FPc.dist %>% relocate(temp_strat, .before=qpAdm_protocol)
TPc.dist <- TPc.dist %>% relocate(temp_strat, .before=qpAdm_protocol)

FPc.prox <- df %>% group_by(qpAdm_protocol, data_quality, genome_size, simulation_replicate_no, .drop=F) %>% filter(conclusion=="FP", positive01=="fitting", temp_strat=="proximal") %>% summarise(FP_count=n())
TPc.prox <- df %>% group_by(qpAdm_protocol, data_quality, genome_size, simulation_replicate_no, .drop=F) %>% filter(conclusion=="TP", positive01=="fitting", temp_strat=="proximal") %>% summarise(TP_count=n())
FPc.prox$temp_strat <- "proximal"
TPc.prox$temp_strat <- "proximal"
FPc.prox <- FPc.prox %>% relocate(temp_strat, .before=qpAdm_protocol)
TPc.prox <- TPc.prox %>% relocate(temp_strat, .before=qpAdm_protocol)

FPc.all <- df %>% group_by(qpAdm_protocol, data_quality, genome_size, simulation_replicate_no, .drop=F) %>% filter(conclusion=="FP", positive01=="fitting") %>% summarise(FP_count=n())
TPc.all <- df %>% group_by(qpAdm_protocol, data_quality, genome_size, simulation_replicate_no, .drop=F) %>% filter(conclusion=="TP", positive01=="fitting") %>% summarise(TP_count=n())
FPc.all$temp_strat <- "dist. & prox."
TPc.all$temp_strat <- "dist. & prox."
FPc.all <- FPc.all %>% relocate(temp_strat, .before=qpAdm_protocol)
TPc.all <- TPc.all %>% relocate(temp_strat, .before=qpAdm_protocol)

FPc <- rbind(FPc.dist, FPc.prox, FPc.all)
TPc <- rbind(TPc.dist, TPc.prox, TPc.all)

FDR <- cbind(FPc, TPc$TP_count)
colnames(FDR) <- c(colnames(FPc), "TP_count")
FDR$FDR <- FDR$FP_count/(FDR$FP_count+FDR$TP_count)
#FDR$FPR <- FDR$FP_count/(858*40*0.823)
#FDR$TPR <- FDR$TP_count/(858*40*(1-0.823))
#min(filter(FDR, temp_strat=="dist. & prox.", qpAdm_protocol=="rotating")$FDR)
#max(filter(FDR, temp_strat=="dist. & prox.", qpAdm_protocol=="rotating")$FDR)
#min(filter(FDR, temp_strat=="dist. & prox.", qpAdm_protocol=="non-rotating")$FDR)
#max(filter(FDR, temp_strat=="dist. & prox.", qpAdm_protocol=="non-rotating")$FDR)
#med <- FDR %>% group_by(temp_strat, qpAdm_protocol, data_quality, genome_size, .drop=F) %>% summarise(median=median(FDR))
#min(filter(med, temp_strat=="proximal", qpAdm_protocol=="non-rotating")$median)
#mean(filter(med, temp_strat=="proximal", qpAdm_protocol=="non-rotating")$median)
#max(filter(med, temp_strat=="proximal", qpAdm_protocol=="non-rotating")$median)
FDR$genome_size <- gsub("300 Mbp, 800 gen.", "300M, 800g", FDR$genome_size)
FDR$genome_size <- gsub("1000 Mbp, 800 gen.", "1000M, 800g", FDR$genome_size)
FDR$genome_size <- gsub("3000 Mbp, 3000 gen.", "3000M, 3000g", FDR$genome_size)

###plotting box plots and violin plots###
stat.test2 <- FDR %>% group_by(qpAdm_protocol, temp_strat, data_quality) %>%  pairwise_wilcox_test(FDR ~ genome_size, paired=F) #non-paired tests should be used for comparisons across genome sizes since these are different simulations
stat.test2 <- stat.test2 %>% add_y_position()

p1 <- ggplot(FDR, aes(x=factor(genome_size, levels=c("300M, 800g","1000M, 800g","3000M, 3000g")), y=FDR))+
  geom_violin(aes(fill=factor(genome_size, levels=c("300M, 800g","1000M, 800g","3000M, 3000g"))), scale="width", trim=T, draw_quantiles=c(0.5))+
  geom_jitter(height = 0, width = 0.2, alpha=0.7)+
  scale_x_discrete(drop=F, na.translate=F)+
  #geom_text(data=filter(summary.data1, complexity_rotating!="not finished"), aes(x = closest_sources_number_rotating, y = -0.05, label=n), color="violetred", nudge_x=1, fontface=2, size=8)+
  scale_fill_manual(values=pal_frontiers("default")(3), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(0,0.25,0.5,0.75,1), limits=c(0,1.4))+
  ggtitle("influence of genome size, qpAdm protocols, and data quality on FDR")+
  xlab("genome size and simulation depth")+
  ylab('qpAdm FDR')+
  theme(plot.title=element_text(size=26, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=22),
        axis.title.x=element_text(size=26, face='bold'),
        axis.title.y=element_text(size=26, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=22, face = "bold"),
        legend.text=element_text(size=22),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.75,  linetype=3),
        panel.grid.major.x=element_blank(),
        strip.text=element_text(size=26, face='bold'),
        plot.margin = margin(5, 15, 20, 5, 'pt'))+
  facet_nested(data_quality ~ factor(qpAdm_protocol, levels=c("rotating","non-rotating")) + factor(temp_strat, levels=c("proximal","distal","dist. & prox.")))+
  stat_pvalue_manual(stat.test2, label = "p.adj.signif", hide.ns=F, tip.length = 0.01, size=7, bracket.size=0.5)

pdf("FDR.pdf", width = 16, height = 13)
print(p1)
dev.off()

FDR$qpAdm_protocol2 <- paste0(FDR$temp_strat,".", FDR$qpAdm_protocol)

###Wilcoxon tests###
stat.test1 <- FDR %>% group_by(qpAdm_protocol2, data_quality) %>%  pairwise_wilcox_test(FDR ~ genome_size, paired=F) #non-paired tests should be used for comparisons across genome sizes since these are different simulations
stat.test2 <- FDR %>% group_by(qpAdm_protocol2, genome_size) %>%  pairwise_wilcox_test(FDR ~ data_quality, paired=F) #non-paired tests should be used for comparisons across data qualities since these are different simulations
stat.test3 <- FDR %>% group_by(data_quality, genome_size) %>%  pairwise_wilcox_test(FDR ~ qpAdm_protocol2, paired=T) #paired tests should be used for comparisons across protocols
