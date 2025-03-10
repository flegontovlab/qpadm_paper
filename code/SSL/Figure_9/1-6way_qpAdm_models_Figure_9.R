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
library(pals)

future::plan(multisession)

setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm/systematic_sampling_Olas_results/6way_and_simpler_models/')

###reading all results###
df1 <- read_delim("results_1way.txt.gz", delim="\t", trim_ws=T)[,c(1:2,4:7)]
df2 <- read_delim("results_2way.txt.gz", delim="\t", trim_ws=T)[,c(1:2,4:8,10,12)]
df3 <- read_delim("results_3way.txt.gz", delim="\t", trim_ws=T)[,c(1:2,4:12)]
df4 <- read_delim("results_4way.txt.gz", delim="\t", trim_ws=T)[,c(1:2,4:14)]
df6 <- read_delim("results_6way.txt.gz", delim="\t", trim_ws=T)[,c(1:2,4:18)]

df1$model <- "1-way"
df2$model <- "2-way"
df3$model <- "3-way"
df4$model <- "4-way"
df6$model <- "6-way"
df1 <- df1 %>% relocate(model, .before=simulation_replicate)
df2 <- df2 %>% relocate(model, .before=simulation_replicate)
df3 <- df3 %>% relocate(model, .before=simulation_replicate)
df4 <- df4 %>% relocate(model, .before=simulation_replicate)
df6 <- df6 %>% relocate(model, .before=simulation_replicate)

###adding model feasibility criteria###
df1$positive02 <- ifelse(df1$p >=0.01, "fitting", "non-fitting")

df2$positive01 <- ifelse(df2$weight.s1<1 & df2$weight.s1>0 & df2$weight.s2<1 & df2$weight.s2>0 & df2$p >=0.01, "fitting", "non-fitting")
df2$positive02 <- ifelse(df2$p >=0.01, "fitting", "non-fitting")
df2$positive03 <- ifelse(df2$weight.s1<1 & df2$weight.s1>0 & df2$weight.s2<1 & df2$weight.s2>0, "fitting", "non-fitting")

df3$positive01 <- ifelse(df3$weight.s1<1 & df3$weight.s1>0 & df3$weight.s2<1 & df3$weight.s2>0 & df3$weight.s3<1 & df3$weight.s3>0 & df3$p >=0.01, "fitting", "non-fitting")
df3$positive02 <- ifelse(df3$p >=0.01, "fitting", "non-fitting")
df3$positive03 <- ifelse(df3$weight.s1<1 & df3$weight.s1>0 & df3$weight.s2<1 & df3$weight.s2>0 & df3$weight.s3<1 & df3$weight.s3>0, "fitting", "non-fitting")

df4$positive01 <- ifelse(df4$weight.s1<1 & df4$weight.s1>0 & df4$weight.s2<1 & df4$weight.s2>0 & df4$weight.s3<1 & df4$weight.s3>0 & df4$weight.s4<1 & df4$weight.s4>0 & df4$p >=0.01, "fitting", "non-fitting")
df4$positive02 <- ifelse(df4$p >=0.01, "fitting", "non-fitting")
df4$positive03 <- ifelse(df4$weight.s1<1 & df4$weight.s1>0 & df4$weight.s2<1 & df4$weight.s2>0 & df4$weight.s3<1 & df4$weight.s3>0 & df4$weight.s4<1 & df4$weight.s4>0, "fitting", "non-fitting")

df6$positive01 <- ifelse(df6$weight.s1<1 & df6$weight.s1>0 & df6$weight.s2<1 & df6$weight.s2>0 & df6$weight.s3<1 & df6$weight.s3>0 & df6$weight.s4<1 & df6$weight.s4>0 & df6$weight.s5<1 & df6$weight.s5>0 & df6$weight.s6<1 & df6$weight.s6>0 & df6$p >=0.01, "fitting", "non-fitting")
df6$positive02 <- ifelse(df6$p >=0.01, "fitting", "non-fitting")
df6$positive03 <- ifelse(df6$weight.s1<1 & df6$weight.s1>0 & df6$weight.s2<1 & df6$weight.s2>0 & df6$weight.s3<1 & df6$weight.s3>0 & df6$weight.s4<1 & df6$weight.s4>0 & df6$weight.s5<1 & df6$weight.s5>0 & df6$weight.s6<1 & df6$weight.s6>0, "fitting", "non-fitting")

###calculating max. difference between estimated admixture fractions and equal fractions (max. |EAF-EF|)###
df1$max_deviation_fromEF <- 1

combs <- c(1:length(df2$target))
outtab = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(df2$weight.s1[x]-1/2), abs(df2$weight.s2[x]-1/2)))
  c("max_deviation_fromEF" = max_deviation_fromEF)
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("max_deviation_fromEF")
df2 <- cbind(df2, outtab)

combs <- c(1:length(df3$target))
outtab = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(df3$weight.s1[x]-1/3), abs(df3$weight.s2[x]-1/3), abs(df3$weight.s3[x]-1/3)))
  c("max_deviation_fromEF" = max_deviation_fromEF)
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("max_deviation_fromEF")
df3 <- cbind(df3, outtab)

combs <- c(1:length(df4$target))
outtab = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(df4$weight.s1[x]-1/4), abs(df4$weight.s2[x]-1/4), abs(df4$weight.s3[x]-1/4), abs(df4$weight.s4[x]-1/4)))
  c("max_deviation_fromEF" = max_deviation_fromEF)
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("max_deviation_fromEF")
df4 <- cbind(df4, outtab)

combs <- c(1:length(df6$target))
outtab = sapply(combs, function(x){
  max_deviation_fromEF = max(c(abs(df6$weight.s1[x]-1/6), abs(df6$weight.s2[x]-1/6), abs(df6$weight.s3[x]-1/6), abs(df6$weight.s4[x]-1/6), abs(df6$weight.s5[x]-1/6), abs(df6$weight.s6[x]-1/6)))
  c("max_deviation_fromEF" = max_deviation_fromEF)
}) %>% as.data.frame() %>% set_rownames(NULL)
colnames(outtab) <- c("max_deviation_fromEF")
df6 <- cbind(df6, outtab)

###further processing###
df <- rbind.fill(df1, df2, df3, df4, df6)
df$log_p <- log10(df$p)
df$log_p[df$log_p==-Inf] <- -303

df$no_closest_sources <- gsub("0", "0 demes", df$no_closest_sources)
df$no_closest_sources <- gsub("1", "1 deme", df$no_closest_sources)
df$no_closest_sources <- gsub("2", "2 demes", df$no_closest_sources)
df$no_closest_sources <- gsub("3", "3 demes", df$no_closest_sources)
df$no_closest_sources <- gsub("4", "4 demes", df$no_closest_sources)
df$no_closest_sources <- gsub("5", "5 demes", df$no_closest_sources)
df$no_closest_sources <- gsub("6", "6 demes", df$no_closest_sources)
df$qpAdm_protocol[df$qpAdm_protocol == "non-rotating"] <- 'distal non-rotating'
df$qpAdm_protocol[df$qpAdm_protocol == "rotating"] <- 'distal rotating'

df1a <- df
df1a$max_deviation_fromEF[df1a$model!="1-way"] <- NA
df1a$log_p[df1a$model!="1-way"] <- NA
df2a <- df
df2a$max_deviation_fromEF[df2a$model!="2-way"] <- NA
df2a$log_p[df2a$model!="2-way"] <- NA
df3a <- df
df3a$max_deviation_fromEF[df3a$model!="3-way"] <- NA
df3a$log_p[df3a$model!="3-way"] <- NA
df4a <- df
df4a$max_deviation_fromEF[df4a$model!="4-way"] <- NA
df4a$log_p[df4a$model!="4-way"] <- NA
df6a <- df
df6a$max_deviation_fromEF[df6a$model!="6-way"] <- NA
df6a$log_p[df6a$model!="6-way"] <- NA

###plotting:
###max. |EAF-EF| vs. p-values###
###density of qpAdm models in this 2D space, stratified by sqpAdm protocol,
###by the number of target's nearest neighbors included in the model, and by model complexity###
p2density2 <- ggplot(df2a, aes(x = log10(max_deviation_fromEF), y = log_p))+
  geom_vline(xintercept=log10(1/2), linewidth=1, color="red", alpha=1)+
  geom_bin_2d(aes(fill=after_stat(count)), bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = pals::parula(20), trans="log", breaks = c(1,5,20,100), labels = c(1,5,20,100), limits = c(1,100), na.value='grey90')+
  geom_hline(yintercept=log10(0.01), linewidth=1, color="darkgreen", alpha=0.5)+
  ggtitle("21,600 2-way models")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=16))+
  theme(plot.title = element_text(size=44, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=38),
        plot.margin = unit(c(0,0,0,63), "pt"),
        strip.background = element_rect(fill = 'darkolivegreen2'), strip.text=element_text(size=42, face='bold'))+
  facet_nested(. ~ factor(qpAdm_protocol, levels=c("distal rotating", "distal non-rotating")) + no_closest_sources, drop=F)

p2density3 <- ggplot(df3a, aes(x = log10(max_deviation_fromEF), y = log_p))+
  geom_vline(xintercept=log10(2/3), linewidth=1, color="red", alpha=1)+
  geom_bin_2d(aes(fill=after_stat(count)), bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = pals::parula(20), trans="log", breaks = c(1,20,100,400), labels = c(1,20,100,400), limits = c(1,400), na.value='grey90')+
  geom_hline(yintercept=log10(0.01), linewidth=1, color="darkgreen", alpha=0.5)+
  ggtitle("112,000 3-way models")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=16))+
  theme(plot.title = element_text(size=44, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=38),
        plot.margin = unit(c(0,0,0,63), "pt"),
        strip.background = element_blank(), strip.text = element_blank())+
  facet_nested(. ~ factor(qpAdm_protocol, levels=c("distal rotating", "distal non-rotating")) + no_closest_sources, drop=F)

p2density4 <- ggplot(df4a, aes(x = log10(max_deviation_fromEF), y = log_p))+
  geom_vline(xintercept=log10(3/4), linewidth=1, color="red", alpha=1)+
  geom_bin_2d(aes(fill=after_stat(count)), bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = pals::parula(20), trans="log", breaks = c(1,20,100,400), labels = c(1,20,100,400), limits = c(1,400), na.value='grey90')+
  geom_hline(yintercept=log10(0.01), linewidth=1, color="darkgreen", alpha=0.5)+
  ggtitle("234,780 4-way models")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=16))+
  theme(plot.title = element_text(size=44, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=38),
        plot.margin = unit(c(0,0,0,63), "pt"),
        strip.background = element_blank(), strip.text = element_blank())+
  facet_nested(. ~ factor(qpAdm_protocol, levels=c("distal rotating", "distal non-rotating")) + no_closest_sources, drop=F)

p2density6 <- ggplot(df6a, aes(x = log10(max_deviation_fromEF), y = log_p))+
  geom_vline(xintercept=log10(5/6), linewidth=1, color="red", alpha=1)+
  geom_bin_2d(aes(fill=after_stat(count)), bins = 50)+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(-5,2.1))+
  scale_y_continuous(expand=c(0.01,0.01), breaks = c(-5,-4,-3,-2,-1,0), limits=c(-5.5,0.1))+
  scale_fill_gradientn(colours = pals::parula(20), trans="log", breaks = c(1,20,400,4000), labels = c(1,20,400,4000), limits = c(1,4000), na.value='grey90')+
  geom_hline(yintercept=log10(0.01), linewidth=1, color="darkgreen", alpha=0.5)+
  ggtitle("1,601,600 6-way models")+
  guides(fill = guide_colorbar(title = NULL, order=1, barwidth=4, barheight=16))+
  theme(plot.title = element_text(size=44, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=28, face='bold'),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "right",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=38),
        plot.margin = unit(c(0,0,0,63), "pt"),
        strip.background = element_blank(), strip.text = element_blank())+
  facet_nested(. ~ factor(qpAdm_protocol, levels=c("distal rotating", "distal non-rotating")) + no_closest_sources, drop=F)

###assembling panels###
p2 <- plot_grid(
  plot_grid(
    p2density2+theme(legend.position="none"),
    p2density3+theme(legend.position="none"),
    p2density4+theme(legend.position="none"),
    p2density6+theme(legend.position="none"),
    ncol = 1, align = "v", rel_heights = c(1.3,1,1,1)),
  plot_grid(
    get_legend(p2density2),
    get_legend(p2density3),
    get_legend(p2density4),
    get_legend(p2density6),
    ncol = 1, align = "v", rel_heights = c(1.3,1,1,1)),
  rel_widths = c(24,2))

p2a <- annotate_figure(p2, bottom=text_grob('log10( max. |estimated admixture fractions - equal fractions| )\n', color="black", face="bold", size=44),
                          left=text_grob('log10(p-value)', color="black", face="bold", size=44, rot=90))


###preparing data for barplots:
summary01f <- df %>% group_by(model, qpAdm_protocol, no_closest_sources, .drop=F) %>% filter(positive01=="fitting") %>% summarise(f=n())
summary01nf <- df %>% group_by(model, qpAdm_protocol, no_closest_sources, .drop=F) %>% filter(positive01=="non-fitting") %>% summarise(nf=n())
summary02f <- df %>% group_by(model, qpAdm_protocol, no_closest_sources, .drop=F) %>% filter(positive02=="fitting") %>% summarise(f=n())
summary02nf <- df %>% group_by(model, qpAdm_protocol, no_closest_sources, .drop=F) %>% filter(positive02=="non-fitting") %>% summarise(nf=n())
summary03f <- df %>% group_by(model, qpAdm_protocol, no_closest_sources, .drop=F) %>% filter(positive03=="fitting") %>% summarise(f=n())
summary03nf <- df %>% group_by(model, qpAdm_protocol, no_closest_sources, .drop=F) %>% filter(positive03=="non-fitting") %>% summarise(nf=n())
summary <- cbind(summary01f, summary01nf$nf, summary02f$f, summary02nf$nf, summary03f$f, summary03nf$nf)
colnames(summary) <- c("model", "qpAdm_protocol", "no_closest_sources", "f01","nf01","f02","nf02","f03","nf03")
summary$PR01 <- summary$f01/(summary$f01+summary$nf01)
summary$PR02 <- summary$f02/(summary$f02+summary$nf02)
summary$PR03 <- summary$f03/(summary$f03+summary$nf03)

summary1 <- summary[,c(1:3,4,5,10)]
summary2 <- summary[,c(1:3,6,7,11)]
summary3 <- summary[,c(1:3,8,9,12)]
colnames(summary1) <- c("model", "qpAdm_protocol", "no_closest_sources", "f", "nf", "PR")
colnames(summary2) <- c("model", "qpAdm_protocol", "no_closest_sources", "f", "nf", "PR")
colnames(summary3) <- c("model", "qpAdm_protocol", "no_closest_sources", "f", "nf", "PR")
summary1$criteria <- "p>=0.01, %(0,1)"
summary2$criteria <- "p>=0.01"
summary3$criteria <- "%(0,1)"
summary_final <- rbind(summary1, summary2, summary3)
summary_final <- filter(summary_final, model!="1-way")

###barplots showing qpAdm positive rate (percentage of models satisfying selected feasibility criteria) stratified by qpAdm protocol,
###by the number of target's nearest neighbors included in the model, and by model complexity
pbar <- ggplot(summary_final, aes(x=factor(criteria, levels=c("p>=0.01, %(0,1)", "p>=0.01", "%(0,1)")), y=PR))+ 
  geom_bar(aes(fill=factor(criteria, levels=c("p>=0.01, %(0,1)", "p>=0.01", "%(0,1)")), color=factor(criteria, levels=c("p>=0.01, %(0,1)", "p>=0.01", "%(0,1)"))), stat="identity", position="dodge", width=.8)+
  geom_text(data=summary_final, aes(x = factor(criteria, levels=c("p>=0.01, %(0,1)", "p>=0.01", "%(0,1)")), y = -0.07, label=round(PR, digits=2)), color="violetred", fontface=2, size=10)+
  scale_fill_manual(values=pal_frontiers("default")(3))+
  scale_color_manual(values=pal_frontiers("default")(3))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(0,0.25,0.5,0.75,1), limits=c(-0.12,1))+
  ggtitle("\nqpAdm positive rate for different feasibility criteria")+
  xlab('feasibility criteria')+
  ylab('positive rate')+
  theme(plot.title = element_text(size=44, face='bold', color="darkmagenta"),
        axis.text.x=element_text(size=34, face='bold', angle=45, hjust=0.95, vjust=0.95),
        axis.text.y=element_text(size=28, face='bold'),
        axis.title.x=element_text(size=44, face='bold'),
        axis.title.y=element_text(size=44, face='bold'),
        legend.position = "none",
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=32, family="mono", face = "bold"),
        legend.text=element_text(size=38),
        plot.margin = unit(c(0,210,0,35), "pt"),
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.75,  linetype=3),
        panel.grid.major.x=element_blank(),
        strip.background = element_rect(fill = 'darkolivegreen2'), strip.text=element_text(size=42, face='bold'))+
  facet_nested(model ~ factor(qpAdm_protocol, levels=c("distal rotating", "distal non-rotating")) + no_closest_sources, drop=F)

###merging all panels:
pall <- plot_grid(p2a, pbar, ncol=1, align = "v", rel_widths = c(1,1))

###This is Figure 9###
pdf("2D_plot_p-values_adm_proportions_1-6way_main_correct.pdf", width=47, height=41)
print(annotate_figure(pall, top=text_grob('6-way and simpler qpAdm models in ideal conditions: systematic experiments on the "10-5 to 10-2" landscapes with all 6 demes from the 1st circle available,\nand 10 demes from the 3rd circle (& 16 300-gen.-old "right" demes from the 3rd circle in the non-rotating setup);\n1,973,020 models stratified by the number of demes from the 1st circle', color="purple3", face="bold", size=44)))
dev.off()
