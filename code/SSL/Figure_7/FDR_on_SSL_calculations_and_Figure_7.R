library(readr)
#library(plyr)
library(dplyr)
library(egg)
library(tidyverse)
library(future)
library(future.apply)
library(magrittr)
future::plan(multisession)

setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm/FDR/')
fit_criteria <- read_delim(file="fit_criteria2.txt", delim="\t")

###reading 2-way and 3-way results###
df2 <- read_bulk(directory = "./Isolation_by_distance/stepping-stone_simulations/qpAdm/input_files_basic_1-4way/", name_contains = "_2way_", fun = readr::read_delim, delim="\t")[,c(1:4,21:56,64,65)]
df2$model <- "2-way"
df2 <- df2 %>% relocate(model, .before=pop_no)
df3 <- read_bulk(directory = "./Isolation_by_distance/stepping-stone_simulations/qpAdm/input_files_basic_1-4way/", name_contains = "_3way_", fun = readr::read_delim, delim="\t")[,c(1:4,32:67,77,78)]
df3$model <- "3-way"
df3 <- df3 %>% relocate(model, .before=pop_no)
df4 <- read_bulk(directory = "./Isolation_by_distance/stepping-stone_simulations/qpAdm/input_files_basic_1-4way/", name_contains = "_4way_", fun = readr::read_delim, delim="\t")[,c(1:4,51:86,98,99)]
df4$model <- "4-way"
df4 <- df4 %>% relocate(model, .before=pop_no)

df <- rbind(df2, df3, df4)
saveRDS(df, file="all_models.rds")
rm(df2)
rm(df3)
rm(df4)

dfa <- filter(df, !is.na(dist_to_corner_maxdist))
dfa$actual_maxdist <- ifelse(dfa$dist_to_corner_maxdist<=5, "optimal", "non-optimal")
dfa$actual_avgdist <- ifelse(dfa$dist_to_corner_avgdist<=4, "optimal", "non-optimal")

crs = c(32,28,22,15,9,3,7) #v. 2


###calculating FDR###
###per simulation replicate:
outtab <- data.frame(model=character(),
                      gene.flow.intensity.per.generation=character(),
                      simulation.replicate=integer(),
                      number=integer(),
                      fit_criteria1=character(),
                      fit_criteria2=character(),
                      fit_criteria3=character(),
                      odds_maxdist=double(),
                      odds_avgdist=double(),
                      FPR_maxdist=double(),
                      FPR_avgdist=double(),
                      FNR_maxdist=double(),
                      FNR_avgdist=double(),
                      FDR_maxdist=double(),
                      FDR_avgdist=double(),
                      stringsAsFactors=FALSE)

for (i in crs) {
  dfb <- dfa[,c(1:5,5+i,44,45)]
  colnames(dfb) <- c("model", "pop_no", "qpAdm.setup", "gene.flow.intensity.per.generation", "simulation.replicate", "positive", "actual_maxdist", "actual_avgdist")  
  FT_maxdist <- dfb %>% group_by(model, gene.flow.intensity.per.generation, simulation.replicate, actual_maxdist, .drop=F) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(model, gene.flow.intensity.per.generation, simulation.replicate, actual_maxdist) %>% mutate(n = replace_na(n, 0))
  FT_avgdist <- dfb %>% group_by(model, gene.flow.intensity.per.generation, simulation.replicate, actual_avgdist, .drop=F) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(model, gene.flow.intensity.per.generation, simulation.replicate, actual_avgdist) %>% mutate(n = replace_na(n, 0))
  confusion_maxdist <- dfb %>% group_by(model, gene.flow.intensity.per.generation, simulation.replicate, positive, actual_maxdist, .drop=F) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(model, gene.flow.intensity.per.generation, simulation.replicate, positive, actual_maxdist) %>% mutate(n = replace_na(n, 0))
  confusion_avgdist <- dfb %>% group_by(model, gene.flow.intensity.per.generation, simulation.replicate, positive, actual_avgdist, .drop=F) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(model, gene.flow.intensity.per.generation, simulation.replicate, positive, actual_avgdist) %>% mutate(n = replace_na(n, 0))
  results <- FT_maxdist[,c(1:3)]
  results <- results[!duplicated(results[, c("model", "gene.flow.intensity.per.generation", "simulation.replicate")]), ]
  results$number <- i
  results <- left_join(results, fit_criteria, by="number")
  results$odds_maxdist <- filter(FT_maxdist, actual_maxdist=="optimal")$n/filter(FT_maxdist, actual_maxdist=="non-optimal")$n
  results$odds_avgdist <- filter(FT_avgdist, actual_avgdist=="optimal")$n/filter(FT_avgdist, actual_avgdist=="non-optimal")$n
  results$FPR_maxdist <- filter(confusion_maxdist, positive=="fitting", actual_maxdist=="non-optimal")$n/(filter(confusion_maxdist, positive=="fitting", actual_maxdist=="non-optimal")$n + filter(confusion_maxdist, positive=="non-fitting", actual_maxdist=="non-optimal")$n)
  results$FPR_avgdist <- filter(confusion_avgdist, positive=="fitting", actual_avgdist=="non-optimal")$n/(filter(confusion_avgdist, positive=="fitting", actual_avgdist=="non-optimal")$n + filter(confusion_avgdist, positive=="non-fitting", actual_avgdist=="non-optimal")$n)
  results$FNR_maxdist <- filter(confusion_maxdist, positive=="non-fitting", actual_maxdist=="optimal")$n/(filter(confusion_maxdist, positive=="fitting", actual_maxdist=="optimal")$n + filter(confusion_maxdist, positive=="non-fitting", actual_maxdist=="optimal")$n)
  results$FNR_avgdist <- filter(confusion_avgdist, positive=="non-fitting", actual_avgdist=="optimal")$n/(filter(confusion_avgdist, positive=="fitting", actual_avgdist=="optimal")$n + filter(confusion_avgdist, positive=="non-fitting", actual_avgdist=="optimal")$n)
  results$FDR_maxdist <- filter(confusion_maxdist, positive=="fitting", actual_maxdist=="non-optimal")$n/(filter(confusion_maxdist, positive=="fitting", actual_maxdist=="optimal")$n + filter(confusion_maxdist, positive=="fitting", actual_maxdist=="non-optimal")$n)
  results$FDR_avgdist <- filter(confusion_avgdist, positive=="fitting", actual_avgdist=="non-optimal")$n/(filter(confusion_avgdist, positive=="fitting", actual_avgdist=="optimal")$n + filter(confusion_avgdist, positive=="fitting", actual_avgdist=="non-optimal")$n)
  outtab <- rbind(outtab, results)
}
write.table(outtab, file="FDR_per_replicate.txt", row.names=F, sep='\t')


###per qpAdm protocol and simulation replicate:
outtab <- data.frame(pop_no=character(),
                     model=character(),
                     gene.flow.intensity.per.generation=character(),
                     qpAdm.setup=character(),
                     simulation.replicate=integer(),
                     number=integer(),
                     fit_criteria1=character(),
                     fit_criteria2=character(),
                     fit_criteria3=character(),
                     odds_maxdist=double(),
                     odds_avgdist=double(),
                     FPR_maxdist=double(),
                     FPR_avgdist=double(),
                     FNR_maxdist=double(),
                     FNR_avgdist=double(),
                     FDR_maxdist=double(),
                     FDR_avgdist=double(),
                     stringsAsFactors=FALSE)

for (i in crs) {
  dfb <- dfa[,c(2,1,4,3,5,5+i,44,45)]
  colnames(dfb) <- c("pop_no", "model", "gene.flow.intensity.per.generation", "qpAdm.setup", "simulation.replicate", "positive", "actual_maxdist", "actual_avgdist")  
  FT_maxdist <- dfb %>% group_by(pop_no, model, gene.flow.intensity.per.generation, qpAdm.setup, simulation.replicate, actual_maxdist, .drop=F) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(pop_no, model, gene.flow.intensity.per.generation, qpAdm.setup, simulation.replicate, actual_maxdist) %>% mutate(n = replace_na(n, 0))
  FT_avgdist <- dfb %>% group_by(pop_no, model, gene.flow.intensity.per.generation, qpAdm.setup, simulation.replicate, actual_avgdist, .drop=F) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(pop_no, model, gene.flow.intensity.per.generation, qpAdm.setup, simulation.replicate, actual_avgdist) %>% mutate(n = replace_na(n, 0))
  confusion_maxdist <- dfb %>% group_by(pop_no, model, gene.flow.intensity.per.generation, qpAdm.setup, simulation.replicate, positive, actual_maxdist, .drop=F) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(pop_no, model, gene.flow.intensity.per.generation, qpAdm.setup, simulation.replicate, positive, actual_maxdist) %>% mutate(n = replace_na(n, 0))
  confusion_avgdist <- dfb %>% group_by(pop_no, model, gene.flow.intensity.per.generation, qpAdm.setup, simulation.replicate, positive, actual_avgdist, .drop=F) %>% 
    summarise(n = n()) %>% ungroup() %>% complete(pop_no, model, gene.flow.intensity.per.generation, qpAdm.setup, simulation.replicate, positive, actual_avgdist) %>% mutate(n = replace_na(n, 0))
  results <- FT_maxdist[,c(1:5)]
  results <- results[!duplicated(results[, c("pop_no", "model", "gene.flow.intensity.per.generation", "qpAdm.setup", "simulation.replicate")]), ]
  results$number <- i
  results <- left_join(results, fit_criteria, by="number")
  results$odds_maxdist <- filter(FT_maxdist, actual_maxdist=="optimal")$n/filter(FT_maxdist, actual_maxdist=="non-optimal")$n
  results$odds_avgdist <- filter(FT_avgdist, actual_avgdist=="optimal")$n/filter(FT_avgdist, actual_avgdist=="non-optimal")$n
  results$FPR_maxdist <- filter(confusion_maxdist, positive=="fitting", actual_maxdist=="non-optimal")$n/(filter(confusion_maxdist, positive=="fitting", actual_maxdist=="non-optimal")$n + filter(confusion_maxdist, positive=="non-fitting", actual_maxdist=="non-optimal")$n)
  results$FPR_avgdist <- filter(confusion_avgdist, positive=="fitting", actual_avgdist=="non-optimal")$n/(filter(confusion_avgdist, positive=="fitting", actual_avgdist=="non-optimal")$n + filter(confusion_avgdist, positive=="non-fitting", actual_avgdist=="non-optimal")$n)
  results$FNR_maxdist <- filter(confusion_maxdist, positive=="non-fitting", actual_maxdist=="optimal")$n/(filter(confusion_maxdist, positive=="fitting", actual_maxdist=="optimal")$n + filter(confusion_maxdist, positive=="non-fitting", actual_maxdist=="optimal")$n)
  results$FNR_avgdist <- filter(confusion_avgdist, positive=="non-fitting", actual_avgdist=="optimal")$n/(filter(confusion_avgdist, positive=="fitting", actual_avgdist=="optimal")$n + filter(confusion_avgdist, positive=="non-fitting", actual_avgdist=="optimal")$n)
  results$FDR_maxdist <- filter(confusion_maxdist, positive=="fitting", actual_maxdist=="non-optimal")$n/(filter(confusion_maxdist, positive=="fitting", actual_maxdist=="optimal")$n + filter(confusion_maxdist, positive=="fitting", actual_maxdist=="non-optimal")$n)
  results$FDR_avgdist <- filter(confusion_avgdist, positive=="fitting", actual_avgdist=="non-optimal")$n/(filter(confusion_avgdist, positive=="fitting", actual_avgdist=="optimal")$n + filter(confusion_avgdist, positive=="fitting", actual_avgdist=="non-optimal")$n)
  outtab <- rbind(outtab, results)
}
write.table(outtab, file="FDR_per_protocol_per_replicate.txt", row.names=F, sep='\t')


###PLOTTING FIGURES###

###per-replicate FDR results for selected composite model feasibility criteria (no. 28 and 3):
df <- filter(read_delim(file="FDR_per_replicate.txt", delim="\t"), number==28 | number==3)
df$fit_criteria3 <- gsub("1. simpler models not checked", "s.m.n.c.", df$fit_criteria3)
df$fit_criteria3 <- gsub("2. 'trailing' simpler models rejected", "t.s.m.r.", df$fit_criteria3)
df$variables <- paste0(df$fit_criteria3, ";", df$fit_criteria1, ";", df$fit_criteria2, ";", df$gene.flow.intensity.per.generation, ";", df$model)
df$variables <- factor(df$variables, levels=c("s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;4-way"))

###plotting:
stat.test2 <- df %>% pairwise_wilcox_test(FDR_maxdist ~ variables, paired=F)
stat.test2a <- stat.test2 %>% relocate(group2, .before = group1)
names(stat.test2a)[names(stat.test2a)=='group1'] <- 'group1a'
names(stat.test2a)[names(stat.test2a)=='group2'] <- 'group1'
names(stat.test2a)[names(stat.test2a)=='group1a'] <- 'group2'
stat.test3 <- rbind(stat.test2, stat.test2a)
stat.test3$group1 <- factor(stat.test3$group1, levels=c("s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;4-way"))
stat.test3$group2 <- factor(stat.test3$group2, levels=c("s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-4;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-4 to 10-3;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-5 to 10-2;4-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;2-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;3-way", "s.m.n.c.;% (0,1);p>=0.01;10-3 to 10-2;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-4;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-4 to 10-3;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-5 to 10-2;4-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;2-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;3-way", "t.s.m.r.;%+/-2SE (0,1);p>=0.05;10-3 to 10-2;4-way"))

p1 <- ggplot(df, aes(x=variables, y=odds_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.04), limits=c(0,0.22))+
  ggtitle("a. Ratio of optimal to non-optimal models in a set of tested models; based on max. ST distance and min. STS angle")+
  ylab('pre-study odds')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=28, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p2 <- ggplot(df, aes(x=variables, y=FPR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.02), limits=c(0,0.07))+
  ggtitle("b. False positive rate or type I error probability")+
  ylab('FPR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=28, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p3 <- ggplot(df, aes(x=variables, y=FNR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.1), limits=c(0.55,1))+
  ggtitle("c. False negative rate = 1 - power")+
  ylab('FNR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=28, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p4 <- ggplot(df, aes(x=variables, y=FDR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("d. False discovery rate = FPR/(optimal/non-optimal - FNR * optimal/non-optimal + FPR)")+
  ylab('FDR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=28, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p5 <- ggplot(stat.test3, aes(group1, group2, fill=p.adj))+
  geom_tile()+
  scale_fill_gradientn(colours = pals::parula(20), breaks=seq(0,1,0.1), na.value='cornsilk')+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_y_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  ggtitle("e. Wilcoxon tests for differences in FDR, adjusted p-values")+
  guides(fill = guide_colorbar(title = 'p-value\n', order=1, barwidth=4, barheight=48))+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=28, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=28, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='right',
        panel.background=element_rect(fill = 'cornsilk'),
        panel.grid=element_blank(),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

pall <- plot_grid(p1, p2, p3, p4, p5, nrow=5, ncol=1, align = "v", rel_heights=c(1,1,1,1,3.5))


###This is Figure 7###
pdf("FDR_per_replicate_criteria28_3.pdf", width = 36, height = 60)
print(annotate_figure(pall, top=text_grob('results stratified by SSL sparsity and model complexity; results shown for two composite feasibility criteria\n', color="black", face="bold", size=50)))
dev.off()



###per-protocol and per-replicate FDR results for seven selected composite model feasibility criteria:

###model feasibility criterion no. 32:
df <- filter(read_delim(file="FDR_per_protocol_per_replicate.txt", delim="\t"), number==32)
df$variables <- paste0(df$gene.flow.intensity.per.generation, ";", df$model, ";", df$pop_no, ";", df$qpAdm.setup)
df$variables <- factor(df$variables, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                              "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###Wilcoxon test for comparing distributions (all vs. all):
stat.test2 <- df %>% pairwise_wilcox_test(FDR_maxdist ~ variables, paired=F)
stat.test2a <- stat.test2 %>% relocate(group2, .before = group1)
names(stat.test2a)[names(stat.test2a)=='group1'] <- 'group1a'
names(stat.test2a)[names(stat.test2a)=='group2'] <- 'group1'
names(stat.test2a)[names(stat.test2a)=='group1a'] <- 'group2'
stat.test3 <- rbind(stat.test2, stat.test2a)
stat.test3$group1 <- factor(stat.test3$group1, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))
stat.test3$group2 <- factor(stat.test3$group2, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###plotting:
p1 <- ggplot(df, aes(x=variables, y=odds_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.05), limits=c(0,0.32))+
  ggtitle("a. Ratio of optimal to non-optimal models in a set of tested models; based on max. ST distance and min. STS angle")+
  ylab('pre-study odds')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p2 <- ggplot(df, aes(x=variables, y=FPR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.05), limits=c(0,0.32))+
  ggtitle("b. False positive rate or type I error probability")+
  ylab('FPR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p3 <- ggplot(df, aes(x=variables, y=FNR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("c. False negative rate = 1 - power")+
  ylab('FNR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p4 <- ggplot(df, aes(x=variables, y=FDR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("d. False discovery rate = FPR/(optimal/non-optimal - FNR * optimal/non-optimal + FPR)")+
  ylab('FDR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p5 <- ggplot(stat.test3, aes(group1, group2, fill=p.adj))+
  geom_tile()+
  scale_fill_gradientn(colours = pals::parula(20), breaks=seq(0,1,0.1), na.value='cornsilk')+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_y_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  ggtitle("e. Wilcoxon tests for differences in FDR, adjusted p-values")+
  guides(fill = guide_colorbar(title = 'p-value\n', order=1, barwidth=4, barheight=48))+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=22, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='right',
        panel.background=element_rect(fill = 'cornsilk'),
        panel.grid=element_blank(),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

pall <- plot_grid(p1, p2, p3, p4, p5, nrow=5, ncol=1, align = "v", rel_heights=c(1,1,1,1,5))


###This is Figure S24a###
pdf("FDR_per_protocol_per_replicate_criterion32.pdf", width = 43, height = 70)
print(annotate_figure(pall, top=text_grob('results stratified by SSL sparsity, model complexity, number of alternative sources per target, and qpAdm protocols\ncomposite feasibility criterion: simpler models not checked, EAF (-0.3,1.3), p-value >= 0.001\n', color="black", face="bold", size=50)))
dev.off()



###model feasibility criterion no. 28:
df <- filter(read_delim(file="FDR_per_protocol_per_replicate.txt", delim="\t"), number==28)
df$variables <- paste0(df$gene.flow.intensity.per.generation, ";", df$model, ";", df$pop_no, ";", df$qpAdm.setup)
df$variables <- factor(df$variables, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                              "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###Wilcoxon test for comparing distributions (all vs. all):
stat.test2 <- df %>% pairwise_wilcox_test(FDR_maxdist ~ variables, paired=F)
stat.test2a <- stat.test2 %>% relocate(group2, .before = group1)
names(stat.test2a)[names(stat.test2a)=='group1'] <- 'group1a'
names(stat.test2a)[names(stat.test2a)=='group2'] <- 'group1'
names(stat.test2a)[names(stat.test2a)=='group1a'] <- 'group2'
stat.test3 <- rbind(stat.test2, stat.test2a)
stat.test3$group1 <- factor(stat.test3$group1, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))
stat.test3$group2 <- factor(stat.test3$group2, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###plotting:
p1 <- ggplot(df, aes(x=variables, y=odds_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.05), limits=c(0,0.32))+
  ggtitle("a. Ratio of optimal to non-optimal models in a set of tested models; based on max. ST distance and min. STS angle")+
  ylab('pre-study odds')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p2 <- ggplot(df, aes(x=variables, y=FPR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.02), limits=c(0,0.14))+
  ggtitle("b. False positive rate or type I error probability")+
  ylab('FPR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p3 <- ggplot(df, aes(x=variables, y=FNR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("c. False negative rate = 1 - power")+
  ylab('FNR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p4 <- ggplot(df, aes(x=variables, y=FDR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("d. False discovery rate = FPR/(optimal/non-optimal - FNR * optimal/non-optimal + FPR)")+
  ylab('FDR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p5 <- ggplot(stat.test3, aes(group1, group2, fill=p.adj))+
  geom_tile()+
  scale_fill_gradientn(colours = pals::parula(20), breaks=seq(0,1,0.1), na.value='cornsilk')+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_y_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  ggtitle("e. Wilcoxon tests for differences in FDR, adjusted p-values")+
  guides(fill = guide_colorbar(title = 'p-value\n', order=1, barwidth=4, barheight=48))+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=22, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='right',
        panel.background=element_rect(fill = 'cornsilk'),
        panel.grid=element_blank(),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

pall <- plot_grid(p1, p2, p3, p4, p5, nrow=5, ncol=1, align = "v", rel_heights=c(1,1,1,1,5))


###This is Figure S24b###
pdf("FDR_per_protocol_per_replicate_criterion28.pdf", width = 43, height = 70)
print(annotate_figure(pall, top=text_grob('results stratified by SSL sparsity, model complexity, number of alternative sources per target, and qpAdm protocols\ncomposite feasibility criterion: simpler models not checked, EAF (0,1), p-value >= 0.01\n', color="black", face="bold", size=50)))
dev.off()



###model feasibility criterion no. 22:
df <- filter(read_delim(file="FDR_per_protocol_per_replicate.txt", delim="\t"), number==22)
df$variables <- paste0(df$gene.flow.intensity.per.generation, ";", df$model, ";", df$pop_no, ";", df$qpAdm.setup)
df$variables <- factor(df$variables, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                              "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###Wilcoxon test for comparing distributions (all vs. all):
stat.test2 <- df %>% pairwise_wilcox_test(FDR_maxdist ~ variables, paired=F)
stat.test2a <- stat.test2 %>% relocate(group2, .before = group1)
names(stat.test2a)[names(stat.test2a)=='group1'] <- 'group1a'
names(stat.test2a)[names(stat.test2a)=='group2'] <- 'group1'
names(stat.test2a)[names(stat.test2a)=='group1a'] <- 'group2'
stat.test3 <- rbind(stat.test2, stat.test2a)
stat.test3$group1 <- factor(stat.test3$group1, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))
stat.test3$group2 <- factor(stat.test3$group2, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###plotting:
p1 <- ggplot(df, aes(x=variables, y=odds_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.05), limits=c(0,0.32))+
  ggtitle("a. Ratio of optimal to non-optimal models in a set of tested models; based on max. ST distance and min. STS angle")+
  ylab('pre-study odds')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p2 <- ggplot(df, aes(x=variables, y=FPR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,0.05,0.01), limits=c(0,0.042))+
  ggtitle("b. False positive rate or type I error probability")+
  ylab('FPR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p3 <- ggplot(df, aes(x=variables, y=FNR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("c. False negative rate = 1 - power")+
  ylab('FNR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p4 <- ggplot(df, aes(x=variables, y=FDR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("d. False discovery rate = FPR/(optimal/non-optimal - FNR * optimal/non-optimal + FPR)")+
  ylab('FDR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p5 <- ggplot(stat.test3, aes(group1, group2, fill=p.adj))+
  geom_tile()+
  scale_fill_gradientn(colours = pals::parula(20), breaks=seq(0,1,0.1), na.value='cornsilk')+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_y_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  ggtitle("e. Wilcoxon tests for differences in FDR, adjusted p-values")+
  guides(fill = guide_colorbar(title = 'p-value\n', order=1, barwidth=4, barheight=48))+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=22, face="bold"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=36, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='right',
        panel.background=element_rect(fill = 'cornsilk'),
        panel.grid=element_blank(),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

pall <- plot_grid(p1, p2, p3, p4, p5, nrow=5, ncol=1, align = "v", rel_heights=c(1,1,1,1,5))


###This is Figure S24c###
pdf("FDR_per_protocol_per_replicate_criterion22.pdf", width = 43, height = 70)
print(annotate_figure(pall, top=text_grob('results stratified by SSL sparsity, model complexity, number of alternative sources per target, and qpAdm protocols\ncomposite feasibility criterion: simpler models not checked, EAF +/- 2SE (0,1), p-value >= 0.001\n', color="black", face="bold", size=50)))
dev.off()



###model feasibility criterion no. 15###
df <- filter(read_delim(file="FDR_per_protocol_per_replicate.txt", delim="\t"), number==15)
df$variables <- paste0(df$gene.flow.intensity.per.generation, ";", df$model, ";", df$pop_no, ";", df$qpAdm.setup)
df$variables <- factor(df$variables, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                              "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###Wilcoxon tests for comparing distributions (all vs. all) do not work do to the presence of missing data
#stat.test2 <- df %>% pairwise_wilcox_test(FDR_maxdist ~ variables, paired=F)
#stat.test2a <- stat.test2 %>% relocate(group2, .before = group1)
#names(stat.test2a)[names(stat.test2a)=='group1'] <- 'group1a'
#names(stat.test2a)[names(stat.test2a)=='group2'] <- 'group1'
#names(stat.test2a)[names(stat.test2a)=='group1a'] <- 'group2'
#stat.test3 <- rbind(stat.test2, stat.test2a)
#stat.test3$group1 <- factor(stat.test3$group1, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
#                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))
#stat.test3$group2 <- factor(stat.test3$group2, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
#                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###plotting:
p1 <- ggplot(df, aes(x=variables, y=odds_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.05), limits=c(0,0.32))+
  ggtitle("a. Ratio of optimal to non-optimal models in a set of tested models; based on max. ST distance and min. STS angle")+
  ylab('pre-study odds')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p2 <- ggplot(df, aes(x=variables, y=FPR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.01), limits=c(0,0.05))+
  ggtitle("b. False positive rate or type I error probability")+
  ylab('FPR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p3 <- ggplot(df, aes(x=variables, y=FNR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("c. False negative rate = 1 - power")+
  ylab('FNR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p4 <- ggplot(df, aes(x=variables, y=FDR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("d. False discovery rate = FPR/(optimal/non-optimal - FNR * optimal/non-optimal + FPR)")+
  ylab('FDR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

#p5 <- ggplot(stat.test3, aes(group1, group2, fill=p.adj))+
#  geom_tile()+
#  scale_fill_gradientn(colours = pals::parula(20), breaks=seq(0,1,0.1), na.value='cornsilk')+
#  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
#  scale_y_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
#  ggtitle("e. Wilcoxon tests for differences in FDR, adjusted p-values")+
#  guides(fill = guide_colorbar(title = 'p-value\n', order=1, barwidth=4, barheight=48))+
#  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
#        axis.line=element_line(linewidth=1),
#        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
#        axis.text.y=element_text(size=22, face="bold"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_blank(),
#        legend.key=element_rect(fill = 'white'),
#        legend.title=element_text(size=36, face = "bold"),
#        legend.text=element_text(size=28),
#        legend.position='right',
#        panel.background=element_rect(fill = 'cornsilk'),
#        panel.grid=element_blank(),
#        plot.margin = margin(10, 5, 10, 5, 'pt'))

pall <- plot_grid(p1, p2, p3, p4, nrow=4, ncol=1, align = "v", rel_heights=c(1,1,1,1))


###This is Figure S24d###
pdf("FDR_per_protocol_per_replicate_criterion15.pdf", width = 43, height = 32)
print(annotate_figure(pall, top=text_grob('results stratified by SSL sparsity, model complexity, number of alternative sources per target, and qpAdm protocols\ncomposite feasibility criterion: "trailing" simpler models rejected, EAF (-0.3,1.3), p-value >= 0.001\n', color="black", face="bold", size=50)))
dev.off()



###model feasibility criterion no. 9:
df <- filter(read_delim(file="FDR_per_protocol_per_replicate.txt", delim="\t"), number==9)
df$variables <- paste0(df$gene.flow.intensity.per.generation, ";", df$model, ";", df$pop_no, ";", df$qpAdm.setup)
df$variables <- factor(df$variables, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                              "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###Wilcoxon tests for comparing distributions (all vs. all) do not work do to the presence of missing data
#stat.test2 <- df %>% pairwise_wilcox_test(FDR_maxdist ~ variables, paired=F)
#stat.test2a <- stat.test2 %>% relocate(group2, .before = group1)
#names(stat.test2a)[names(stat.test2a)=='group1'] <- 'group1a'
#names(stat.test2a)[names(stat.test2a)=='group2'] <- 'group1'
#names(stat.test2a)[names(stat.test2a)=='group1a'] <- 'group2'
#stat.test3 <- rbind(stat.test2, stat.test2a)
#stat.test3$group1 <- factor(stat.test3$group1, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
#                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))
#stat.test3$group2 <- factor(stat.test3$group2, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
#                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

p1 <- ggplot(df, aes(x=variables, y=odds_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.05), limits=c(0,0.32))+
  ggtitle("a. Ratio of optimal to non-optimal models in a set of tested models; based on max. ST distance and min. STS angle")+
  ylab('pre-study odds')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p2 <- ggplot(df, aes(x=variables, y=FPR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.004), limits=c(0,0.020))+
  ggtitle("b. False positive rate or type I error probability")+
  ylab('FPR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p3 <- ggplot(df, aes(x=variables, y=FNR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("c. False negative rate = 1 - power")+
  ylab('FNR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p4 <- ggplot(df, aes(x=variables, y=FDR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5))+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("d. False discovery rate = FPR/(optimal/non-optimal - FNR * optimal/non-optimal + FPR)")+
  ylab('FDR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

#p5 <- ggplot(stat.test3, aes(group1, group2, fill=p.adj))+
#  geom_tile()+
#  scale_fill_gradientn(colours = pals::parula(20), breaks=seq(0,1,0.1), na.value='cornsilk')+
#  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
#  scale_y_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
#  ggtitle("e. Wilcoxon tests for differences in FDR, adjusted p-values")+
#  guides(fill = guide_colorbar(title = 'p-value\n', order=1, barwidth=4, barheight=48))+
#  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
#        axis.line=element_line(linewidth=1),
#        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
#        axis.text.y=element_text(size=22, face="bold"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_blank(),
#        legend.key=element_rect(fill = 'white'),
#        legend.title=element_text(size=36, face = "bold"),
#        legend.text=element_text(size=28),
#        legend.position='right',
#        panel.background=element_rect(fill = 'cornsilk'),
#        panel.grid=element_blank(),
#        plot.margin = margin(10, 5, 10, 5, 'pt'))

pall <- plot_grid(p1, p2, p3, p4, nrow=4, ncol=1, align = "v", rel_heights=c(1,1,1,1))


###This is Figure S24e###
pdf("FDR_per_protocol_per_replicate_criterion09.pdf", width = 43, height = 32)
print(annotate_figure(pall, top=text_grob('results stratified by SSL sparsity, model complexity, number of alternative sources per target, and qpAdm protocols\ncomposite feasibility criterion: "trailing" simpler models rejected, EAF (0,1), p-value >= 0.01\n', color="black", face="bold", size=50)))
dev.off()



###model feasibility criterion no. 3:
df <- filter(read_delim(file="FDR_per_protocol_per_replicate.txt", delim="\t"), number==3)
df$variables <- paste0(df$gene.flow.intensity.per.generation, ";", df$model, ";", df$pop_no, ";", df$qpAdm.setup)
df$variables <- factor(df$variables, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
                                              "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###Wilcoxon tests for comparing distributions (all vs. all) do not work do to the presence of missing data
#stat.test2 <- df %>% pairwise_wilcox_test(FDR_maxdist ~ variables, paired=F)
#stat.test2a <- stat.test2 %>% relocate(group2, .before = group1)
#names(stat.test2a)[names(stat.test2a)=='group1'] <- 'group1a'
#names(stat.test2a)[names(stat.test2a)=='group2'] <- 'group1'
#names(stat.test2a)[names(stat.test2a)=='group1a'] <- 'group2'
#stat.test3 <- rbind(stat.test2, stat.test2a)
#stat.test3$group1 <- factor(stat.test3$group1, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
#                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))
#stat.test3$group2 <- factor(stat.test3$group2, levels=c("10-5 to 10-4;2-way;13 demes;dist. rot.", "10-5 to 10-4;2-way;13 demes;prox. rot.", "10-5 to 10-4;2-way;13 demes;dist. non-rot.", "10-5 to 10-4;2-way;13 demes;prox. non-rot.", "10-5 to 10-4;2-way;18 demes;dist. rot.", "10-5 to 10-4;2-way;18 demes;prox. rot.", "10-5 to 10-4;2-way;18 demes;dist. non-rot.", "10-5 to 10-4;2-way;18 demes;prox. non-rot.", "10-5 to 10-4;3-way;13 demes;dist. rot.", "10-5 to 10-4;3-way;13 demes;prox. rot.", "10-5 to 10-4;3-way;13 demes;dist. non-rot.", "10-5 to 10-4;3-way;13 demes;prox. non-rot.", "10-5 to 10-4;3-way;18 demes;dist. rot.", "10-5 to 10-4;3-way;18 demes;prox. rot.", "10-5 to 10-4;3-way;18 demes;dist. non-rot.", "10-5 to 10-4;3-way;18 demes;prox. non-rot.", "10-5 to 10-4;4-way;13 demes;dist. rot.", "10-5 to 10-4;4-way;13 demes;prox. rot.", "10-5 to 10-4;4-way;13 demes;dist. non-rot.", "10-5 to 10-4;4-way;13 demes;prox. non-rot.", "10-5 to 10-4;4-way;18 demes;dist. rot.", "10-5 to 10-4;4-way;18 demes;prox. rot.", "10-5 to 10-4;4-way;18 demes;dist. non-rot.", "10-5 to 10-4;4-way;18 demes;prox. non-rot.", "10-4 to 10-3;2-way;13 demes;dist. rot.", "10-4 to 10-3;2-way;13 demes;prox. rot.", "10-4 to 10-3;2-way;13 demes;dist. non-rot.", "10-4 to 10-3;2-way;13 demes;prox. non-rot.", "10-4 to 10-3;2-way;18 demes;dist. rot.", "10-4 to 10-3;2-way;18 demes;prox. rot.", "10-4 to 10-3;2-way;18 demes;dist. non-rot.", "10-4 to 10-3;2-way;18 demes;prox. non-rot.", "10-4 to 10-3;3-way;13 demes;dist. rot.", "10-4 to 10-3;3-way;13 demes;prox. rot.", "10-4 to 10-3;3-way;13 demes;dist. non-rot.", "10-4 to 10-3;3-way;13 demes;prox. non-rot.", "10-4 to 10-3;3-way;18 demes;dist. rot.", "10-4 to 10-3;3-way;18 demes;prox. rot.", "10-4 to 10-3;3-way;18 demes;dist. non-rot.", "10-4 to 10-3;3-way;18 demes;prox. non-rot.", "10-4 to 10-3;4-way;13 demes;dist. rot.", "10-4 to 10-3;4-way;13 demes;prox. rot.", "10-4 to 10-3;4-way;13 demes;dist. non-rot.", "10-4 to 10-3;4-way;13 demes;prox. non-rot.", "10-4 to 10-3;4-way;18 demes;dist. rot.", "10-4 to 10-3;4-way;18 demes;prox. rot.", "10-4 to 10-3;4-way;18 demes;dist. non-rot.", "10-4 to 10-3;4-way;18 demes;prox. non-rot.", "10-5 to 10-2;2-way;13 demes;dist. rot.", "10-5 to 10-2;2-way;13 demes;prox. rot.", "10-5 to 10-2;2-way;13 demes;dist. non-rot.", 
#                                                        "10-5 to 10-2;2-way;13 demes;prox. non-rot.", "10-5 to 10-2;2-way;18 demes;dist. rot.", "10-5 to 10-2;2-way;18 demes;prox. rot.", "10-5 to 10-2;2-way;18 demes;dist. non-rot.", "10-5 to 10-2;2-way;18 demes;prox. non-rot.", "10-5 to 10-2;3-way;13 demes;dist. rot.", "10-5 to 10-2;3-way;13 demes;prox. rot.", "10-5 to 10-2;3-way;13 demes;dist. non-rot.", "10-5 to 10-2;3-way;13 demes;prox. non-rot.", "10-5 to 10-2;3-way;18 demes;dist. rot.", "10-5 to 10-2;3-way;18 demes;prox. rot.", "10-5 to 10-2;3-way;18 demes;dist. non-rot.", "10-5 to 10-2;3-way;18 demes;prox. non-rot.", "10-5 to 10-2;4-way;13 demes;dist. rot.", "10-5 to 10-2;4-way;13 demes;prox. rot.", "10-5 to 10-2;4-way;13 demes;dist. non-rot.", "10-5 to 10-2;4-way;13 demes;prox. non-rot.", "10-5 to 10-2;4-way;18 demes;dist. rot.", "10-5 to 10-2;4-way;18 demes;prox. rot.", "10-5 to 10-2;4-way;18 demes;dist. non-rot.", "10-5 to 10-2;4-way;18 demes;prox. non-rot.", "10-3 to 10-2;2-way;13 demes;dist. rot.", "10-3 to 10-2;2-way;13 demes;prox. rot.", "10-3 to 10-2;2-way;13 demes;dist. non-rot.", "10-3 to 10-2;2-way;13 demes;prox. non-rot.", "10-3 to 10-2;2-way;18 demes;dist. rot.", "10-3 to 10-2;2-way;18 demes;prox. rot.", "10-3 to 10-2;2-way;18 demes;dist. non-rot.", "10-3 to 10-2;2-way;18 demes;prox. non-rot.", "10-3 to 10-2;3-way;13 demes;dist. rot.", "10-3 to 10-2;3-way;13 demes;prox. rot.", "10-3 to 10-2;3-way;13 demes;dist. non-rot.", "10-3 to 10-2;3-way;13 demes;prox. non-rot.", "10-3 to 10-2;3-way;18 demes;dist. rot.", "10-3 to 10-2;3-way;18 demes;prox. rot.", "10-3 to 10-2;3-way;18 demes;dist. non-rot.", "10-3 to 10-2;3-way;18 demes;prox. non-rot.", "10-3 to 10-2;4-way;13 demes;dist. rot.", "10-3 to 10-2;4-way;13 demes;prox. rot.", "10-3 to 10-2;4-way;13 demes;dist. non-rot.", "10-3 to 10-2;4-way;13 demes;prox. non-rot.", "10-3 to 10-2;4-way;18 demes;dist. rot.", "10-3 to 10-2;4-way;18 demes;prox. rot.", "10-3 to 10-2;4-way;18 demes;dist. non-rot.", "10-3 to 10-2;4-way;18 demes;prox. non-rot."))

###plotting:
p1 <- ggplot(df, aes(x=variables, y=odds_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5), drop=F)+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.001,0.001), breaks=seq(0,1,0.05), limits=c(0,0.32))+
  ggtitle("a. Ratio of optimal to non-optimal models in a set of tested models; based on max. ST distance and min. STS angle")+
  ylab('pre-study odds')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p2 <- ggplot(df, aes(x=variables, y=FPR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5), drop=F)+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.0005,0.0005), breaks=seq(0,1,0.004), limits=c(0,0.016))+
  ggtitle("b. False positive rate or type I error probability")+
  ylab('FPR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p3 <- ggplot(df, aes(x=variables, y=FNR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5), drop=F)+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("c. False negative rate = 1 - power")+
  ylab('FNR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

p4 <- ggplot(df, aes(x=variables, y=FDR_maxdist))+
  geom_violin(aes(fill=variables), scale="width", trim=T, draw_quantiles=c(0.5), drop=F)+
  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
  scale_fill_manual(values=c(stepped(24), stepped(24), stepped(24), stepped(24)), drop=F)+
  scale_y_continuous(expand=c(0.01,0.01), breaks=seq(0,1,0.2), limits=c(0,1))+
  ggtitle("d. False discovery rate = FPR/(optimal/non-optimal - FNR * optimal/non-optimal + FPR)")+
  ylab('FDR')+
  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
        axis.line=element_line(linewidth=1),
        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
        axis.text.y=element_text(size=30, face='bold'),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=46, face='bold'),
        legend.key=element_rect(fill = 'white'),
        legend.title=element_text(size=28, face = "bold"),
        legend.text=element_text(size=28),
        legend.position='none',
        panel.background=element_rect(fill = 'white'),
        panel.grid.major.x=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.5,  linetype=3),
        plot.margin = margin(10, 5, 10, 5, 'pt'))

#p5 <- ggplot(stat.test3, aes(group1, group2, fill=p.adj))+
#  geom_tile()+
#  scale_fill_gradientn(colours = pals::parula(20), breaks=seq(0,1,0.1), na.value='cornsilk')+
#  scale_x_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
#  scale_y_discrete(guide = guide_axis_nested(delim = ";", inv=T))+
#  ggtitle("e. Wilcoxon tests for differences in FDR, adjusted p-values")+
#  guides(fill = guide_colorbar(title = 'p-value\n', order=1, barwidth=4, barheight=48))+
#  theme(plot.title=element_text(size=46, face='bold', color="purple3"),
#        axis.line=element_line(linewidth=1),
#        axis.text.x=element_text(size=22, face="bold", angle=90, hjust=0, vjust=0.5),
#        axis.text.y=element_text(size=22, face="bold"),
#        axis.title.x=element_blank(),
#        axis.title.y=element_blank(),
#        legend.key=element_rect(fill = 'white'),
#        legend.title=element_text(size=36, face = "bold"),
#        legend.text=element_text(size=28),
#        legend.position='right',
#        panel.background=element_rect(fill = 'cornsilk'),
#        panel.grid=element_blank(),
#        plot.margin = margin(10, 5, 10, 5, 'pt'))

pall <- plot_grid(p1, p2, p3, p4, nrow=4, ncol=1, align = "v", rel_heights=c(1,1,1,1))


###This is Figure S24f###
pdf("FDR_per_protocol_per_replicate_criterion03.pdf", width = 43, height = 32)
print(annotate_figure(pall, top=text_grob('results stratified by SSL sparsity, model complexity, number of alternative sources per target, and qpAdm protocols\ncomposite feasibility criterion: "trailing" simpler models rejected, EAF +/- 2SE (0,1), p-value >= 0.05\n', color="black", face="bold", size=50)))
dev.off()

