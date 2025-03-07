library(readr)
library(dplyr)
library(ggplot2)
library(pals)
library(ggpubr)
library(gridExtra)
library(ggsci)
library(rstatix)
library(plyr)

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm/systematic_sampling_Olas_results')
df <- read_delim(file="input_data_Figure_10.txt", delim="\t", trim_ws=T)

###summarizing results for the distal rotating qpAdm protocol:
df %>% count(closest_sources_number_in_experiment, complexity_rotating) %>%
  group_by(closest_sources_number_in_experiment) %>%
  mutate(freq = n / sum(n)) -> summary.data1a

df %>% count(closest_sources_number_in_experiment, complexity_rotating, simulation_replicate) %>%
  group_by(closest_sources_number_in_experiment, simulation_replicate) %>%
  mutate(freq = n / sum(n)) -> summary.data1b

summary.data1c <- summary.data1b %>% group_by(closest_sources_number_in_experiment, complexity_rotating) %>% summarise_at(vars(freq), list(SD=sd, mean=mean))
summary.data1d <- merge(summary.data1a, summary.data1c, by = c("closest_sources_number_in_experiment", "complexity_rotating"))


###making barplots for the distal rotating qpAdm protocol:
p1a <- ggplot(summary.data1d, aes(x=complexity_rotating, y=freq))+ 
  geom_bar(aes(fill=complexity_rotating, color=complexity_rotating), stat="identity", position="dodge", width=.8)+
  geom_errorbar(aes(ymin=freq-SD, ymax=freq+SD), width=0, linewidth=1, colour="grey30")+
  geom_text(data=summary.data1d, aes(x = complexity_rotating, y = 0.95, label=n), color="violetred", fontface=2, size=7)+
  scale_fill_manual(values=pal_frontiers("default")(4))+
  scale_color_manual(values=pal_frontiers("default")(4))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(0,0.25,0.5,0.75,1))+
  coord_cartesian(ylim = c(0,1))+
  ggtitle("10-5 to 10-2 landscape, 10 sim. repl., distal rotating protocol;\nexperiments grouped by number of target's nearest neighbors available")+
  xlab('model complexity level at which experiment ends')+
  ylab('% experiments ending')+
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
  facet_grid(. ~ closest_sources_number_in_experiment)


###pairwise Wilcoxon tests for the distal rotating qpAdm protocol:
df %>% count(complexity_rotating, closest_sources_number_rotating) %>%
  group_by(complexity_rotating, closest_sources_number_rotating) -> summary.data1

stat.test1 <- df %>% filter(complexity_rotating!="not finished") %>% group_by(complexity_rotating) %>%  pairwise_wilcox_test(p_rotating ~ closest_sources_number_rotating, paired=F)
stat.test1 <- stat.test1 %>% add_y_position()


###making violin plots for the distal rotating qpAdm protocol:
p1 <- ggplot(filter(df, complexity_rotating!="not finished"), aes(x=factor(closest_sources_number_rotating, levels=c(0,1,2,3,4)), y=p_rotating))+
  geom_violin(aes(fill=factor(closest_sources_number_rotating, levels=c(0,1,2,3,4))), scale="width", trim=T, linewidth=1, color="grey20", draw_quantiles = c(0.5))+
  geom_jitter(height = 0, width = 0.2, alpha=0.5)+
  scale_x_discrete(drop=F, na.translate=F)+
  geom_text(data=filter(summary.data1, complexity_rotating!="not finished"), aes(x = closest_sources_number_rotating, y = -0.02, label=n), color="violetred", nudge_x=1, fontface=2, size=8)+
  scale_fill_manual(values=pal_flatui("default")(5), drop=F)+
  scale_y_continuous(expand = c(0.01,0.01), breaks=c(0,0.25,0.5,0.75,1), limits=c(-0.05,1.17))+
  ggtitle("grouped by model complexity level at which experiment ends")+
  xlab("number of target's nearest neighbors in a model")+
  ylab('max. p-value per experiment & complexity level')+
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
  facet_grid(. ~ complexity_rotating)+
  stat_pvalue_manual(stat.test1, label = "p.adj.signif", y.position=c(1, 1.05, 1.1, 1.15, 1), hide.ns=T, tip.length = 0.01, size=10, bracket.size=1)




###summarizing results for the distal non-rotating qpAdm protocol:
df %>% count(closest_sources_number_in_experiment, complexity_non.rotating) %>%
  group_by(closest_sources_number_in_experiment) %>%
  mutate(freq = n / sum(n)) -> summary.data2a

df %>% count(closest_sources_number_in_experiment, complexity_non.rotating, simulation_replicate) %>%
  group_by(closest_sources_number_in_experiment, simulation_replicate) %>%
  mutate(freq = n / sum(n)) -> summary.data2b

summary.data2c <- summary.data2b %>% group_by(closest_sources_number_in_experiment, complexity_non.rotating) %>% summarise_at(vars(freq), list(SD=sd, mean=mean))
summary.data2d <- merge(summary.data2a, summary.data2c, by = c("closest_sources_number_in_experiment", "complexity_non.rotating"))


###making barplots for the distal non-rotating qpAdm protocol:
p2a <- ggplot(summary.data2d, aes(x=complexity_non.rotating, y=freq))+ 
  geom_bar(aes(fill=complexity_non.rotating, color=complexity_non.rotating), stat="identity", position="dodge", width=.8)+
  geom_errorbar(aes(ymin=freq-SD, ymax=freq+SD), width=0, linewidth=1, colour="grey30")+
  geom_text(data=summary.data2d, aes(x = complexity_non.rotating, y = 0.95, label=n), color="violetred", fontface=2, size=7)+
  scale_fill_manual(values=pal_frontiers("default")(4))+
  scale_color_manual(values=pal_frontiers("default")(4))+
  scale_y_continuous(expand=c(0.01,0.01), breaks=c(0,0.25,0.5,0.75,1))+
  coord_cartesian(ylim = c(0,1))+
  ggtitle("10-5 to 10-2 landscape, 10 sim. repl., distal non-rotating protocol;\nexperiments grouped by number of target's nearest neighbors available")+
  xlab('model complexity level at which experiment ends')+
  ylab('% experiments ending')+
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
  facet_grid(. ~ closest_sources_number_in_experiment)


###pairwise Wilcoxon tests for the distal non-rotating qpAdm protocol:
df %>% count(complexity_non.rotating, closest_sources_number_non.rotating) %>%
  group_by(complexity_non.rotating, closest_sources_number_non.rotating) -> summary.data2

stat.test2 <- df %>% filter(complexity_non.rotating!="not finished") %>% group_by(complexity_non.rotating) %>%  pairwise_wilcox_test(p_non.rotating ~ closest_sources_number_non.rotating, paired=F)
stat.test2 <- stat.test2 %>% add_y_position()


###making violin plots for the distal non-rotating qpAdm protocol:
p2 <- ggplot(filter(df, complexity_non.rotating!="not finished"), aes(x=factor(closest_sources_number_non.rotating, levels=c(0,1,2,3,4)), y=p_non.rotating))+
  geom_violin(aes(fill=factor(closest_sources_number_non.rotating, levels=c(0,1,2,3,4))), scale="width", trim=T, linewidth=1, color="grey20", draw_quantiles = c(0.5))+
  geom_jitter(height = 0, width = 0.2, alpha=0.5)+
  scale_x_discrete(drop=F, na.translate=F)+
  geom_text(data=filter(summary.data2, complexity_non.rotating!="not finished"), aes(x = closest_sources_number_non.rotating, y = -0.02, label=n), color="violetred", nudge_x=1, fontface=2, size=8)+
  scale_fill_manual(values=pal_flatui("default")(5), drop=F)+
  scale_y_continuous(expand = c(0.01,0.01), breaks=c(0,0.25,0.5,0.75,1), limits=c(-0.05,1.17))+
  ggtitle("grouped by model complexity level at which experiment ends")+
  xlab("number of target's nearest neighbors in a model")+
  ylab('max. p-value per experiment & complexity level')+
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
  facet_grid(. ~ complexity_non.rotating)+
  stat_pvalue_manual(stat.test2, label = "p.adj.signif", y.position=c(1, 1, 1.05, 1.1, 1), hide.ns=T, tip.length = 0.01, size=10, bracket.size=1)


###This is Figure 10###
pdf("systematic_sampling_2D_plots_best_models_p4_simplified_v2_SD.pdf", width = 30, height = 15)
grid.arrange(arrangeGrob(p1a, p2a, p1, p2, nrow=2, ncol=2, heights=c(3, 4)))
dev.off()





