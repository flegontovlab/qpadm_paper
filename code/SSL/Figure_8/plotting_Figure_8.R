library(readr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(pals)
library(ggpubr)
library(ggh4x)

setwd('./Isolation_by_distance/stepping-stone_simulations/qpAdm/barplots_experiments_basic_1-4way_replicates/')

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

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

###defining colors:
lab=c("best available 1-way models NOT rejected"="darkolivegreen",
      "best available 1-way models rejected"="firebrick",
      "best available 2-way models NOT rejected, all simpler models rejected"="olivedrab1",
      "best available 2-way models rejected, all simpler models rejected"="darkorange1",
      "best available 3-way models NOT rejected, all simpler models rejected"="springgreen",
      "best available 3-way models rejected, all simpler models rejected"="violetred1",
      "best available 4-way models NOT rejected, all simpler models rejected"="lightseagreen",
      "best available 4-way models rejected, all simpler models rejected"="magenta3")
leg=c("best available 1-way models NOT rejected",
      "best available 1-way models rejected",
      "best available 2-way models NOT rejected, all simpler models rejected",
      "best available 2-way models rejected, all simpler models rejected",
      "best available 3-way models NOT rejected, all simpler models rejected",
      "best available 3-way models rejected, all simpler models rejected",
      "best available 4-way models NOT rejected, all simpler models rejected",
      "best available 4-way models rejected, all simpler models rejected")


###based on avg. distance and last 15 feasibility criteria: part 2
gene_flows <- c("10-5 to 10-4", "10-4 to 10-3", "10-3 to 10-2", "10-5 to 10-2")
for (gf in gene_flows) {
  df <- filter(read_delim("data_for_Figure_8_v3.txt", delim="\t", col_types = "dcccccdccccd", trim_ws=T), metric!="false discovery rate" & metric!="positive rate" & metric!="percentage of experiments ending per complexity level" 
               & simpler_models_treatment=="simpler models not checked" & gene_flow==gf & fit_criteria2!="0.01" & fit_criteria2!="0.1" & fit_criteria2!="0.001,   0.05")
  ending2 <- filter(read_delim("data_for_Figure_8_v3.txt", delim="\t", col_types = "dcccccdccccd", trim_ws=T), metric=="percentage of experiments ending per complexity level"
                   & simpler_models_treatment=="simpler models not checked" & gene_flow==gf & fit_criteria2!="0.01" & fit_criteria2!="0.1" & fit_criteria2!="0.001,   0.05")
  fdr <- filter(read_delim("data_for_Figure_8_v3.txt", delim="\t", col_types = "dcccccdccccd", trim_ws=T), metric=="false discovery rate"
                & simpler_models_treatment=="simpler models not checked" & gene_flow==gf & fit_criteria2!="0.01" & fit_criteria2!="0.1" & fit_criteria2!="0.001,   0.05")
  df2 <- filter(df, basic_metric!="max. distance" & basic_metric!="avg. distance" & basic_metric!="dist_to_corner_maxdist" & basic_metric!="dist_to_corner_avgdist_diffsqrt8")
  #ending2 <- filter(ending, basic_metric!="max. distance" & basic_metric!="avg. distance" & basic_metric!="dist_to_corner_maxdist")
  fdr2 <- filter(fdr, basic_metric!="max. distance" & basic_metric!="avg. distance" & basic_metric!="dist_to_corner_maxdist" & basic_metric!="dist_to_corner_avgdist_diffsqrt8")

  df2 <- data_summary(df2, varname="value", groupnames=c("no_outgroups_proxies", "basic_metric", "no_sources", "metric", "gene_flow", "qpAdm_setup", "simpler_models_treatment", "fit_criteria1", "fit_criteria2"))
  ending2 <- data_summary(ending2, varname="value", groupnames=c("no_outgroups_proxies", "basic_metric", "no_sources", "metric", "gene_flow", "qpAdm_setup", "simpler_models_treatment", "fit_criteria1", "fit_criteria2"))
  fdr2 <- data_summary(fdr2, varname="value", groupnames=c("no_outgroups_proxies", "basic_metric", "no_sources", "metric", "gene_flow", "qpAdm_setup", "simpler_models_treatment", "fit_criteria1", "fit_criteria2"))
  
  p1<-ggplot(ending2, aes(x=factor(fit_criteria2, level=c("0.001", "0.05", "0.5", "0.001,   0.5")), y=value))+ 
    geom_bar(position="stack", stat="identity", fill="coral", width=.8)+
    geom_errorbar(aes(ymin=value-3*sd/sqrt(10), ymax=value+3*sd/sqrt(10)), width=0, linewidth=1.5, colour="grey30")+
    scale_y_continuous(expand=c(0.01,0.01), breaks=c(0.2,0.4,0.6,0.8,1), minor_breaks = seq(0,1,0.1))+
    coord_cartesian(ylim = c(0,1))+
    xlab('')+
    ylab('frac. of experiments ending per complexity level')+
    ggtitle('results grouped by conditions on admixture proportions, qpAdm protocols, and sampling density; error bars: 3 SE on 10 sim. iters.')+
    theme(
      panel.ontop=T,
      axis.text.x=element_blank(),
      axis.title.x=element_text(size=44, face='bold', colour="black"), #x-axis title font
      axis.ticks.x=element_line(linewidth=1),
      axis.text.y=element_text(size=36, face='bold'), #y-axis font
      axis.title.y=element_text(size=44, face='bold', colour="black"), #y-axis title font
      axis.ticks.y=element_line(linewidth=1),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.75,  linetype=3),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_line(color="#8ccde3",  linewidth=0.75,  linetype=3),
      panel.background=element_rect(fill=NA),
      legend.position="none",
      legend.key=element_rect(fill = 'white'),
      legend.title=element_text(size=44, face = "bold"),
      legend.text=element_text(size=40),
      plot.title=element_text(size=50, face='bold', color="darkred"),
      plot.margin=margin(t=5, r=5, b=5, l=35, unit="pt"),
      strip.text=element_text(size=40, face='bold'))+
    facet_nested(no_sources ~ factor(fit_criteria1, level=c("% (-0.3,1.3)", "% (0,1)", "%+/-2SE (0,1)")) + factor(qpAdm_setup, levels=c("dist. rot.", "prox. rot.", "dist. non-rot.", "prox. non-rot.")) + no_outgroups_proxies)
  
  p2<-ggplot(fdr2, aes(x=factor(fit_criteria2, level=c("0.001", "0.05", "0.5", "0.001,   0.5")), y=value))+ 
    geom_bar(position="stack", stat="identity", fill="bisque4", width=.8)+
    geom_errorbar(aes(ymin=value-3*sd/sqrt(10), ymax=value+3*sd/sqrt(10)), width=0, linewidth=1.5, colour="grey30")+
    scale_y_continuous(expand=c(0.01,0.01), breaks=c(0.2,0.4,0.6,0.8,1), minor_breaks = seq(0,1,0.1))+
    coord_cartesian(ylim = c(0,1))+
    xlab('')+
    ylab('fraction of experiments with undesirable outcomes')+
    theme(
      panel.ontop=T,
      axis.text.x=element_blank(),
      axis.title.x=element_text(size=44, face='bold', colour="black"), #x-axis title font
      axis.ticks.x=element_line(linewidth=1),
      axis.text.y=element_text(size=36, face='bold'), #y-axis font
      axis.title.y=element_text(size=44, face='bold', colour="black"), #y-axis title font
      axis.ticks.y=element_line(linewidth=1),
      panel.grid.major.x=element_blank(),
      panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.75,  linetype=3),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_line(color="#8ccde3",  linewidth=0.75,  linetype=3),
      panel.background=element_rect(fill=NA),
      legend.position="none",
      legend.key=element_rect(fill = 'white'),
      legend.title=element_text(size=44, face = "bold"),
      legend.text=element_text(size=40),
      plot.title=element_text(size=44, face='bold', color="darkred"),
      plot.margin=margin(t=5, r=5, b=5, l=35, unit="pt"),
      strip.text=element_text(size=40, face='bold'))+
    facet_nested(no_sources ~ factor(fit_criteria1, level=c("% (-0.3,1.3)", "% (0,1)", "%+/-2SE (0,1)")) + factor(qpAdm_setup, levels=c("dist. rot.", "prox. rot.", "dist. non-rot.", "prox. non-rot.")) + no_outgroups_proxies)
  
  p3<-ggplot(df2, aes(fill=metric, x=factor(fit_criteria2, level=c("0.001", "0.05", "0.5", "0.001,   0.5")), y=value))+ 
    geom_bar(position="stack", stat="identity", width=.8)+
    scale_fill_manual(values=lab, breaks=leg) + #provide a user-defined color palette here
    scale_y_continuous(expand=c(0.001,0.001), breaks=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1), limits=c(0,1))+
    xlab('p-value thresholds')+
    ylab('fraction of all experiments')+
    guides(fill = guide_legend(title = 'at least one (1, 2, 3, 4-way) model for the target is fitting and:', order=1, keywidth=4, keyheight=4, override.aes = list(size = 6)))+
    theme(
      panel.ontop=T,
      axis.text.x=element_text(size=36, face='bold', angle=45, hjust=0.95, vjust=0.95),
      axis.title.x=element_text(size=48, face='bold', colour="black"), #x-axis title font
      axis.ticks.x=element_line(linewidth=1),
      axis.text.y=element_text(size=36, face='bold'), #y-axis font
      axis.title.y=element_text(size=44, face='bold', colour="black"), #y-axis title font
      axis.ticks.y=element_line(linewidth=1),
      panel.grid.major.y=element_line(color="#8ccde3",  linewidth=0.75,  linetype=3),
      panel.grid.major.x=element_blank(),
      panel.grid.minor=element_blank(),
      panel.background=element_rect(fill=NA),
      legend.position="right",
      legend.key=element_rect(fill = 'white'),
      legend.title=element_text(size=48, face = "bold"),
      legend.text=element_text(size=44),
      plot.title=element_text(size=44, face='bold', color="darkred"),
      plot.margin=margin(t=5, r=52, b=5, l=35, unit="pt"),
      strip.text=element_text(size=40, face='bold'))+
    facet_nested(. ~ factor(fit_criteria1, level=c("% (-0.3,1.3)", "% (0,1)", "%+/-2SE (0,1)")) + factor(qpAdm_setup, levels=c("dist. rot.", "prox. rot.", "dist. non-rot.", "prox. non-rot.")) + no_outgroups_proxies)
  
  mylegend<-g_legend(p3)
  
  pall <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                   p2 + theme(legend.position="none"),
                                   p3 + theme(legend.position="none"),
                                   nrow=3, heights=c(1.5, 1.5, 2)),
                                   mylegend, nrow=2, heights=c(5, 0.8))
  
  ###The plot for the "10-3 to 10-2" SSL is Figure 8###
  pdf(paste0("experiments_stacked_histogram_dist_to_corner_maxdist_diffsqrt8_errorbars3se_v5_p2_",gf,".pdf"), width=48, height=66)
  print(annotate_figure(pall, top=text_grob(paste0('qpAdm performance metrics at the level of experiments (based on Euclidean distance to the ideal symmetric model);\n13 or 18 proxy and "right" groups randomly picked in each experiment; 400 experiments per 10 sim. iters., protocol, and sampling density;\n"trailing" simpler models were not checked as part of the feasibility criteria; gene-flow intensities: ',gf,'\n'), color="darkmagenta", face="bold", size=50)))
  dev.off()
}
