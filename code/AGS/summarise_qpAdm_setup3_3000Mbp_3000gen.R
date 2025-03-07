library(tidyverse)
library(admixtools)
library(magrittr)

args = commandArgs(trailingOnly = TRUE)
sim = args[1]

outs = "./projects/qpadm/results/qpadm/30chr_deep"

outdir = file.path(outs, sim)

for (replc in 1:10){
  out1 = readRDS(paste0(outdir, '/qpadm_sim_', replc, "_1way_rotating.rds"))
  out2 = readRDS(paste0(outdir, '/qpadm_sim_', replc, "_2way_rotating.rds"))
  
  combs = names(out2)
  
  outtab = sapply(combs, function(comb){
    target = strsplit(comb, ";")[[1]][1]
    s = strsplit(comb, ";")[[1]][2]
    s1 = strsplit(gsub(" ", "", s), ",")[[1]][1]
    s2 = strsplit(gsub(" ", "", s), ",")[[1]][2]
    
    out_s1 = names(out1)[startsWith(names(out1), paste0(target, "; ", s1, "; "))]
    out_s2 = names(out1)[startsWith(names(out1), paste0(target, "; ", s2, "; "))]
    out_s1_s2 = names(out1)[startsWith(names(out1), paste0(s1, "; ", s2, "; "))]
    
    weight_s1 = out2[[comb]]$weights %>%
      filter(target == target, left == s1) %>%
      pull(weight)
    
    se_s1 = out2[[comb]]$weights %>%
      filter(target == target, left == s1) %>%
      pull(se)
    
    weight_s2 = out2[[comb]]$weights %>%
      filter(target == target, left == s2) %>%
      pull(weight)
    
    se_s2 = out2[[comb]]$weights %>%
      filter(target == target, left == s2) %>%
      pull(se)
    
    c("target" = target, "s1" = s1, "s2" = s2, 
      "weight_t.s1" = weight_s1, "se_t.s1" = se_s1, 
      "weight_t.s2" = weight_s2, "se_t.s2" = se_s2,  
      "p_2way" = out2[[comb]]$rankdrop$p[1],
      "p_1way_t.s1" = out1[[out_s1]]$rankdrop$p, 
      "p_1way_t.s2" = out1[[out_s2]]$rankdrop$p,
      "p_1way_s1.s2" = out1[[out_s1_s2]]$rankdrop$p)
  }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
  
  write.csv(outtab, file=paste0(outdir, "/qpadm_summary_rotating_", sim, "_", replc, ".csv"), row.names=FALSE)
}

