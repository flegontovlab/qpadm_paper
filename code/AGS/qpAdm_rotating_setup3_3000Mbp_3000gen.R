library(tidyverse)
library(admixtools)

args = commandArgs(trailingOnly = TRUE)
sim = args[1]
subit = args[2]

basedir = "./projects/qpadm/results/f2/30chr_deep"
outdir = "./projects/qpadm/results/qpadm/30chr_deep"

f=file.path(basedir, sim, subit)
f2_blocks = f2_from_precomp(f)
pops = rownames(f2_blocks[,,1])

qpadm_results_1 = list()
qpadm_results_2 = list()
for (target in pops){
  remaining_pops = setdiff(pops, target)
  all_left_pairs_2 = t(combn(remaining_pops, 2))
  all_left_pairs_1 = t(combn(remaining_pops, 1))
  
  qpadm_out_2 = apply(all_left_pairs_2, 1, function(x) {
    right_pops = setdiff(pops, c(target, x))
    
    out = qpadm(data = f2_blocks, 
                target = target, left = x, 
                right=right_pops, 
                boot = F, getcov = T, cpp = T, verbose = T)
    nm = paste0(target, "; ", paste(x, collapse=", "), "; ", paste(right_pops, collapse=", "))
    mylist = list()
    mylist[[nm]] = out
    mylist
  }) %>% unlist(recursive = F)
  qpadm_results_2 = append(qpadm_results_2, qpadm_out_2)
  
  qpadm_out_1 = apply(all_left_pairs_1, 1, function(x) {
    right_pops = setdiff(pops, c(target, x))
    
    out = qpadm(data = f2_blocks, 
                target = target, left = x, 
                right=right_pops, 
                boot = F, getcov = T, cpp = T, verbose = T)
    nm = paste0(target, "; ", paste(x, collapse=", "), "; ", paste(right_pops, collapse=", "))
    mylist = list()
    mylist[[nm]] = out
    mylist
  }) %>% unlist(recursive = F)
  qpadm_results_1 = append(qpadm_results_1, qpadm_out_1)
}

sf1 = file.path(outdir, sim, paste0(subit, "_1way_rotating.rds"))
sf2 = file.path(outdir, sim, paste0(subit, "_2way_rotating.rds"))
saveRDS(qpadm_results_1, file=sf1)
saveRDS(qpadm_results_2, file=sf2)
