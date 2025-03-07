library(tidyverse)
library(admixtools)

args = commandArgs(trailingOnly = TRUE)
sim = args[1]
subit = args[2]

basedir = "./projects/qpadm/results/f2/30chr_deep"
outdir = "./projects/qpadm/results/qpadm/30chr_deep"

right_groups = list(
  sim_1 = c("A","B","C","E","F","L"),
  sim_2 = c("B","C","E","G","H","I"),
  sim_3 = c("D","E","G","I","J","M"),
  sim_4 = c("A","B","D","E","L","M"),
  sim_5 = c("A","B","D","H","J","L"),
  sim_6 = c("B","E","I","J","K","M"),
  sim_7 = c("A","B","E","F","J","M"),
  sim_8 = c("A","C","E","H","J","L"),
  sim_9 = c("A","B","D","E","F","K"),
  sim_10 = c("A","C","F","H","I","K"),
  sim_11 = c("A","B","G","H","I","J"),
  sim_12 = c("B","D","F","G","H","J"),
  sim_13 = c("A","C","H","I","J","L"),
  sim_14 = c("A","D","F","H","I","L"),
  sim_15 = c("B","D","F","H","K","M"),
  sim_16 = c("C","D","F","G","I","L"),
  sim_17 = c("B","C","D","I","J","L"),
  sim_18 = c("B","E","F","G","J","L"),
  sim_19 = c("E","F","G","K","L","M"),
  sim_20 = c("A","D","F","G","J","L"),
  sim_21 = c("A","B","C","G","K","M"),
  sim_22 = c("B","E","F","G","I","K"),
  sim_23 = c("A","B","E","G","J","K"),
  sim_24 = c("A","E","F","H","J","K"),
  sim_25 = c("C","D","F","G","J","M"),
  sim_26 = c("B","C","D","H","I","M"),
  sim_27 = c("B","C","D","F","I","K"),
  sim_28 = c("C","D","E","F","H","L"),
  sim_29 = c("C","E","I","J","K","M"),
  sim_30 = c("A","C","E","H","K","L"),
  sim_31 = c("B","C","E","G","J","M"),
  sim_32 = c("E","H","I","J","K","M"),
  sim_33 = c("A","E","G","J","L","M"),
  sim_34 = c("A","B","C","E","H","I"),
  sim_35 = c("A","D","I","J","K","M"),
  sim_36 = c("B","D","E","F","I","J"),
  sim_37 = c("A","D","E","G","H","L"),
  sim_38 = c("A","D","F","G","H","K"),
  sim_39 = c("A","C","D","G","I","J"),
  sim_40 = c("A","B","D","F","I","M")
)

right_pops = right_groups[[sim]]
f = file.path(basedir, sim, subit)

f2_blocks = f2_from_precomp(f)
allpops = rownames(f2_blocks[,,1])

pops = setdiff(allpops, right_pops)

qpadm_results_1 = list()
qpadm_results_2 = list()
for (target in pops){
  remaining_pops = setdiff(pops, target)
  all_left_pairs_2 = t(combn(remaining_pops, 2))
  all_left_pairs_1 = t(combn(remaining_pops, 1))
  
  qpadm_out_2 = apply(all_left_pairs_2, 1, function(x) {
    # right_pops = setdiff(pops, c(target, x))
    
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
    # right_pops = setdiff(pops, c(target, x))
    
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

sf1 = file.path(outdir, sim, paste0(subit, "_1way_nonrotating.rds"))
sf2 = file.path(outdir, sim, paste0(subit, "_2way_nonrotating.rds"))
saveRDS(qpadm_results_1, file=sf1)
saveRDS(qpadm_results_2, file=sf2)
