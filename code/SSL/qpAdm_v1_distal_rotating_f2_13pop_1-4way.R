#devtools::install_github("uqrmaie1/admixtools")
library(tidyverse)
library(readr)
library(admixtools)
library(future)
library(future.apply)
library(magrittr)

setwd('./stepping-stone_simulations_4way')

###randomized rotating distal qpAdm###
future::plan(multicore, workers=10)

args = commandArgs(trailingOnly=TRUE)

poplist <- read.table(paste0("./stepping-stone_simulations_4way/",args[1],".ind"), header=F, sep="\t")$V3 %>% unique()
rotating <- grep("_present", poplist, value=T, invert=T)
targets <- grep("_present", poplist, value=T)
selected_targets <- sample(targets, 20, replace=F) %>% sort()

for (t in selected_targets){
  iters <- c(1:2)
  for (i in iters){
    selected_rotating <- sample(rotating, 13, replace=F) %>% sort()
    proxies1 = t(combn(selected_rotating, 1))
    proxies2 = t(combn(selected_rotating, 2))
    proxies3 = t(combn(selected_rotating, 3))
    proxies4 = t(combn(selected_rotating, 4))
    extract_f2(pref=paste0("./stepping-stone_simulations_4way/",args[1]), outdir=paste0(args[1],"_distal_rotating_13pop_",t,"_",i), 
               pops=c(t, selected_rotating), auto_only=F, maxmiss=0, minac2=F, 
               blgsize=4000000, adjust_pseudohaploid=F, maxmem = 50000, overwrite=T, verbose=T)
    f2_blocks = f2_from_precomp(paste0(args[1],"_distal_rotating_13pop_",t,"_",i))
    
    ###1-way###
    test1a = future_apply(proxies1, 1, function(x) {
      out = qpadm(data = f2_blocks, 
              target = t, left = x, 
              right = setdiff(selected_rotating, x), 
              fudge_twice=T)
      nm = paste0(t, "; ", paste(x, collapse=", "), "; ", paste(setdiff(selected_rotating, x), collapse=", "))
      mylist = list()
      mylist[[nm]] = out
      mylist
    }) %>% unlist(recursive = F)

    test1b = future_apply(proxies2, 1, function(x) {
      out = qpadm(data = f2_blocks, 
              target = x[1], left = x[2], 
              right = setdiff(selected_rotating, x), 
              fudge_twice=T)
      nm = paste0(x[1], "; ", paste(x[2], collapse=", "), "; ", paste(setdiff(selected_rotating, x), collapse=", "))
      mylist = list()
      mylist[[nm]] = out
      mylist
    }) %>% unlist(recursive = F)

    ###2-way###
    test2a = future_apply(proxies2, 1, function(x) {
      out = qpadm(data = f2_blocks, 
              target = t, left = x, 
              right = setdiff(selected_rotating, x), 
              fudge_twice=T)
      nm = paste0(t, "; ", paste(x, collapse=", "), "; ", paste(setdiff(selected_rotating, x), collapse=", "))
      mylist = list()
      mylist[[nm]] = out
      mylist
      }) %>% unlist(recursive = F)
    
    test2b = future_apply(proxies3, 1, function(x) {
      out = qpadm(data = f2_blocks, 
              target = x[1], left = c(x[2],x[3]), 
              right = setdiff(selected_rotating, x), 
              fudge_twice=T)
      nm = paste0(paste(x, collapse=", "), "; ", paste(setdiff(selected_rotating, x), collapse=", "))
      mylist = list()
      mylist[[nm]] = out
      mylist
    }) %>% unlist(recursive = F)
    
    ###3-way###
    test3a = future_apply(proxies3, 1, function(x) {
      out = qpadm(data = f2_blocks, 
              target = t, left = x, 
              right = setdiff(selected_rotating, x), 
              fudge_twice=T)
      nm = paste0(t, "; ", paste(x, collapse=", "), "; ", paste(setdiff(selected_rotating, x), collapse=", "))
      mylist = list()
      mylist[[nm]] = out
      mylist
    }) %>% unlist(recursive = F)
    
    test3b = future_apply(proxies4, 1, function(x) {
      out = qpadm(data = f2_blocks, 
                  target = x[1], left = c(x[2],x[3],x[4]), 
                  right = setdiff(selected_rotating, x), 
                  fudge_twice=T)
      nm = paste0(paste(x, collapse=", "), "; ", paste(setdiff(selected_rotating, x), collapse=", "))
      mylist = list()
      mylist[[nm]] = out
      mylist
    }) %>% unlist(recursive = F)
    
    ###4-way###
    test4 = future_apply(proxies4, 1, function(x) {
      out = qpadm(data = f2_blocks, 
                  target = t, left = x, 
                  right = setdiff(selected_rotating, x), 
                  fudge_twice=T)
      nm = paste0(t, "; ", paste(x, collapse=", "), "; ", paste(setdiff(selected_rotating, x), collapse=", "))
      mylist = list()
      mylist[[nm]] = out
      mylist
    }) %>% unlist(recursive = F)

    ###summarizing 1-way results###
    combs = names(test1a)
    outtab = sapply(combs, function(comb){
      target = strsplit(comb, ";")[[1]][1]
      s = strsplit(comb, ";")[[1]][2]
      c("target" = target, "s" = s, 
        "p_1way" = test1a[[comb]]$rankdrop$p[1],
        "rotating" = paste(selected_rotating, collapse=", "))
    }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
    write.table(outtab, file=paste0(args[1],"_distal_rotating_13pop_",t,"_",i,"_1way.txt"), sep='\t', row.names=F)
    
    ###summarizing 2-way results###
    combs = names(test2a)
    outtab = sapply(combs, function(comb){
      target = strsplit(comb, ";")[[1]][1]
      s = strsplit(comb, ";")[[1]][2]
      s1 = strsplit(gsub(" ", "", s), ",")[[1]][1]
      s2 = strsplit(gsub(" ", "", s), ",")[[1]][2]
      out_t_s1 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s1, "; "))]
      out_t_s2 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s2, "; "))]
      out_s1_s2 = names(test1b)[startsWith(names(test1b), paste0(s1, "; ", s2, "; "))]
      weight_s1 = test2a[[comb]]$weights %>% filter(target == target, left == s1) %>% pull(weight)
      se_s1 = test2a[[comb]]$weights %>% filter(target == target, left == s1) %>% pull(se)
      weight_s2 = test2a[[comb]]$weights %>% filter(target == target, left == s2) %>% pull(weight)
      se_s2 = test2a[[comb]]$weights %>% filter(target == target, left == s2) %>% pull(se)
      c("target" = target, "s1" = s1, "s2" = s2, 
        "weight.s1" = weight_s1, "se.s1" = se_s1, 
        "weight.s2" = weight_s2, "se.s2" = se_s2,  
        "p_2way" = test2a[[comb]]$rankdrop$p[1],
        "p_1way_t.s1" = test1a[[out_t_s1]]$rankdrop$p, 
        "p_1way_t.s2" = test1a[[out_t_s2]]$rankdrop$p,
        "p_1way_s1.s2" = test1b[[out_s1_s2]]$rankdrop$p,
        "rotating" = paste(selected_rotating, collapse=", "))
    }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
    write.table(outtab, file=paste0(args[1],"_distal_rotating_13pop_",t,"_",i,"_2way.txt"), sep='\t', row.names=F)
    
    ###summarizing 3-way results###
    combs = names(test3a)
    outtab = sapply(combs, function(comb){
      target = strsplit(comb, ";")[[1]][1]
      s = strsplit(comb, ";")[[1]][2]
      s1 = strsplit(gsub(" ", "", s), ",")[[1]][1]
      s2 = strsplit(gsub(" ", "", s), ",")[[1]][2]
      s3 = strsplit(gsub(" ", "", s), ",")[[1]][3]
      out_t_s1 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s1, "; "))]
      out_t_s2 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s2, "; "))]
      out_t_s3 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s3, "; "))]
      out_s1_s2 = names(test1b)[startsWith(names(test1b), paste0(s1, "; ", s2, "; "))]
      out_s1_s3 = names(test1b)[startsWith(names(test1b), paste0(s1, "; ", s3, "; "))]
      out_s2_s3 = names(test1b)[startsWith(names(test1b), paste0(s2, "; ", s3, "; "))]
      out_t_s1_s2 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s1, ", ", s2, "; "))] #pay attention to the punctuation
      out_t_s1_s3 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s1, ", ", s3, "; "))] #pay attention to the punctuation
      out_t_s2_s3 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s2, ", ", s3, "; "))] #pay attention to the punctuation
      out_s1_s2_s3 = names(test2b)[startsWith(names(test2b), paste0(s1, ", ", s2, ", ", s3, "; "))] #pay attention to the punctuation
      weight_s1 = test3a[[comb]]$weights %>% filter(target == target, left == s1) %>% pull(weight)
      se_s1 = test3a[[comb]]$weights %>% filter(target == target, left == s1) %>% pull(se)
      weight_s2 = test3a[[comb]]$weights %>% filter(target == target, left == s2) %>% pull(weight)
      se_s2 = test3a[[comb]]$weights %>% filter(target == target, left == s2) %>% pull(se)
      weight_s3 = test3a[[comb]]$weights %>% filter(target == target, left == s3) %>% pull(weight)
      se_s3 = test3a[[comb]]$weights %>% filter(target == target, left == s3) %>% pull(se)
      c("target" = target, "s1" = s1, "s2" = s2, "s3" = s3, 
        "weight.s1" = weight_s1, "se.s1" = se_s1, 
        "weight.s2" = weight_s2, "se.s2" = se_s2,  
        "weight.s3" = weight_s3, "se.s3" = se_s3,  
        "p_3way" = test3a[[comb]]$rankdrop$p[1],
        "p_1way_t.s1" = test1a[[out_t_s1]]$rankdrop$p, 
        "p_1way_t.s2" = test1a[[out_t_s2]]$rankdrop$p,
        "p_1way_t.s3" = test1a[[out_t_s3]]$rankdrop$p,
        "p_1way_s1.s2" = test1b[[out_s1_s2]]$rankdrop$p,
        "p_1way_s1.s3" = test1b[[out_s1_s3]]$rankdrop$p,
        "p_1way_s2.s3" = test1b[[out_s2_s3]]$rankdrop$p,
        "p_2way_t.s1.s2" = test2a[[out_t_s1_s2]]$rankdrop$p[1], 
        "p_2way_t.s1.s3" = test2a[[out_t_s1_s3]]$rankdrop$p[1],
        "p_2way_t.s2.s3" = test2a[[out_t_s2_s3]]$rankdrop$p[1],
        "p_2way_s1.s2.s3" = test2b[[out_s1_s2_s3]]$rankdrop$p[1],
        "rotating" = paste(selected_rotating, collapse=", "))
    }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
    write.table(outtab, file=paste0(args[1],"_distal_rotating_13pop_",t,"_",i,"_3way.txt"), sep='\t', row.names=F)
    
    ###summarizing 4-way results###
    combs = names(test4)
    outtab = sapply(combs, function(comb){
      target = strsplit(comb, ";")[[1]][1]
      s = strsplit(comb, ";")[[1]][2]
      s1 = strsplit(gsub(" ", "", s), ",")[[1]][1]
      s2 = strsplit(gsub(" ", "", s), ",")[[1]][2]
      s3 = strsplit(gsub(" ", "", s), ",")[[1]][3]
      s4 = strsplit(gsub(" ", "", s), ",")[[1]][4]
      out_t_s1 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s1, "; "))]
      out_t_s2 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s2, "; "))]
      out_t_s3 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s3, "; "))]
      out_t_s4 = names(test1a)[startsWith(names(test1a), paste0(target, "; ", s4, "; "))]
      out_s1_s2 = names(test1b)[startsWith(names(test1b), paste0(s1, "; ", s2, "; "))]
      out_s1_s3 = names(test1b)[startsWith(names(test1b), paste0(s1, "; ", s3, "; "))]
      out_s1_s4 = names(test1b)[startsWith(names(test1b), paste0(s1, "; ", s4, "; "))]
      out_s2_s3 = names(test1b)[startsWith(names(test1b), paste0(s2, "; ", s3, "; "))]
      out_s2_s4 = names(test1b)[startsWith(names(test1b), paste0(s2, "; ", s4, "; "))]
      out_s3_s4 = names(test1b)[startsWith(names(test1b), paste0(s3, "; ", s4, "; "))]
      out_t_s1_s2 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s1, ", ", s2, "; "))] #pay attention to the punctuation
      out_t_s1_s3 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s1, ", ", s3, "; "))] #pay attention to the punctuation
      out_t_s1_s4 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s1, ", ", s4, "; "))] #pay attention to the punctuation
      out_t_s2_s3 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s2, ", ", s3, "; "))] #pay attention to the punctuation
      out_t_s2_s4 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s2, ", ", s4, "; "))] #pay attention to the punctuation
      out_t_s3_s4 = names(test2a)[startsWith(names(test2a), paste0(target, "; ", s3, ", ", s4, "; "))] #pay attention to the punctuation
      out_s1_s2_s3 = names(test2b)[startsWith(names(test2b), paste0(s1, ", ", s2, ", ", s3, "; "))] #pay attention to the punctuation
      out_s1_s2_s4 = names(test2b)[startsWith(names(test2b), paste0(s1, ", ", s2, ", ", s4, "; "))] #pay attention to the punctuation
      out_s1_s3_s4 = names(test2b)[startsWith(names(test2b), paste0(s1, ", ", s3, ", ", s4, "; "))] #pay attention to the punctuation
      out_s2_s3_s4 = names(test2b)[startsWith(names(test2b), paste0(s2, ", ", s3, ", ", s4, "; "))] #pay attention to the punctuation
      out_t_s1_s2_s3 = names(test3a)[startsWith(names(test3a), paste0(target, "; ", s1, ", ", s2, ", ", s3, "; "))] #pay attention to the punctuation
      out_t_s1_s2_s4 = names(test3a)[startsWith(names(test3a), paste0(target, "; ", s1, ", ", s2, ", ", s4, "; "))] #pay attention to the punctuation
      out_t_s1_s3_s4 = names(test3a)[startsWith(names(test3a), paste0(target, "; ", s1, ", ", s3, ", ", s4, "; "))] #pay attention to the punctuation
      out_t_s2_s3_s4 = names(test3a)[startsWith(names(test3a), paste0(target, "; ", s2, ", ", s3, ", ", s4, "; "))] #pay attention to the punctuation
      out_s1_s2_s3_s4 = names(test3b)[startsWith(names(test3b), paste0(s1, ", ", s2, ", ", s3, ", ", s4, "; "))] #pay attention to the punctuation
      weight_s1 = test4[[comb]]$weights %>% filter(target == target, left == s1) %>% pull(weight)
      se_s1 = test4[[comb]]$weights %>% filter(target == target, left == s1) %>% pull(se)
      weight_s2 = test4[[comb]]$weights %>% filter(target == target, left == s2) %>% pull(weight)
      se_s2 = test4[[comb]]$weights %>% filter(target == target, left == s2) %>% pull(se)
      weight_s3 = test4[[comb]]$weights %>% filter(target == target, left == s3) %>% pull(weight)
      se_s3 = test4[[comb]]$weights %>% filter(target == target, left == s3) %>% pull(se)
      weight_s4 = test4[[comb]]$weights %>% filter(target == target, left == s4) %>% pull(weight)
      se_s4 = test4[[comb]]$weights %>% filter(target == target, left == s4) %>% pull(se)
      c("target" = target, "s1" = s1, "s2" = s2, "s3" = s3, "s4" = s4,
        "weight.s1" = weight_s1, "se.s1" = se_s1, 
        "weight.s2" = weight_s2, "se.s2" = se_s2,  
        "weight.s3" = weight_s3, "se.s3" = se_s3,  
        "weight.s4" = weight_s4, "se.s4" = se_s4,  
        "p_4way" = test4[[comb]]$rankdrop$p[1],
        "p_1way_t.s1" = test1a[[out_t_s1]]$rankdrop$p, 
        "p_1way_t.s2" = test1a[[out_t_s2]]$rankdrop$p,
        "p_1way_t.s3" = test1a[[out_t_s3]]$rankdrop$p,
        "p_1way_t.s4" = test1a[[out_t_s4]]$rankdrop$p,
        "p_1way_s1.s2" = test1b[[out_s1_s2]]$rankdrop$p,
        "p_1way_s1.s3" = test1b[[out_s1_s3]]$rankdrop$p,
        "p_1way_s1.s4" = test1b[[out_s1_s4]]$rankdrop$p,
        "p_1way_s2.s3" = test1b[[out_s2_s3]]$rankdrop$p,
        "p_1way_s2.s4" = test1b[[out_s2_s4]]$rankdrop$p,
        "p_1way_s3.s4" = test1b[[out_s3_s4]]$rankdrop$p,
        "p_2way_t.s1.s2" = test2a[[out_t_s1_s2]]$rankdrop$p[1], 
        "p_2way_t.s1.s3" = test2a[[out_t_s1_s3]]$rankdrop$p[1],
        "p_2way_t.s1.s4" = test2a[[out_t_s1_s4]]$rankdrop$p[1],
        "p_2way_t.s2.s3" = test2a[[out_t_s2_s3]]$rankdrop$p[1],
        "p_2way_t.s2.s4" = test2a[[out_t_s2_s4]]$rankdrop$p[1],
        "p_2way_t.s3.s4" = test2a[[out_t_s3_s4]]$rankdrop$p[1],
        "p_2way_s1.s2.s3" = test2b[[out_s1_s2_s3]]$rankdrop$p[1],
        "p_2way_s1.s2.s4" = test2b[[out_s1_s2_s4]]$rankdrop$p[1],
        "p_2way_s1.s3.s4" = test2b[[out_s1_s3_s4]]$rankdrop$p[1],
        "p_2way_s2.s3.s4" = test2b[[out_s2_s3_s4]]$rankdrop$p[1],
        "p_3way_t.s1.s2.s3" = test3a[[out_t_s1_s2_s3]]$rankdrop$p[1],
        "p_3way_t.s1.s2.s4" = test3a[[out_t_s1_s2_s4]]$rankdrop$p[1],
        "p_3way_t.s1.s3.s4" = test3a[[out_t_s1_s3_s4]]$rankdrop$p[1],
        "p_3way_t.s2.s3.s4" = test3a[[out_t_s2_s3_s4]]$rankdrop$p[1],
        "p_3way_s1.s2.s3.s4" = test3b[[out_s1_s2_s3_s4]]$rankdrop$p[1],
        "rotating" = paste(selected_rotating, collapse=", "))
    }) %>% t() %>% as.data.frame() %>% set_rownames(NULL)
    write.table(outtab, file=paste0(args[1],"_distal_rotating_13pop_",t,"_",i,"_4way.txt"), sep='\t', row.names=F)
  }
}

