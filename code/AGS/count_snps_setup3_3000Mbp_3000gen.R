library(admixtools)
library(tidyverse)
library(reshape2)

num_polymorphic = function(f){
  ## Calculate number of polymorphic sites for given simulation
  afdat = admixtools:::anygeno_to_aftable(f, format = NULL,
                                          adjust_pseudohaploid = TRUE, verbose = FALSE)
  
  afdat2 = admixtools:::discard_from_aftable(afdat, maxmiss = 0, minmaf = 0, maxmaf = 0.5, minac2 = FALSE,
                                             outpop = NULL,
                                             transitions = TRUE, transversions = TRUE,
                                             keepsnps = NULL, auto_only = TRUE, poly_only = FALSE)
  
  snpfile = afdat2$snpfile %>% mutate(poly = as.logical(admixtools:::cpp_is_polymorphic(afdat2$afs)))
  sum(snpfile$poly)
}

inpref = "./projects/qpadm/results/simout/30chr_deep"
outf = "./projects/qpadm/results/snp_counts/snp_counts_30chr_deep_all.rds"

mylist = list()
for (sim in 1:40){
  for (subit in 1:10){
    inpf = paste0(inpref, '/sim_', sim, '/qpadm_sim_', subit)
    c = num_polymorphic(inpf)
    
    mylist[[paste0('sim_', sim)]][[paste0('subit_', subit)]] = c
  }
}

mydf = melt(mylist, measure.vars=c()) %>% 
  transmute(sim_no = L1, subit = L2, count = value)

saveRDS(mydf, outf)
