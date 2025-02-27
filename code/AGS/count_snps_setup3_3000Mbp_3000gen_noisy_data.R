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

inpref = "./projects/qpadm/results/simout/30chr_deep_subsampled"
outf = "./projects/qpadm/results/snp_counts/snp_counts_30chr_deep_subsampled.rds"

mylist = apply(as.matrix(expand.grid(1:40, 1:10)), 1, function(x){
      sim = x[1]
      subit = x[2]
      inpf = paste0(inpref, '/sim_', sim,  '/qpadm_sim_1_sub', subit)
      n = num_polymorphic(inpf)
      c(paste0('sim_', sim), paste0('subit_', subit), n)
					     })

mydf = mylist %>%
	t() %>%
	as.data.frame() %>%
	magrittr::set_colnames(c('sim', 'subit', 'count'))

saveRDS(mydf, outf)



