library(admixtools)
library(tidyverse)

nleaves = 13
nadmix = 10
gen_time = 25
max_depth = 20000 / gen_time
seq_length = 1e8
nchr = 3
neff_range = c(1000, 10000)
admix_range = c(0.1, 0.5)
time_range = c(500, 3000) / gen_time
samples = 10
mu = 1.25e-8
recomb_rate = 2e-8

num = 100
for (n in 21:num){
  pref = paste0("./projects/qpadm/results/simout/qpadm_sim_", n, ".py")
  
  out = random_sim(nleaf = nleaves, nadmix = nadmix, outpref = pref, max_depth = max_depth, ind_per_pop = samples,
                   mutation_rate = mu, admix_weights = admix_range, neff = neff_range, time = time_range, 
                   nchr = nchr, seq_length = seq_length, recomb_rate = recomb_rate, run = FALSE)
  
  saveRDS(out, file=paste0("./projects/qpadm/results/simmeta/qpadm_sim_", n, ".rds"))
  
  p = plot_graph(out$edges, dates=out$dates, hide_weights=T, textsize = 3, scale_y = FALSE,
                 title=paste0("qpadm_sim_", n), fix = T)
  ggsave(paste0("./projects/qpadm/results/topos/qpadm_sim_", n, ".png"), p, width = 12, height = 8)
  
}

