library(tidyverse)
library(admixtools)

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

run_subsample = function(inp, out, subsample_script, python_){
  ## Run subsampling script until a subsampled replicate with number of polymorphic sites greater than 20k is obtained
  nump = 0
  while (nump < 66000){
    cmd = paste0(python_, " ", subsample_script, " --in ", inp, " --out ", out)
    system(cmd)
    
    nump = num_polymorphic(out)
    if (nump < 66000){
      print(paste0("Num of polymorphic sites = ", nump, ". Removing files.."))
      print("Retrying..")
      system(paste0("rm ", out, ".geno"))
      system(paste0("rm ", out, ".snp"))
      system(paste0("rm ", out, ".ind"))
    } else {
      print(paste0("Found a subsampled iteration with ", nump, " polymorphic site."))
    }
  }
}


# inp_pref = "./Projects/qpadm/data/simout/qpadm_sim_"
# out_pref = "./Projects/qpadm/data/simout/subsampled/qpadm_sim_"
# subsample_script = "./Projects/qpadm/scripts/generate_noisy_data.py"
# python = "./opt/anaconda3/bin/python"

inp_pref = "./projects/qpadm/results/simout/30chr_deep/sim_"
out_pref = "./projects/qpadm/results/simout/30chr_deep_subsampled/sim_"
subsample_script = "./projects/qpadm/scripts/generate_noisy_data.py"
python = "./miniconda3/bin/python"
sims = 1:40
sub_iters = 1:10

for (sim in sims){
  inp_file = paste0(inp_pref, sim, '/qpadm_sim_1')
  out_file = paste0(out_pref, sim)
  for (sub_iter in sub_iters){
    out_file_2 = paste0(out_file, "/qpadm_sim_1_sub", sub_iter)
    print(out_file_2)
    run_subsample(inp = inp_file, out = out_file_2, subsample_script = subsample_script, python_ = python)
  }
}

