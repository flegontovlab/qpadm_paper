import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("K", "A", "E", "J", "H", "B", "G", "C", "L", "D", "F", "I", "M"), (466, 200, 642, 616, 616, 33, 98, 0, 56, 0, 466, 56, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 1554, name = "R")
demo.add_population(initial_size = 3358, name = "Rl")
demo.add_population(initial_size = 3962, name = "Rll")
demo.add_population(initial_size = 6605, name = "Rlll")
demo.add_population(initial_size = 2441, name = "Rllll")
demo.add_population(initial_size = 6043, name = "K")
demo.add_population(initial_size = 6872, name = "Rlllr")
demo.add_population(initial_size = 6880, name = "Rlllrl")
demo.add_population(initial_size = 5316, name = "Rlllrlr")
demo.add_population(initial_size = 2466, name = "A")
demo.add_population(initial_size = 6087, name = "Rllr")
demo.add_population(initial_size = 5395, name = "E")
demo.add_population(initial_size = 2432, name = "Rr")
demo.add_population(initial_size = 8376, name = "Rrl")
demo.add_population(initial_size = 9735, name = "J")
demo.add_population(initial_size = 2949, name = "H")
demo.add_population(initial_size = 6829, name = "Rllllr")
demo.add_population(initial_size = 8517, name = "Rllllrd")
demo.add_population(initial_size = 5178, name = "admix")
demo.add_population(initial_size = 4295, name = "B")
demo.add_population(initial_size = 9789, name = "Rlllrlry")
demo.add_population(initial_size = 2394, name = "G")
demo.add_population(initial_size = 3261, name = "admixw")
demo.add_population(initial_size = 4276, name = "C")
demo.add_population(initial_size = 4702, name = "Rllra")
demo.add_population(initial_size = 5547, name = "admixi")
demo.add_population(initial_size = 1592, name = "Rlllls")
demo.add_population(initial_size = 8660, name = "admixy")
demo.add_population(initial_size = 1680, name = "Rllllrs")
demo.add_population(initial_size = 4553, name = "L")
demo.add_population(initial_size = 2525, name = "admixk")
demo.add_population(initial_size = 9995, name = "D")
demo.add_population(initial_size = 5159, name = "admixiu")
demo.add_population(initial_size = 1720, name = "admixd")
demo.add_population(initial_size = 8354, name = "Rrp")
demo.add_population(initial_size = 1513, name = "admixu")
demo.add_population(initial_size = 4575, name = "F")
demo.add_population(initial_size = 4878, name = "Rllllrf")
demo.add_population(initial_size = 1849, name = "admixwh")
demo.add_population(initial_size = 2011, name = "I")
demo.add_population(initial_size = 5296, name = "Rlllrlrm")
demo.add_population(initial_size = 7809, name = "admixe")
demo.add_population(initial_size = 6195, name = "Rllllrsl")
demo.add_population(initial_size = 2136, name = "admixz")
demo.add_population(initial_size = 2842, name = "M")
demo.add_population_split(time=56, ancestral="admix", derived=["B"])
demo.add_population_split(time=363, ancestral="admixi", derived=["admixiu"])
demo.add_population_split(time=279, ancestral="admixiu", derived=["Rllllr"])
demo.add_population_split(time=33, ancestral="admixk", derived=["D"])
demo.add_population_split(time=570, ancestral="admixu", derived=["F"])
demo.add_population_split(time=33, ancestral="admixw", derived=["C"])
demo.add_population_split(time=98, ancestral="admixwh", derived=["I"])
demo.add_population_split(time=363, ancestral="admixy", derived=["Rlllrlr"])
demo.add_population_split(time=33, ancestral="admixz", derived=["M"])
demo.add_population_split(time=735, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=664, ancestral="Rl", derived=["Rll","E"])
demo.add_population_split(time=642, ancestral="Rll", derived=["Rlll","Rllr"])
demo.add_population_split(time=616, ancestral="Rlll", derived=["Rllll","Rlllr"])
demo.add_population_split(time=570, ancestral="Rllll", derived=["K","Rlllls"])
demo.add_population_split(time=200, ancestral="Rllllr", derived=["Rllllrd","Rllllrf"])
demo.add_population_split(time=160, ancestral="Rllllrf", derived=["Rllllrs"])
demo.add_population_split(time=98, ancestral="Rllllrs", derived=["L","Rllllrsl"])
demo.add_population_split(time=570, ancestral="Rlllr", derived=["Rlllrl"])
demo.add_population_split(time=279, ancestral="Rlllrlr", derived=["A","Rlllrlrm"])
demo.add_population_split(time=200, ancestral="Rlllrlrm", derived=["Rlllrlry"])
demo.add_population_split(time=160, ancestral="Rlllrlry", derived=["G"])
demo.add_population_split(time=616, ancestral="Rllr", derived=["Rllra"])
demo.add_population_split(time=664, ancestral="Rr", derived=["Rrl","Rrp"])
demo.add_population_split(time=642, ancestral="Rrl", derived=["J","H"])
demo.add_admixture(time=56, derived="admix", ancestral=["admixe","Rlllr"], proportions=[0.5,0.5])
demo.add_admixture(time=200, derived="admixd", ancestral=["admixiu","Rlllrl"], proportions=[0.36,0.64])
demo.add_admixture(time=98, derived="admixe", ancestral=["Rllllrd","Rlllrlrm"], proportions=[0.3,0.7])
demo.add_admixture(time=363, derived="admixi", ancestral=["Rlllls","Rllra"], proportions=[0.33,0.67])
demo.add_admixture(time=33, derived="admixk", ancestral=["Rllllrsl","Rllllrd"], proportions=[0.76,0.24])
demo.add_admixture(time=570, derived="admixu", ancestral=["Rllr","Rrp"], proportions=[0.22,0.78])
demo.add_admixture(time=33, derived="admixw", ancestral=["Rlllrlry","admixd"], proportions=[0.6,0.4])
demo.add_admixture(time=98, derived="admixwh", ancestral=["Rllllrf","Rllra"], proportions=[0.2,0.8])
demo.add_admixture(time=363, derived="admixy", ancestral=["Rlllrl","Rlllls"], proportions=[0.77,0.23])
demo.add_admixture(time=33, derived="admixz", ancestral=["Rllllrsl","Rrp"], proportions=[0.25,0.75])
demo.sort_events()

r_chrom = 2e-08
r_break = math.log(2)
chrom_positions = [0, 1e+08, 2e+08, 3e+08]
map_positions = [0, 1e+08, 100000001, 2e+08, 200000001, 3e+08]
rates = [r_chrom, r_break, r_chrom, r_break, r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)

ts = msprime.sim_ancestry(samples=samples, model=[msprime.DiscreteTimeWrightFisher(duration=25), msprime.StandardCoalescent()], demography=demo, recombination_rate=rate_map)
tree_sequence = msprime.sim_mutations(ts, rate=1.25e-08, model='binary')
hap_gt = tree_sequence.genotype_matrix()
gt = hap_gt[:, range(0,nhap,2)] + hap_gt[:, range(1,nhap,2)]
nsnps = gt.shape[0]
ts_chroms = numpy.searchsorted(numpy.array(chrom_positions[1:]), tree_sequence.tables.sites.position, side='right') + 1

numpy.savetxt('./projects/qpadm/results/simout/3chr/sim_24/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/3chr/sim_24/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/3chr/sim_24/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

