import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("A", "B", "I", "H", "L", "C", "E", "K", "J", "G", "F", "M", "D"), (63, 63, 377, 87, 162, 401, 20, 87, 20, 377, 123, 0, 123) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 8592, name = "Rll")
demo.add_population(initial_size = 4456, name = "A")
demo.add_population(initial_size = 7024, name = "B")
demo.add_population(initial_size = 1837, name = "Rl")
demo.add_population(initial_size = 9838, name = "Rlr")
demo.add_population(initial_size = 6932, name = "I")
demo.add_population(initial_size = 7235, name = "R")
demo.add_population(initial_size = 5269, name = "Rr")
demo.add_population(initial_size = 7705, name = "Rrl")
demo.add_population(initial_size = 8581, name = "Rrll")
demo.add_population(initial_size = 4709, name = "Rrlll")
demo.add_population(initial_size = 7934, name = "H")
demo.add_population(initial_size = 4546, name = "Rrlr")
demo.add_population(initial_size = 8765, name = "Rrlrl")
demo.add_population(initial_size = 6381, name = "Rrlrr")
demo.add_population(initial_size = 7124, name = "L")
demo.add_population(initial_size = 1683, name = "Rrr")
demo.add_population(initial_size = 8516, name = "C")
demo.add_population(initial_size = 3838, name = "Rk")
demo.add_population(initial_size = 9304, name = "admix")
demo.add_population(initial_size = 1858, name = "E")
demo.add_population(initial_size = 7765, name = "Rlz")
demo.add_population(initial_size = 1434, name = "admixo")
demo.add_population(initial_size = 5809, name = "K")
demo.add_population(initial_size = 3698, name = "Rrllz")
demo.add_population(initial_size = 6840, name = "admixv")
demo.add_population(initial_size = 5174, name = "J")
demo.add_population(initial_size = 6951, name = "Rrrk")
demo.add_population(initial_size = 1099, name = "G")
demo.add_population(initial_size = 9681, name = "admixs")
demo.add_population(initial_size = 6910, name = "Rrlrre")
demo.add_population(initial_size = 8018, name = "F")
demo.add_population(initial_size = 5025, name = "admixh")
demo.add_population(initial_size = 7001, name = "Rrlllr")
demo.add_population(initial_size = 8611, name = "admixc")
demo.add_population(initial_size = 9710, name = "M")
demo.add_population(initial_size = 7907, name = "Rrlly")
demo.add_population(initial_size = 1670, name = "admixsg")
demo.add_population(initial_size = 1102, name = "D")
demo.add_population(initial_size = 8307, name = "Rrlk")
demo.add_population(initial_size = 6989, name = "admixvc")
demo.add_population(initial_size = 9274, name = "Rrli")
demo.add_population(initial_size = 5939, name = "admixp")
demo.add_population(initial_size = 2598, name = "Rrlio")
demo.add_population(initial_size = 8717, name = "admixpb")
demo.add_population_split(time=63, ancestral="admix", derived=["E"])
demo.add_population_split(time=20, ancestral="admixc", derived=["M"])
demo.add_population_split(time=123, ancestral="admixh", derived=["Rll"])
demo.add_population_split(time=123, ancestral="admixo", derived=["K"])
demo.add_population_split(time=280, ancestral="admixp", derived=["Rrlly"])
demo.add_population_split(time=162, ancestral="admixsg", derived=["D"])
demo.add_population_split(time=63, ancestral="admixv", derived=["J"])
demo.add_population_split(time=627, ancestral="R", derived=["Rr","Rk"])
demo.add_population_split(time=585, ancestral="Rk", derived=["Rl"])
demo.add_population_split(time=504, ancestral="Rl", derived=["Rlr","Rlz"])
demo.add_population_split(time=87, ancestral="Rll", derived=["A","B"])
demo.add_population_split(time=401, ancestral="Rlr", derived=["I"])
demo.add_population_split(time=585, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=504, ancestral="Rrl", derived=["Rrll","Rrli"])
demo.add_population_split(time=401, ancestral="Rrli", derived=["Rrlk","Rrlio"])
demo.add_population_split(time=377, ancestral="Rrlk", derived=["Rrlr"])
demo.add_population_split(time=123, ancestral="Rrlll", derived=["H","Rrlllr"])
demo.add_population_split(time=241, ancestral="Rrlly", derived=["Rrllz"])
demo.add_population_split(time=162, ancestral="Rrllz", derived=["Rrlll"])
demo.add_population_split(time=280, ancestral="Rrlr", derived=["Rrlrl","Rrlrr"])
demo.add_population_split(time=241, ancestral="Rrlrr", derived=["L","Rrlrre"])
demo.add_population_split(time=162, ancestral="Rrlrre", derived=["F"])
demo.add_population_split(time=504, ancestral="Rrr", derived=["C","Rrrk"])
demo.add_population_split(time=401, ancestral="Rrrk", derived=["G"])
demo.add_admixture(time=63, derived="admix", ancestral=["Rrlllr","Rk"], proportions=[0.86,0.14])
demo.add_admixture(time=20, derived="admixc", ancestral=["admixpb","Rrll"], proportions=[0.27,0.73])
demo.add_admixture(time=123, derived="admixh", ancestral=["Rrlrre","Rlz"], proportions=[0.86,0.14])
demo.add_admixture(time=123, derived="admixo", ancestral=["admixs","Rlz"], proportions=[0.36,0.64])
demo.add_admixture(time=280, derived="admixp", ancestral=["Rrlio","Rrll"], proportions=[0.41,0.59])
demo.add_admixture(time=63, derived="admixpb", ancestral=["Rrlllr","Rrlio"], proportions=[0.79,0.21])
demo.add_admixture(time=162, derived="admixs", ancestral=["Rrlrl","Rrrk"], proportions=[0.18,0.82])
demo.add_admixture(time=162, derived="admixsg", ancestral=["Rrlly","Rlr"], proportions=[0.53,0.47])
demo.add_admixture(time=63, derived="admixv", ancestral=["admixvc","Rrlrl"], proportions=[0.2,0.8])
demo.add_admixture(time=123, derived="admixvc", ancestral=["Rrllz","Rrlk"], proportions=[0.83,0.17])
demo.sort_events()

r_chrom = 2e-08
r_break = math.log(2)
chrom_positions = [0, 1e+08, 2e+08, 3e+08, 4e+08, 5e+08, 6e+08, 7e+08, 8e+08, 9e+08, 1e+09]
map_positions = [0, 1e+08, 100000001, 2e+08, 200000001, 3e+08, 300000001, 4e+08, 400000001, 5e+08, 500000001, 6e+08, 600000001, 7e+08, 700000001, 8e+08, 800000001, 9e+08, 900000001, 1e+09]
rates = [r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)

ts = msprime.sim_ancestry(samples=samples, model=[msprime.DiscreteTimeWrightFisher(duration=25), msprime.StandardCoalescent()], demography=demo, recombination_rate=rate_map)
tree_sequence = msprime.sim_mutations(ts, rate=1.25e-08, model='binary')
hap_gt = tree_sequence.genotype_matrix()
gt = hap_gt[:, range(0,nhap,2)] + hap_gt[:, range(1,nhap,2)]
nsnps = gt.shape[0]
ts_chroms = numpy.searchsorted(numpy.array(chrom_positions[1:]), tree_sequence.tables.sites.position, side='right') + 1

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_16/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_16/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_16/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

