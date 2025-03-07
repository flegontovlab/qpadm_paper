import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("H", "E", "D", "J", "I", "M", "A", "L", "B", "F", "G", "K", "C"), (431, 431, 582, 582, 461, 0, 210, 266, 528, 461, 210, 210, 114) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 8972, name = "R")
demo.add_population(initial_size = 2680, name = "Rl")
demo.add_population(initial_size = 2953, name = "Rll")
demo.add_population(initial_size = 9970, name = "Rlll")
demo.add_population(initial_size = 6791, name = "Rlr")
demo.add_population(initial_size = 6060, name = "Rlrl")
demo.add_population(initial_size = 1889, name = "H")
demo.add_population(initial_size = 1992, name = "E")
demo.add_population(initial_size = 4166, name = "Rr")
demo.add_population(initial_size = 5341, name = "Rrl")
demo.add_population(initial_size = 6919, name = "D")
demo.add_population(initial_size = 9846, name = "J")
demo.add_population(initial_size = 5842, name = "Rrr")
demo.add_population(initial_size = 6678, name = "Rrrl")
demo.add_population(initial_size = 9479, name = "Rlld")
demo.add_population(initial_size = 2127, name = "I")
demo.add_population(initial_size = 1000, name = "admix")
demo.add_population(initial_size = 5824, name = "Rllll")
demo.add_population(initial_size = 3200, name = "Rlllrt")
demo.add_population(initial_size = 4197, name = "Rlllr")
demo.add_population(initial_size = 8816, name = "admixh")
demo.add_population(initial_size = 8297, name = "M")
demo.add_population(initial_size = 9885, name = "Rlrd")
demo.add_population(initial_size = 5309, name = "admixe")
demo.add_population(initial_size = 7273, name = "A")
demo.add_population(initial_size = 4858, name = "Rllly")
demo.add_population(initial_size = 4600, name = "admixm")
demo.add_population(initial_size = 4328, name = "L")
demo.add_population(initial_size = 8931, name = "Rrrg")
demo.add_population(initial_size = 9769, name = "B")
demo.add_population(initial_size = 7485, name = "admixv")
demo.add_population(initial_size = 4636, name = "Rlo")
demo.add_population(initial_size = 2655, name = "admixvt")
demo.add_population(initial_size = 1557, name = "F")
demo.add_population(initial_size = 2053, name = "Rlrdd")
demo.add_population(initial_size = 4402, name = "admixi")
demo.add_population(initial_size = 8748, name = "G")
demo.add_population(initial_size = 8325, name = "Rlllrw")
demo.add_population(initial_size = 7674, name = "admixq")
demo.add_population(initial_size = 1116, name = "Rloc")
demo.add_population(initial_size = 4826, name = "admixn")
demo.add_population(initial_size = 3392, name = "Rlllrwi")
demo.add_population(initial_size = 5634, name = "K")
demo.add_population(initial_size = 2518, name = "admixy")
demo.add_population(initial_size = 2449, name = "C")
demo.add_population_split(time=431, ancestral="admix", derived=["Rllll"])
demo.add_population_split(time=266, ancestral="admixe", derived=["A"])
demo.add_population_split(time=114, ancestral="admixh", derived=["M"])
demo.add_population_split(time=266, ancestral="admixi", derived=["G"])
demo.add_population_split(time=371, ancestral="admixm", derived=["L"])
demo.add_population_split(time=528, ancestral="admixvt", derived=["F"])
demo.add_population_split(time=210, ancestral="admixy", derived=["C"])
demo.add_population_split(time=775, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=715, ancestral="Rl", derived=["Rlr","Rlo"])
demo.add_population_split(time=582, ancestral="Rll", derived=["Rlll","Rlld"])
demo.add_population_split(time=528, ancestral="Rlld", derived=["I"])
demo.add_population_split(time=528, ancestral="Rlll", derived=["Rlllrt","Rllly"])
demo.add_population_split(time=431, ancestral="Rlllr", derived=["Rlllrw"])
demo.add_population_split(time=461, ancestral="Rlllrt", derived=["Rlllr"])
demo.add_population_split(time=371, ancestral="Rlllrw", derived=["Rlllrwi"])
demo.add_population_split(time=266, ancestral="Rlllrwi", derived=["K"])
demo.add_population_split(time=693, ancestral="Rlo", derived=["Rll","Rloc"])
demo.add_population_split(time=693, ancestral="Rlr", derived=["Rlrd"])
demo.add_population_split(time=582, ancestral="Rlrd", derived=["Rlrdd"])
demo.add_population_split(time=528, ancestral="Rlrdd", derived=["Rlrl"])
demo.add_population_split(time=461, ancestral="Rlrl", derived=["H","E"])
demo.add_population_split(time=715, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=693, ancestral="Rrl", derived=["D","J"])
demo.add_population_split(time=693, ancestral="Rrr", derived=["Rrrl","Rrrg"])
demo.add_population_split(time=582, ancestral="Rrrg", derived=["B"])
demo.add_admixture(time=431, derived="admix", ancestral=["Rllly","Rlld"], proportions=[0.6,0.4])
demo.add_admixture(time=266, derived="admixe", ancestral=["Rllll","Rlrd"], proportions=[0.46,0.54])
demo.add_admixture(time=114, derived="admixh", ancestral=["admixq","Rllll"], proportions=[0.12,0.88])
demo.add_admixture(time=266, derived="admixi", ancestral=["Rlllr","admixn"], proportions=[0.78,0.22])
demo.add_admixture(time=371, derived="admixm", ancestral=["admixv","Rrrl"], proportions=[0.44,0.56])
demo.add_admixture(time=461, derived="admixn", ancestral=["Rlrdd","Rloc"], proportions=[0.81,0.19])
demo.add_admixture(time=266, derived="admixq", ancestral=["Rlllrw","Rlllrt"], proportions=[0.51,0.49])
demo.add_admixture(time=431, derived="admixv", ancestral=["Rllly","Rrrg"], proportions=[0.38,0.62])
demo.add_admixture(time=528, derived="admixvt", ancestral=["Rrrl","Rloc"], proportions=[0.8,0.2])
demo.add_admixture(time=210, derived="admixy", ancestral=["Rlllrwi","Rlr"], proportions=[0.43,0.57])
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

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_36/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_36/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_36/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

