import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("D", "A", "L", "G", "B", "J", "K", "E", "H", "F", "I", "C", "M"), (453, 383, 383, 248, 453, 552, 129, 52, 383, 129, 129, 0, 284) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 9843, name = "Rl")
demo.add_population(initial_size = 8783, name = "Rll")
demo.add_population(initial_size = 7739, name = "D")
demo.add_population(initial_size = 9157, name = "Rllr")
demo.add_population(initial_size = 3499, name = "Rlr")
demo.add_population(initial_size = 5080, name = "Rlrl")
demo.add_population(initial_size = 8257, name = "A")
demo.add_population(initial_size = 5617, name = "Rlrr")
demo.add_population(initial_size = 2724, name = "L")
demo.add_population(initial_size = 4678, name = "R")
demo.add_population(initial_size = 1409, name = "Rr")
demo.add_population(initial_size = 6456, name = "Rrl")
demo.add_population(initial_size = 3876, name = "Rrll")
demo.add_population(initial_size = 3396, name = "G")
demo.add_population(initial_size = 7298, name = "Rrllr")
demo.add_population(initial_size = 1886, name = "Rrlr")
demo.add_population(initial_size = 5941, name = "B")
demo.add_population(initial_size = 3593, name = "J")
demo.add_population(initial_size = 6242, name = "Rlrlh")
demo.add_population(initial_size = 7919, name = "admix")
demo.add_population(initial_size = 2124, name = "admixt")
demo.add_population(initial_size = 7728, name = "K")
demo.add_population(initial_size = 8514, name = "Rlrrq")
demo.add_population(initial_size = 8929, name = "Rlrrh")
demo.add_population(initial_size = 2869, name = "admixu")
demo.add_population(initial_size = 1014, name = "E")
demo.add_population(initial_size = 1544, name = "Rg")
demo.add_population(initial_size = 1255, name = "admixv")
demo.add_population(initial_size = 6323, name = "Rllrd")
demo.add_population(initial_size = 1004, name = "H")
demo.add_population(initial_size = 1474, name = "admixl")
demo.add_population(initial_size = 7197, name = "Rllrdn")
demo.add_population(initial_size = 5731, name = "admixz")
demo.add_population(initial_size = 5036, name = "F")
demo.add_population(initial_size = 4431, name = "Rllrdno")
demo.add_population(initial_size = 4023, name = "admixg")
demo.add_population(initial_size = 8585, name = "Rrllrz")
demo.add_population(initial_size = 8542, name = "I")
demo.add_population(initial_size = 9177, name = "admixf")
demo.add_population(initial_size = 8402, name = "C")
demo.add_population(initial_size = 5964, name = "Rlrlhm")
demo.add_population(initial_size = 9988, name = "M")
demo.add_population(initial_size = 8096, name = "admixp")
demo.add_population(initial_size = 2306, name = "Rla")
demo.add_population(initial_size = 2686, name = "admixh")
demo.add_population_split(time=341, ancestral="admix", derived=["Rrll"])
demo.add_population_split(time=52, ancestral="admixf", derived=["C"])
demo.add_population_split(time=284, ancestral="admixp", derived=["Rlrrq"])
demo.add_population_split(time=153, ancestral="admixt", derived=["K"])
demo.add_population_split(time=129, ancestral="admixu", derived=["E"])
demo.add_population_split(time=153, ancestral="admixz", derived=["F"])
demo.add_population_split(time=757, ancestral="R", derived=["Rr","Rg"])
demo.add_population_split(time=665, ancestral="Rg", derived=["Rl"])
demo.add_population_split(time=552, ancestral="Rl", derived=["Rll","Rla"])
demo.add_population_split(time=501, ancestral="Rla", derived=["Rlr"])
demo.add_population_split(time=501, ancestral="Rll", derived=["D","Rllr"])
demo.add_population_split(time=453, ancestral="Rllr", derived=["Rllrd"])
demo.add_population_split(time=428, ancestral="Rllrd", derived=["H","Rllrdn"])
demo.add_population_split(time=383, ancestral="Rllrdn", derived=["Rllrdno"])
demo.add_population_split(time=453, ancestral="Rlr", derived=["Rlrl","Rlrr"])
demo.add_population_split(time=428, ancestral="Rlrl", derived=["A","Rlrlh"])
demo.add_population_split(time=383, ancestral="Rlrlh", derived=["Rlrlhm"])
demo.add_population_split(time=341, ancestral="Rlrlhm", derived=["M"])
demo.add_population_split(time=428, ancestral="Rlrr", derived=["L","Rlrrh"])
demo.add_population_split(time=665, ancestral="Rr", derived=["Rrl","J"])
demo.add_population_split(time=552, ancestral="Rrl", derived=["Rrlr"])
demo.add_population_split(time=284, ancestral="Rrll", derived=["G","Rrllr"])
demo.add_population_split(time=248, ancestral="Rrllr", derived=["Rrllrz"])
demo.add_population_split(time=153, ancestral="Rrllrz", derived=["I"])
demo.add_population_split(time=501, ancestral="Rrlr", derived=["B"])
demo.add_admixture(time=341, derived="admix", ancestral=["Rlrlh","Rrl"], proportions=[0.37,0.63])
demo.add_admixture(time=52, derived="admixf", ancestral=["Rrllrz","Rlrrq"], proportions=[0.87,0.13])
demo.add_admixture(time=153, derived="admixg", ancestral=["Rrllr","Rllrdno"], proportions=[0.72,0.28])
demo.add_admixture(time=248, derived="admixh", ancestral=["admixv","Rla"], proportions=[0.12,0.88])
demo.add_admixture(time=341, derived="admixl", ancestral=["Rllrdn","Rrlr"], proportions=[0.4,0.6])
demo.add_admixture(time=284, derived="admixp", ancestral=["Rlrlhm","Rlrrh"], proportions=[0.13,0.87])
demo.add_admixture(time=153, derived="admixt", ancestral=["Rlrrq","Rllr"], proportions=[0.81,0.19])
demo.add_admixture(time=129, derived="admixu", ancestral=["admixg","Rlrrh"], proportions=[0.35,0.65])
demo.add_admixture(time=284, derived="admixv", ancestral=["admixl","Rg"], proportions=[0.18,0.82])
demo.add_admixture(time=153, derived="admixz", ancestral=["admixh","Rllrdno"], proportions=[0.46,0.54])
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

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_5/qpadm_sim_' + sim_id + 'geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_5/qpadm_sim_' + sim_id + 'snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_5/qpadm_sim_' + sim_id + 'ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

