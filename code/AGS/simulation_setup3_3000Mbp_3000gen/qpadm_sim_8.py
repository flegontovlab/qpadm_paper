import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("J", "F", "A", "C", "H", "I", "K", "E", "M", "L", "G", "D", "B"), (1132, 1132, 1894, 1894, 1894, 0, 0, 1365, 1132, 1132, 371, 371, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 7135, name = "Rll")
demo.add_population(initial_size = 8689, name = "Rlll")
demo.add_population(initial_size = 5549, name = "J")
demo.add_population(initial_size = 7274, name = "F")
demo.add_population(initial_size = 2247, name = "R")
demo.add_population(initial_size = 1955, name = "Rr")
demo.add_population(initial_size = 9131, name = "Rrl")
demo.add_population(initial_size = 2359, name = "A")
demo.add_population(initial_size = 4666, name = "Rrr")
demo.add_population(initial_size = 3858, name = "Rrrl")
demo.add_population(initial_size = 8032, name = "C")
demo.add_population(initial_size = 6451, name = "H")
demo.add_population(initial_size = 3865, name = "Rrrrl")
demo.add_population(initial_size = 4557, name = "I")
demo.add_population(initial_size = 7882, name = "K")
demo.add_population(initial_size = 7788, name = "Rrrr")
demo.add_population(initial_size = 3552, name = "Rrrrr")
demo.add_population(initial_size = 1102, name = "E")
demo.add_population(initial_size = 9421, name = "Rrrrrj")
demo.add_population(initial_size = 5747, name = "M")
demo.add_population(initial_size = 7060, name = "admix")
demo.add_population(initial_size = 2632, name = "admixk")
demo.add_population(initial_size = 8902, name = "Rlr")
demo.add_population(initial_size = 8524, name = "Rrn")
demo.add_population(initial_size = 4288, name = "admixl")
demo.add_population(initial_size = 9700, name = "Rm")
demo.add_population(initial_size = 5127, name = "Rl")
demo.add_population(initial_size = 7901, name = "admixkn")
demo.add_population(initial_size = 4674, name = "Rlrj")
demo.add_population(initial_size = 6064, name = "L")
demo.add_population(initial_size = 3796, name = "admixp")
demo.add_population(initial_size = 7987, name = "Rmv")
demo.add_population(initial_size = 7733, name = "admixj")
demo.add_population(initial_size = 7555, name = "Rz")
demo.add_population(initial_size = 1725, name = "Rllr")
demo.add_population(initial_size = 6256, name = "admixh")
demo.add_population(initial_size = 1259, name = "G")
demo.add_population(initial_size = 7013, name = "Rli")
demo.add_population(initial_size = 6127, name = "admixknv")
demo.add_population(initial_size = 6058, name = "Rllrt")
demo.add_population(initial_size = 5668, name = "D")
demo.add_population(initial_size = 6877, name = "admixr")
demo.add_population(initial_size = 9605, name = "B")
demo.add_population(initial_size = 1037, name = "Rln")
demo.add_population(initial_size = 1395, name = "admixkk")
demo.add_population_split(time=2678, ancestral="R", derived=["Rr","Rz"])
demo.add_population_split(time=2044, ancestral="Rl", derived=["Rli","Rln"])
demo.add_population_split(time=1894, ancestral="Rli", derived=["Rll"])
demo.add_population_split(time=1541, ancestral="Rll", derived=["Rlll"])
demo.add_population_split(time=1365, ancestral="Rlll", derived=["J","F"])
demo.add_population_split(time=600, ancestral="Rllr", derived=["Rllrt"])
demo.add_population_split(time=495, ancestral="Rllrt", derived=["D"])
demo.add_population_split(time=1894, ancestral="Rln", derived=["Rlr"])
demo.add_population_split(time=1541, ancestral="Rlr", derived=["Rlrj"])
demo.add_population_split(time=1365, ancestral="Rlrj", derived=["L"])
demo.add_population_split(time=2182, ancestral="Rm", derived=["Rl","Rmv"])
demo.add_population_split(time=2550, ancestral="Rr", derived=["Rrr","Rrn"])
demo.add_population_split(time=2044, ancestral="Rrl", derived=["A"])
demo.add_population_split(time=2182, ancestral="Rrn", derived=["Rrl"])
demo.add_population_split(time=2182, ancestral="Rrr", derived=["Rrrl"])
demo.add_population_split(time=2044, ancestral="Rrrl", derived=["C","H"])
demo.add_population_split(time=1894, ancestral="Rrrr", derived=["Rrrrr"])
demo.add_population_split(time=371, ancestral="Rrrrl", derived=["I","K"])
demo.add_population_split(time=1541, ancestral="Rrrrr", derived=["E","Rrrrrj"])
demo.add_population_split(time=1365, ancestral="Rrrrrj", derived=["M"])
demo.add_population_split(time=2550, ancestral="Rz", derived=["Rm"])
demo.add_population_split(time=495, ancestral="admixh", derived=["G"])
demo.add_population_split(time=731, ancestral="admixkk", derived=["Rllr"])
demo.add_population_split(time=2044, ancestral="admixl", derived=["Rrrr"])
demo.add_population_split(time=495, ancestral="admixp", derived=["Rrrrl"])
demo.add_population_split(time=371, ancestral="admixr", derived=["B"])
demo.add_admixture(time=945, derived="admix", ancestral=["admixknv","Rll"], proportions=[0.19,0.81])
demo.add_admixture(time=495, derived="admixh", ancestral=["Rllr","Rz"], proportions=[0.31,0.69])
demo.add_admixture(time=1132, derived="admixj", ancestral=["Rlrj","Rmv"], proportions=[0.45,0.55])
demo.add_admixture(time=1365, derived="admixk", ancestral=["Rlr","Rrrr"], proportions=[0.29,0.71])
demo.add_admixture(time=731, derived="admixkk", ancestral=["admix","Rln"], proportions=[0.42,0.58])
demo.add_admixture(time=1132, derived="admixkn", ancestral=["admixk","Rmv"], proportions=[0.32,0.68])
demo.add_admixture(time=1132, derived="admixknv", ancestral=["Rrrrrj","Rli"], proportions=[0.74,0.26])
demo.add_admixture(time=2044, derived="admixl", ancestral=["Rrr","Rrn"], proportions=[0.66,0.34])
demo.add_admixture(time=495, derived="admixp", ancestral=["admixkn","admixj"], proportions=[0.57,0.43])
demo.add_admixture(time=371, derived="admixr", ancestral=["Rllrt","Rrl"], proportions=[0.2,0.8])
demo.sort_events()

r_chrom = 2e-08
r_break = math.log(2)
chrom_positions = [0, 1e+08, 2e+08, 3e+08, 4e+08, 5e+08, 6e+08, 7e+08, 8e+08, 9e+08, 1e+09, 1.1e+09, 1.2e+09, 1.3e+09, 1.4e+09, 1.5e+09, 1.6e+09, 1.7e+09, 1.8e+09, 1.9e+09, 2e+09, 2.1e+09, 2.2e+09, 2.3e+09, 2.4e+09, 2.5e+09, 2.6e+09, 2.7e+09, 2.8e+09, 2.9e+09, 3e+09]
map_positions = [0, 1e+08, 100000001, 2e+08, 200000001, 3e+08, 300000001, 4e+08, 400000001, 5e+08, 500000001, 6e+08, 600000001, 7e+08, 700000001, 8e+08, 800000001, 9e+08, 900000001, 1e+09, 1000000001, 1.1e+09, 1100000001, 1.2e+09, 1200000001, 1.3e+09, 1300000001, 1.4e+09, 1400000001, 1.5e+09, 1500000001, 1.6e+09, 1600000001, 1.7e+09, 1700000001, 1.8e+09, 1800000001, 1.9e+09, 1900000001, 2e+09, 2000000001, 2.1e+09, 2100000001, 2.2e+09, 2200000001, 2.3e+09, 2300000001, 2.4e+09, 2400000001, 2.5e+09, 2500000001, 2.6e+09, 2600000001, 2.7e+09, 2700000001, 2.8e+09, 2800000001, 2.9e+09, 2900000001, 3e+09]
rates = [r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)

ts = msprime.sim_ancestry(samples=samples, model=[msprime.DiscreteTimeWrightFisher(duration=25), msprime.StandardCoalescent()], demography=demo, recombination_rate=rate_map)
tree_sequence = msprime.sim_mutations(ts, rate=1.25e-08, model='binary')
hap_gt = tree_sequence.genotype_matrix()
gt = hap_gt[:, range(0,nhap,2)] + hap_gt[:, range(1,nhap,2)]
nsnps = gt.shape[0]
ts_chroms = numpy.searchsorted(numpy.array(chrom_positions[1:]), tree_sequence.tables.sites.position, side='right') + 1

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_8/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_8/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_8/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

