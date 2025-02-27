import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("A", "B", "G", "M", "I", "L", "E", "D", "J", "C", "F", "H", "K"), (1721, 1721, 589, 589, 184, 761, 2164, 589, 0, 761, 1721, 184, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 6837, name = "Rl")
demo.add_population(initial_size = 2282, name = "Rlr")
demo.add_population(initial_size = 1473, name = "A")
demo.add_population(initial_size = 1782, name = "B")
demo.add_population(initial_size = 9731, name = "R")
demo.add_population(initial_size = 4649, name = "Rr")
demo.add_population(initial_size = 6170, name = "Rrlrl")
demo.add_population(initial_size = 1333, name = "G")
demo.add_population(initial_size = 8302, name = "Rrlr")
demo.add_population(initial_size = 9960, name = "Rrlrr")
demo.add_population(initial_size = 8688, name = "Rrlrrl")
demo.add_population(initial_size = 8934, name = "Rrlrrll")
demo.add_population(initial_size = 5930, name = "M")
demo.add_population(initial_size = 6706, name = "Rrlrrllr")
demo.add_population(initial_size = 2055, name = "I")
demo.add_population(initial_size = 8383, name = "L")
demo.add_population(initial_size = 2215, name = "E")
demo.add_population(initial_size = 5812, name = "Ru")
demo.add_population(initial_size = 9946, name = "Rup")
demo.add_population(initial_size = 5254, name = "Rrll")
demo.add_population(initial_size = 8549, name = "admixg")
demo.add_population(initial_size = 8211, name = "D")
demo.add_population(initial_size = 8204, name = "Ruu")
demo.add_population(initial_size = 2580, name = "admixy")
demo.add_population(initial_size = 6204, name = "J")
demo.add_population(initial_size = 8974, name = "admix")
demo.add_population(initial_size = 8074, name = "admixr")
demo.add_population(initial_size = 2033, name = "admixl")
demo.add_population(initial_size = 7103, name = "C")
demo.add_population(initial_size = 4570, name = "admixrj")
demo.add_population(initial_size = 5442, name = "admixlf")
demo.add_population(initial_size = 5014, name = "Rrl")
demo.add_population(initial_size = 6802, name = "Rrlh")
demo.add_population(initial_size = 7504, name = "admixa")
demo.add_population(initial_size = 4998, name = "admixrb")
demo.add_population(initial_size = 1666, name = "admixq")
demo.add_population(initial_size = 2285, name = "Rln")
demo.add_population(initial_size = 6136, name = "F")
demo.add_population(initial_size = 5290, name = "admixk")
demo.add_population(initial_size = 6519, name = "Rrc")
demo.add_population(initial_size = 2686, name = "admixu")
demo.add_population(initial_size = 7816, name = "Rrlrlj")
demo.add_population(initial_size = 9723, name = "H")
demo.add_population(initial_size = 5443, name = "admixd")
demo.add_population(initial_size = 7224, name = "K")
demo.add_population_split(time=2835, ancestral="R", derived=["Rr","Ru"])
demo.add_population_split(time=2029, ancestral="Rl", derived=["Rlr","Rln"])
demo.add_population_split(time=1924, ancestral="Rln", derived=["F"])
demo.add_population_split(time=1924, ancestral="Rlr", derived=["A","B"])
demo.add_population_split(time=2426, ancestral="Rr", derived=["E","Rrc"])
demo.add_population_split(time=2164, ancestral="Rrc", derived=["Rrl"])
demo.add_population_split(time=2029, ancestral="Rrl", derived=["Rrlh"])
demo.add_population_split(time=1485, ancestral="Rrlr", derived=["Rrlrr"])
demo.add_population_split(time=761, ancestral="Rrlrl", derived=["G","Rrlrlj"])
demo.add_population_split(time=589, ancestral="Rrlrlj", derived=["H"])
demo.add_population_split(time=1069, ancestral="Rrlrr", derived=["Rrlrrl"])
demo.add_population_split(time=982, ancestral="Rrlrrl", derived=["Rrlrrll","L"])
demo.add_population_split(time=761, ancestral="Rrlrrll", derived=["M","Rrlrrllr"])
demo.add_population_split(time=589, ancestral="Rrlrrllr", derived=["I"])
demo.add_population_split(time=2426, ancestral="Ru", derived=["Rup","Ruu"])
demo.add_population_split(time=2164, ancestral="Rup", derived=["Rl"])
demo.add_population_split(time=1721, ancestral="admix", derived=["admixr"])
demo.add_population_split(time=184, ancestral="admixd", derived=["K"])
demo.add_population_split(time=761, ancestral="admixg", derived=["D"])
demo.add_population_split(time=982, ancestral="admixl", derived=["C"])
demo.add_population_split(time=982, ancestral="admixlf", derived=["Rrlrl"])
demo.add_population_split(time=1069, ancestral="admixq", derived=["Rrll"])
demo.add_population_split(time=1631, ancestral="admixr", derived=["Rrlr","admixrb"])
demo.add_population_split(time=1485, ancestral="admixrb", derived=["admixrj"])
demo.add_population_split(time=184, ancestral="admixy", derived=["J"])
demo.add_admixture(time=1721, derived="admix", ancestral=["Rrlh","admixu"], proportions=[0.79,0.21])
demo.add_admixture(time=1721, derived="admixa", ancestral=["Rrlh","Ruu"], proportions=[0.42,0.58])
demo.add_admixture(time=184, derived="admixd", ancestral=["Rrlrrllr","Rrlrlj"], proportions=[0.63,0.37])
demo.add_admixture(time=761, derived="admixg", ancestral=["Rrll","Rup"], proportions=[0.26,0.74])
demo.add_admixture(time=1631, derived="admixk", ancestral=["admixa","Rln"], proportions=[0.5,0.5])
demo.add_admixture(time=982, derived="admixl", ancestral=["Rrlrr","admixrj"], proportions=[0.53,0.47])
demo.add_admixture(time=982, derived="admixlf", ancestral=["admixrj","Rrlr"], proportions=[0.11,0.89])
demo.add_admixture(time=1069, derived="admixq", ancestral=["admixrb","Rrl"], proportions=[0.3,0.7])
demo.add_admixture(time=2029, derived="admixu", ancestral=["Ruu","Rrc"], proportions=[0.81,0.19])
demo.add_admixture(time=184, derived="admixy", ancestral=["Rrll","admixk"], proportions=[0.65,0.35])
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

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_1/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_1/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_1/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

