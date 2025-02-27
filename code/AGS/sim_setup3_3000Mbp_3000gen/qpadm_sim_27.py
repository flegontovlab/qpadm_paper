import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("C", "K", "E", "D", "F", "I", "M", "J", "B", "L", "H", "G", "A"), (1399, 1560, 1226, 1312, 1312, 1312, 0, 788, 1312, 896, 788, 788, 525) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 4640, name = "R")
demo.add_population(initial_size = 8518, name = "Rl")
demo.add_population(initial_size = 9288, name = "Rll")
demo.add_population(initial_size = 3394, name = "Rlll")
demo.add_population(initial_size = 6818, name = "C")
demo.add_population(initial_size = 9812, name = "Rllr")
demo.add_population(initial_size = 2309, name = "Rlr")
demo.add_population(initial_size = 4196, name = "K")
demo.add_population(initial_size = 1620, name = "Rrl")
demo.add_population(initial_size = 8324, name = "Rrlr")
demo.add_population(initial_size = 2569, name = "E")
demo.add_population(initial_size = 5457, name = "Rr")
demo.add_population(initial_size = 7913, name = "Rrr")
demo.add_population(initial_size = 2565, name = "Rrrl")
demo.add_population(initial_size = 8421, name = "D")
demo.add_population(initial_size = 5998, name = "F")
demo.add_population(initial_size = 9357, name = "Rrrr")
demo.add_population(initial_size = 1517, name = "Rllrx")
demo.add_population(initial_size = 8644, name = "I")
demo.add_population(initial_size = 8106, name = "admix")
demo.add_population(initial_size = 6054, name = "M")
demo.add_population(initial_size = 5350, name = "Rrf")
demo.add_population(initial_size = 3028, name = "admixo")
demo.add_population(initial_size = 7500, name = "Rlrb")
demo.add_population(initial_size = 2998, name = "admixn")
demo.add_population(initial_size = 8954, name = "J")
demo.add_population(initial_size = 5054, name = "Rrfm")
demo.add_population(initial_size = 9083, name = "admixz")
demo.add_population(initial_size = 5150, name = "Rlllk")
demo.add_population(initial_size = 8005, name = "B")
demo.add_population(initial_size = 7426, name = "admixy")
demo.add_population(initial_size = 3602, name = "L")
demo.add_population(initial_size = 4907, name = "Rrfmv")
demo.add_population(initial_size = 7161, name = "admixd")
demo.add_population(initial_size = 2377, name = "Rrfmvp")
demo.add_population(initial_size = 1501, name = "admixou")
demo.add_population(initial_size = 7858, name = "H")
demo.add_population(initial_size = 6676, name = "Rx")
demo.add_population(initial_size = 6687, name = "admixod")
demo.add_population(initial_size = 9464, name = "Rrlt")
demo.add_population(initial_size = 3314, name = "admixr")
demo.add_population(initial_size = 8254, name = "admixrx")
demo.add_population(initial_size = 3958, name = "G")
demo.add_population(initial_size = 9144, name = "admixm")
demo.add_population(initial_size = 2933, name = "A")
demo.add_population_split(time=2618, ancestral="R", derived=["Rl","Rx"])
demo.add_population_split(time=2310, ancestral="Rl", derived=["Rll","Rlr"])
demo.add_population_split(time=1976, ancestral="Rll", derived=["Rlll","Rllr"])
demo.add_population_split(time=1560, ancestral="Rlll", derived=["C","Rlllk"])
demo.add_population_split(time=1399, ancestral="Rlllk", derived=["B"])
demo.add_population_split(time=1560, ancestral="Rllr", derived=["Rllrx"])
demo.add_population_split(time=1399, ancestral="Rllrx", derived=["I"])
demo.add_population_split(time=1976, ancestral="Rlr", derived=["K","Rlrb"])
demo.add_population_split(time=1976, ancestral="Rr", derived=["Rrr","Rrf"])
demo.add_population_split(time=1560, ancestral="Rrf", derived=["Rrl","Rrfm"])
demo.add_population_split(time=1399, ancestral="Rrfm", derived=["Rrfmv"])
demo.add_population_split(time=1312, ancestral="Rrfmv", derived=["Rrfmvp"])
demo.add_population_split(time=1399, ancestral="Rrl", derived=["Rrlr","Rrlt"])
demo.add_population_split(time=1312, ancestral="Rrlr", derived=["E"])
demo.add_population_split(time=1560, ancestral="Rrr", derived=["Rrrl","Rrrr"])
demo.add_population_split(time=1399, ancestral="Rrrl", derived=["D","F"])
demo.add_population_split(time=2310, ancestral="Rx", derived=["Rr"])
demo.add_population_split(time=120, ancestral="admix", derived=["M"])
demo.add_population_split(time=788, ancestral="admixm", derived=["A"])
demo.add_population_split(time=896, ancestral="admixn", derived=["J"])
demo.add_population_split(time=896, ancestral="admixou", derived=["H"])
demo.add_population_split(time=1226, ancestral="admixr", derived=["admixrx"])
demo.add_population_split(time=896, ancestral="admixrx", derived=["G"])
demo.add_population_split(time=1226, ancestral="admixy", derived=["L"])
demo.add_admixture(time=120, derived="admix", ancestral=["admixo","Rlrb"], proportions=[0.19,0.81])
demo.add_admixture(time=1226, derived="admixd", ancestral=["Rrfmv","Rlrb"], proportions=[0.88,0.12])
demo.add_admixture(time=788, derived="admixm", ancestral=["admixrx","Rrrr"], proportions=[0.33,0.67])
demo.add_admixture(time=896, derived="admixn", ancestral=["admixd","Rrrr"], proportions=[0.29,0.71])
demo.add_admixture(time=439, derived="admixo", ancestral=["admixod","Rllrx"], proportions=[0.39,0.61])
demo.add_admixture(time=896, derived="admixod", ancestral=["Rrfmvp","Rx"], proportions=[0.37,0.63])
demo.add_admixture(time=896, derived="admixou", ancestral=["Rrfmvp","admixz"], proportions=[0.63,0.37])
demo.add_admixture(time=1226, derived="admixr", ancestral=["Rrlr","Rrlt"], proportions=[0.17,0.83])
demo.add_admixture(time=1226, derived="admixy", ancestral=["Rrlt","Rlllk"], proportions=[0.57,0.43])
demo.add_admixture(time=1312, derived="admixz", ancestral=["Rrfm","Rllr"], proportions=[0.68,0.32])
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

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_27/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_27/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_27/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

