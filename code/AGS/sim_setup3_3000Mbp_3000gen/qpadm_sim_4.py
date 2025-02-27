import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("B", "I", "G", "L", "M", "E", "H", "D", "C", "J", "A", "K", "F"), (2722, 1148, 622, 1395, 1395, 1508, 154, 1508, 154, 330, 2306, 0, 240) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 6851, name = "R")
demo.add_population(initial_size = 6824, name = "B")
demo.add_population(initial_size = 1531, name = "Rr")
demo.add_population(initial_size = 6124, name = "Rrl")
demo.add_population(initial_size = 6437, name = "Rrll")
demo.add_population(initial_size = 6171, name = "Rrllr")
demo.add_population(initial_size = 5830, name = "Rrllrl")
demo.add_population(initial_size = 9471, name = "I")
demo.add_population(initial_size = 4563, name = "Rrlr")
demo.add_population(initial_size = 6792, name = "Rrlrl")
demo.add_population(initial_size = 6085, name = "Rrlrll")
demo.add_population(initial_size = 1033, name = "Rrlrlll")
demo.add_population(initial_size = 5315, name = "G")
demo.add_population(initial_size = 9069, name = "Rrlrlr")
demo.add_population(initial_size = 5943, name = "Rrlrlrl")
demo.add_population(initial_size = 1596, name = "L")
demo.add_population(initial_size = 4170, name = "Rrlrllj")
demo.add_population(initial_size = 9713, name = "M")
demo.add_population(initial_size = 9749, name = "admix")
demo.add_population(initial_size = 3303, name = "Rrlry")
demo.add_population(initial_size = 6025, name = "E")
demo.add_population(initial_size = 3532, name = "admixk")
demo.add_population(initial_size = 7938, name = "Rrlrlllk")
demo.add_population(initial_size = 8672, name = "admixd")
demo.add_population(initial_size = 3491, name = "H")
demo.add_population(initial_size = 3922, name = "Rrlrlrlp")
demo.add_population(initial_size = 7542, name = "admixw")
demo.add_population(initial_size = 7594, name = "Rrlrv")
demo.add_population(initial_size = 1268, name = "admixm")
demo.add_population(initial_size = 5305, name = "D")
demo.add_population(initial_size = 8272, name = "Rrlrlllkm")
demo.add_population(initial_size = 3784, name = "admixx")
demo.add_population(initial_size = 7078, name = "C")
demo.add_population(initial_size = 1821, name = "admixj")
demo.add_population(initial_size = 5296, name = "admixwr")
demo.add_population(initial_size = 3901, name = "J")
demo.add_population(initial_size = 1648, name = "Rrn")
demo.add_population(initial_size = 6153, name = "A")
demo.add_population(initial_size = 8239, name = "admixl")
demo.add_population(initial_size = 1467, name = "K")
demo.add_population(initial_size = 7783, name = "Rrlrlrw")
demo.add_population(initial_size = 5819, name = "admixs")
demo.add_population(initial_size = 5277, name = "F")
demo.add_population(initial_size = 2786, name = "Rrllrls")
demo.add_population(initial_size = 3283, name = "admixf")
demo.add_population_split(time=2884, ancestral="R", derived=["B","Rr"])
demo.add_population_split(time=2722, ancestral="Rr", derived=["Rrl","Rrn"])
demo.add_population_split(time=2411, ancestral="Rrl", derived=["Rrll","Rrlr"])
demo.add_population_split(time=2306, ancestral="Rrll", derived=["Rrllr"])
demo.add_population_split(time=1395, ancestral="Rrllrl", derived=["I","Rrllrls"])
demo.add_population_split(time=2306, ancestral="Rrlr", derived=["Rrlrl","Rrlrv"])
demo.add_population_split(time=1931, ancestral="Rrlrl", derived=["Rrlrll","Rrlrlr"])
demo.add_population_split(time=1740, ancestral="Rrlrll", derived=["Rrlrllj"])
demo.add_population_split(time=1508, ancestral="Rrlrllj", derived=["M"])
demo.add_population_split(time=881, ancestral="Rrlrlll", derived=["G","Rrlrlllk"])
demo.add_population_split(time=622, ancestral="Rrlrlllk", derived=["Rrlrlllkm"])
demo.add_population_split(time=1740, ancestral="Rrlrlr", derived=["Rrlrlrl","Rrlrlrw"])
demo.add_population_split(time=1508, ancestral="Rrlrlrl", derived=["L","Rrlrlrlp"])
demo.add_population_split(time=1931, ancestral="Rrlrv", derived=["Rrlry"])
demo.add_population_split(time=1740, ancestral="Rrlry", derived=["E"])
demo.add_population_split(time=2411, ancestral="Rrn", derived=["A"])
demo.add_population_split(time=1395, ancestral="admix", derived=["admixj"])
demo.add_population_split(time=240, ancestral="admixd", derived=["H"])
demo.add_population_split(time=1508, ancestral="admixk", derived=["Rrllrl"])
demo.add_population_split(time=154, ancestral="admixl", derived=["K"])
demo.add_population_split(time=1740, ancestral="admixm", derived=["D"])
demo.add_population_split(time=330, ancestral="admixs", derived=["F"])
demo.add_population_split(time=1148, ancestral="admixw", derived=["Rrlrlll"])
demo.add_population_split(time=622, ancestral="admixwr", derived=["J"])
demo.add_population_split(time=240, ancestral="admixx", derived=["C"])
demo.add_admixture(time=1395, derived="admix", ancestral=["Rrlrllj","Rrllr"], proportions=[0.8,0.2])
demo.add_admixture(time=240, derived="admixd", ancestral=["Rrlrlllkm","Rrlrlrw"], proportions=[0.36,0.64])
demo.add_admixture(time=881, derived="admixf", ancestral=["Rrllrls","Rrn"], proportions=[0.48,0.52])
demo.add_admixture(time=1508, derived="admixk", ancestral=["Rrlry","Rrllr"], proportions=[0.88,0.12])
demo.add_admixture(time=154, derived="admixl", ancestral=["admixf","Rrlrlrlp"], proportions=[0.4,0.6])
demo.add_admixture(time=1740, derived="admixm", ancestral=["Rrlrv","Rrll"], proportions=[0.82,0.18])
demo.add_admixture(time=330, derived="admixs", ancestral=["Rrlrlllk","Rrlrlrw"], proportions=[0.3,0.7])
demo.add_admixture(time=1148, derived="admixw", ancestral=["Rrlrlrlp","Rrlrll"], proportions=[0.52,0.48])
demo.add_admixture(time=622, derived="admixwr", ancestral=["admixj","Rrllrls"], proportions=[0.29,0.71])
demo.add_admixture(time=240, derived="admixx", ancestral=["Rrlrlllkm","admixj"], proportions=[0.53,0.47])
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

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_4/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_4/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_4/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

