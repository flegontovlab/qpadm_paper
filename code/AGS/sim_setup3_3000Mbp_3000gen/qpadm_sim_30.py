import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("L", "K", "D", "C", "E", "I", "H", "J", "M", "G", "B", "F", "A"), (1822, 1822, 1470, 2186, 1822, 0, 1635, 1129, 686, 412, 1129, 686, 2186) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 3803, name = "Rl")
demo.add_population(initial_size = 2017, name = "Rll")
demo.add_population(initial_size = 5116, name = "L")
demo.add_population(initial_size = 5456, name = "Rlr")
demo.add_population(initial_size = 6241, name = "K")
demo.add_population(initial_size = 4656, name = "R")
demo.add_population(initial_size = 1902, name = "Rr")
demo.add_population(initial_size = 4926, name = "Rrl")
demo.add_population(initial_size = 5458, name = "Rrlr")
demo.add_population(initial_size = 7223, name = "Rrlrl")
demo.add_population(initial_size = 6891, name = "D")
demo.add_population(initial_size = 7527, name = "Rrr")
demo.add_population(initial_size = 2412, name = "Rrrl")
demo.add_population(initial_size = 4626, name = "C")
demo.add_population(initial_size = 4404, name = "Rrrr")
demo.add_population(initial_size = 9210, name = "Rrrrl")
demo.add_population(initial_size = 3246, name = "E")
demo.add_population(initial_size = 3618, name = "Rrlrh")
demo.add_population(initial_size = 4897, name = "admix")
demo.add_population(initial_size = 2605, name = "Rrrx")
demo.add_population(initial_size = 5886, name = "admixy")
demo.add_population(initial_size = 1663, name = "I")
demo.add_population(initial_size = 6528, name = "Rlla")
demo.add_population(initial_size = 3684, name = "H")
demo.add_population(initial_size = 7471, name = "admixu")
demo.add_population(initial_size = 3331, name = "Rg")
demo.add_population(initial_size = 1310, name = "admixa")
demo.add_population(initial_size = 1636, name = "J")
demo.add_population(initial_size = 1727, name = "Rrlrw")
demo.add_population(initial_size = 9143, name = "M")
demo.add_population(initial_size = 7624, name = "admixi")
demo.add_population(initial_size = 1653, name = "G")
demo.add_population(initial_size = 2670, name = "Rq")
demo.add_population(initial_size = 6537, name = "admixt")
demo.add_population(initial_size = 9410, name = "Rrlrlq")
demo.add_population(initial_size = 9599, name = "B")
demo.add_population(initial_size = 7149, name = "admixv")
demo.add_population(initial_size = 6756, name = "Rqv")
demo.add_population(initial_size = 2533, name = "admixd")
demo.add_population(initial_size = 3155, name = "F")
demo.add_population(initial_size = 1602, name = "Rrrrlb")
demo.add_population(initial_size = 5198, name = "admixk")
demo.add_population(initial_size = 8954, name = "Rrlz")
demo.add_population(initial_size = 7817, name = "A")
demo.add_population(initial_size = 1129, name = "admixis")
demo.add_population_split(time=2985, ancestral="R", derived=["Rr","Rq"])
demo.add_population_split(time=2610, ancestral="Rg", derived=["Rl"])
demo.add_population_split(time=2344, ancestral="Rl", derived=["Rll","Rlr"])
demo.add_population_split(time=2186, ancestral="Rll", derived=["L","Rlla"])
demo.add_population_split(time=1822, ancestral="Rlla", derived=["H"])
demo.add_population_split(time=2186, ancestral="Rlr", derived=["K"])
demo.add_population_split(time=2704, ancestral="Rq", derived=["Rg","Rqv"])
demo.add_population_split(time=2704, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=2610, ancestral="Rrl", derived=["Rrlr","Rrlz"])
demo.add_population_split(time=1822, ancestral="Rrlrh", derived=["Rrlrl"])
demo.add_population_split(time=1635, ancestral="Rrlrl", derived=["D","Rrlrlq"])
demo.add_population_split(time=1470, ancestral="Rrlrlq", derived=["B"])
demo.add_population_split(time=1129, ancestral="Rrlrw", derived=["M"])
demo.add_population_split(time=2344, ancestral="Rrlz", derived=["A"])
demo.add_population_split(time=2610, ancestral="Rrr", derived=["Rrrl","Rrrx"])
demo.add_population_split(time=2344, ancestral="Rrrl", derived=["C"])
demo.add_population_split(time=2186, ancestral="Rrrr", derived=["Rrrrl","E"])
demo.add_population_split(time=1822, ancestral="Rrrrl", derived=["Rrrrlb"])
demo.add_population_split(time=2344, ancestral="Rrrx", derived=["Rrrr"])
demo.add_population_split(time=1470, ancestral="admixa", derived=["J"])
demo.add_population_split(time=1129, ancestral="admixd", derived=["F"])
demo.add_population_split(time=686, ancestral="admixi", derived=["G"])
demo.add_population_split(time=1470, ancestral="admixk", derived=["Rrlrw"])
demo.add_population_split(time=2186, ancestral="admixt", derived=["Rrlrh"])
demo.add_population_split(time=412, ancestral="admixy", derived=["I"])
demo.add_admixture(time=686, derived="admix", ancestral=["admixv","Rrrl"], proportions=[0.48,0.52])
demo.add_admixture(time=1470, derived="admixa", ancestral=["Rrrrl","admixis"], proportions=[0.77,0.23])
demo.add_admixture(time=1129, derived="admixd", ancestral=["admixu","Rqv"], proportions=[0.35,0.65])
demo.add_admixture(time=686, derived="admixi", ancestral=["Rrlrw","Rlr"], proportions=[0.7,0.3])
demo.add_admixture(time=2186, derived="admixis", ancestral=["Rrlz","Rg"], proportions=[0.68,0.32])
demo.add_admixture(time=1470, derived="admixk", ancestral=["Rrrrlb","Rrlr"], proportions=[0.55,0.45])
demo.add_admixture(time=2186, derived="admixt", ancestral=["Rrlr","Rqv"], proportions=[0.69,0.31])
demo.add_admixture(time=1470, derived="admixu", ancestral=["Rrrrlb","Rlla"], proportions=[0.59,0.41])
demo.add_admixture(time=1129, derived="admixv", ancestral=["Rrlrlq","Rrlrh"], proportions=[0.78,0.22])
demo.add_admixture(time=412, derived="admixy", ancestral=["admix","Rrrx"], proportions=[0.41,0.59])
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

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_30/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_30/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_30/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

