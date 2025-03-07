import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("A", "L", "J", "K", "I", "B", "C", "M", "E", "H", "F", "D", "G"), (1132, 938, 1924, 1924, 158, 1800, 1039, 0, 1132, 0, 514, 158, 1132) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 4971, name = "R")
demo.add_population(initial_size = 7627, name = "Rl")
demo.add_population(initial_size = 2558, name = "Rll")
demo.add_population(initial_size = 1268, name = "Rlr")
demo.add_population(initial_size = 2469, name = "Rlrl")
demo.add_population(initial_size = 2668, name = "A")
demo.add_population(initial_size = 2004, name = "Rr")
demo.add_population(initial_size = 9860, name = "Rrl")
demo.add_population(initial_size = 1405, name = "Rrll")
demo.add_population(initial_size = 2930, name = "Rrllr")
demo.add_population(initial_size = 5598, name = "Rrlr")
demo.add_population(initial_size = 5713, name = "L")
demo.add_population(initial_size = 2799, name = "Rrr")
demo.add_population(initial_size = 7095, name = "J")
demo.add_population(initial_size = 8915, name = "K")
demo.add_population(initial_size = 5329, name = "Rllr")
demo.add_population(initial_size = 2474, name = "Rllru")
demo.add_population(initial_size = 1661, name = "admix")
demo.add_population(initial_size = 8667, name = "I")
demo.add_population(initial_size = 2751, name = "Rrlrd")
demo.add_population(initial_size = 3294, name = "admixt")
demo.add_population(initial_size = 4459, name = "Rlle")
demo.add_population(initial_size = 3966, name = "B")
demo.add_population(initial_size = 7206, name = "admixq")
demo.add_population(initial_size = 9300, name = "Rrlrdy")
demo.add_population(initial_size = 5136, name = "admixj")
demo.add_population(initial_size = 8118, name = "Rlrls")
demo.add_population(initial_size = 9941, name = "C")
demo.add_population(initial_size = 8756, name = "admixz")
demo.add_population(initial_size = 2100, name = "M")
demo.add_population(initial_size = 9849, name = "Rlll")
demo.add_population(initial_size = 9041, name = "admixp")
demo.add_population(initial_size = 1917, name = "Rrllrb")
demo.add_population(initial_size = 2080, name = "E")
demo.add_population(initial_size = 3336, name = "admixa")
demo.add_population(initial_size = 4333, name = "Rllrue")
demo.add_population(initial_size = 1415, name = "admixqx")
demo.add_population(initial_size = 4203, name = "H")
demo.add_population(initial_size = 9890, name = "admixqy")
demo.add_population(initial_size = 4868, name = "F")
demo.add_population(initial_size = 3357, name = "admixe")
demo.add_population(initial_size = 1699, name = "D")
demo.add_population(initial_size = 8252, name = "Rrlrdb")
demo.add_population(initial_size = 5762, name = "admixzv")
demo.add_population(initial_size = 3086, name = "G")
demo.add_population_split(time=2719, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=2408, ancestral="Rl", derived=["Rll"])
demo.add_population_split(time=2070, ancestral="Rll", derived=["Rlle","Rlll"])
demo.add_population_split(time=1924, ancestral="Rlle", derived=["B"])
demo.add_population_split(time=1132, ancestral="Rllr", derived=["Rllru"])
demo.add_population_split(time=1039, ancestral="Rllru", derived=["Rllrue"])
demo.add_population_split(time=1594, ancestral="Rlr", derived=["Rlrl"])
demo.add_population_split(time=1414, ancestral="Rlrl", derived=["A","Rlrls"])
demo.add_population_split(time=1132, ancestral="Rlrls", derived=["C"])
demo.add_population_split(time=2408, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=2070, ancestral="Rrl", derived=["Rrlrd"])
demo.add_population_split(time=1800, ancestral="Rrll", derived=["Rrllr"])
demo.add_population_split(time=1594, ancestral="Rrllr", derived=["Rrllrb"])
demo.add_population_split(time=1414, ancestral="Rrllrb", derived=["E"])
demo.add_population_split(time=1039, ancestral="Rrlr", derived=["L"])
demo.add_population_split(time=1924, ancestral="Rrlrd", derived=["Rrll","Rrlrdb"])
demo.add_population_split(time=1800, ancestral="Rrlrdb", derived=["Rrlrdy"])
demo.add_population_split(time=2070, ancestral="Rrr", derived=["J","K"])
demo.add_population_split(time=514, ancestral="admix", derived=["I"])
demo.add_population_split(time=1132, ancestral="admixa", derived=["Rrlr"])
demo.add_population_split(time=514, ancestral="admixe", derived=["D"])
demo.add_population_split(time=1800, ancestral="admixp", derived=["Rlr"])
demo.add_population_split(time=1039, ancestral="admixq", derived=["admixqy"])
demo.add_population_split(time=158, ancestral="admixqx", derived=["H"])
demo.add_population_split(time=938, ancestral="admixqy", derived=["F"])
demo.add_population_split(time=1414, ancestral="admixt", derived=["Rllr"])
demo.add_population_split(time=158, ancestral="admixz", derived=["M"])
demo.add_population_split(time=1414, ancestral="admixzv", derived=["G"])
demo.add_admixture(time=514, derived="admix", ancestral=["Rllrue","Rrll"], proportions=[0.29,0.71])
demo.add_admixture(time=1132, derived="admixa", ancestral=["Rrllrb","Rrl"], proportions=[0.61,0.39])
demo.add_admixture(time=514, derived="admixe", ancestral=["admixqy","admixj"], proportions=[0.55,0.45])
demo.add_admixture(time=1414, derived="admixj", ancestral=["Rrllr","Rrlrdy"], proportions=[0.47,0.53])
demo.add_admixture(time=1800, derived="admixp", ancestral=["Rlll","Rl"], proportions=[0.71,0.29])
demo.add_admixture(time=1039, derived="admixq", ancestral=["Rllr","Rlle"], proportions=[0.59,0.41])
demo.add_admixture(time=158, derived="admixqx", ancestral=["Rllrue","Rrlr"], proportions=[0.25,0.75])
demo.add_admixture(time=1414, derived="admixt", ancestral=["Rrlrdy","Rlll"], proportions=[0.6,0.4])
demo.add_admixture(time=158, derived="admixz", ancestral=["Rllru","Rlrls"], proportions=[0.39,0.61])
demo.add_admixture(time=1414, derived="admixzv", ancestral=["Rlr","Rrlrdb"], proportions=[0.2,0.8])
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

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_23/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_23/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_23/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

