import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("C", "A", "K", "M", "G", "J", "D", "L", "E", "I", "B", "F", "H"), (608, 2006, 829, 1361, 608, 608, 345, 191, 608, 345, 608, 98, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 4517, name = "R")
demo.add_population(initial_size = 2114, name = "Rl")
demo.add_population(initial_size = 5393, name = "Rll")
demo.add_population(initial_size = 1645, name = "Rllr")
demo.add_population(initial_size = 1962, name = "Rllrl")
demo.add_population(initial_size = 7548, name = "C")
demo.add_population(initial_size = 3041, name = "A")
demo.add_population(initial_size = 2229, name = "Rr")
demo.add_population(initial_size = 5673, name = "Rrl")
demo.add_population(initial_size = 5913, name = "Rrll")
demo.add_population(initial_size = 8148, name = "Rrlll")
demo.add_population(initial_size = 4350, name = "Rrlllr")
demo.add_population(initial_size = 4857, name = "K")
demo.add_population(initial_size = 1840, name = "M")
demo.add_population(initial_size = 3675, name = "Rrlr")
demo.add_population(initial_size = 7758, name = "Rrlrr")
demo.add_population(initial_size = 9631, name = "G")
demo.add_population(initial_size = 7251, name = "J")
demo.add_population(initial_size = 6737, name = "Rrlrf")
demo.add_population(initial_size = 2095, name = "admix")
demo.add_population(initial_size = 1501, name = "D")
demo.add_population(initial_size = 5864, name = "Rllrd")
demo.add_population(initial_size = 5235, name = "admixp")
demo.add_population(initial_size = 7725, name = "L")
demo.add_population(initial_size = 8052, name = "Rllb")
demo.add_population(initial_size = 9640, name = "E")
demo.add_population(initial_size = 7268, name = "admixv")
demo.add_population(initial_size = 7562, name = "Rrlrrs")
demo.add_population(initial_size = 9478, name = "admixw")
demo.add_population(initial_size = 1641, name = "Rrlh")
demo.add_population(initial_size = 8187, name = "admixj")
demo.add_population(initial_size = 6705, name = "Rrlhx")
demo.add_population(initial_size = 6549, name = "admixpb")
demo.add_population(initial_size = 1414, name = "I")
demo.add_population(initial_size = 6770, name = "Rly")
demo.add_population(initial_size = 2821, name = "admixn")
demo.add_population(initial_size = 8272, name = "admixjw")
demo.add_population(initial_size = 4670, name = "B")
demo.add_population(initial_size = 7838, name = "admixe")
demo.add_population(initial_size = 4883, name = "admixem")
demo.add_population(initial_size = 5390, name = "admixi")
demo.add_population(initial_size = 9272, name = "F")
demo.add_population(initial_size = 6874, name = "Rrlllrm")
demo.add_population(initial_size = 2862, name = "admixu")
demo.add_population(initial_size = 5500, name = "H")
demo.add_population_split(time=2809, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=2411, ancestral="Rl", derived=["A","Rly"])
demo.add_population_split(time=1751, ancestral="Rll", derived=["Rllr"])
demo.add_population_split(time=829, ancestral="Rllb", derived=["E"])
demo.add_population_split(time=1361, ancestral="Rllr", derived=["Rllrd"])
demo.add_population_split(time=1132, ancestral="Rllrd", derived=["Rllrl"])
demo.add_population_split(time=829, ancestral="Rllrl", derived=["C"])
demo.add_population_split(time=2006, ancestral="Rly", derived=["Rll"])
demo.add_population_split(time=2411, ancestral="Rr", derived=["Rrl"])
demo.add_population_split(time=2006, ancestral="Rrl", derived=["Rrll","Rrlh"])
demo.add_population_split(time=1751, ancestral="Rrlh", derived=["Rrlrrs","Rrlhx"])
demo.add_population_split(time=1751, ancestral="Rrll", derived=["Rrlll","M"])
demo.add_population_split(time=1361, ancestral="Rrlll", derived=["Rrlllr"])
demo.add_population_split(time=1132, ancestral="Rrlllr", derived=["K","Rrlllrm"])
demo.add_population_split(time=1132, ancestral="Rrlr", derived=["Rrlrr","Rrlrf"])
demo.add_population_split(time=829, ancestral="Rrlrr", derived=["G","J"])
demo.add_population_split(time=1361, ancestral="Rrlrrs", derived=["Rrlr"])
demo.add_population_split(time=608, ancestral="admix", derived=["D"])
demo.add_population_split(time=608, ancestral="admixe", derived=["admixem"])
demo.add_population_split(time=191, ancestral="admixi", derived=["F"])
demo.add_population_split(time=1132, ancestral="admixj", derived=["admixjw"])
demo.add_population_split(time=829, ancestral="admixjw", derived=["B"])
demo.add_population_split(time=345, ancestral="admixp", derived=["L"])
demo.add_population_split(time=608, ancestral="admixpb", derived=["I"])
demo.add_population_split(time=98, ancestral="admixu", derived=["H"])
demo.add_population_split(time=1132, ancestral="admixw", derived=["Rllb"])
demo.add_admixture(time=608, derived="admix", ancestral=["Rrlrf","Rrlllrm"], proportions=[0.78,0.22])
demo.add_admixture(time=608, derived="admixe", ancestral=["admixjw","Rly"], proportions=[0.14,0.86])
demo.add_admixture(time=191, derived="admixi", ancestral=["admixem","Rllrl"], proportions=[0.19,0.81])
demo.add_admixture(time=1132, derived="admixj", ancestral=["Rllr","Rrlhx"], proportions=[0.49,0.51])
demo.add_admixture(time=191, derived="admixn", ancestral=["admixem","Rrlll"], proportions=[0.35,0.65])
demo.add_admixture(time=345, derived="admixp", ancestral=["admixv","Rr"], proportions=[0.19,0.81])
demo.add_admixture(time=608, derived="admixpb", ancestral=["Rrlrf","Rrlhx"], proportions=[0.55,0.45])
demo.add_admixture(time=98, derived="admixu", ancestral=["admixn","Rrlllrm"], proportions=[0.5,0.5])
demo.add_admixture(time=608, derived="admixv", ancestral=["Rllb","Rllrd"], proportions=[0.28,0.72])
demo.add_admixture(time=1132, derived="admixw", ancestral=["Rrlrrs","Rll"], proportions=[0.75,0.25])
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

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_21/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_21/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_21/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

