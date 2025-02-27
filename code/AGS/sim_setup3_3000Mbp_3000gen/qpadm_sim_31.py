import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("J", "G", "B", "L", "M", "E", "D", "K", "C", "F", "A", "H", "I"), (1901, 2325, 701, 169, 1418, 1552, 169, 0, 1028, 169, 0, 169, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 1442, name = "R")
demo.add_population(initial_size = 6849, name = "Rl")
demo.add_population(initial_size = 3908, name = "Rll")
demo.add_population(initial_size = 1611, name = "J")
demo.add_population(initial_size = 5095, name = "G")
demo.add_population(initial_size = 3083, name = "Rr")
demo.add_population(initial_size = 5339, name = "Rrr")
demo.add_population(initial_size = 4558, name = "Rrrl")
demo.add_population(initial_size = 7058, name = "Rrrlr")
demo.add_population(initial_size = 2788, name = "Rrrlrl")
demo.add_population(initial_size = 7303, name = "B")
demo.add_population(initial_size = 8683, name = "Rrrr")
demo.add_population(initial_size = 6608, name = "Rrrrl")
demo.add_population(initial_size = 1264, name = "L")
demo.add_population(initial_size = 7353, name = "Rrrrr")
demo.add_population(initial_size = 6504, name = "M")
demo.add_population(initial_size = 4403, name = "Rllc")
demo.add_population(initial_size = 3499, name = "E")
demo.add_population(initial_size = 1327, name = "admix")
demo.add_population(initial_size = 4582, name = "Rrrlrlg")
demo.add_population(initial_size = 4033, name = "D")
demo.add_population(initial_size = 5598, name = "admixo")
demo.add_population(initial_size = 5809, name = "K")
demo.add_population(initial_size = 4790, name = "Rrrrv")
demo.add_population(initial_size = 9179, name = "admixp")
demo.add_population(initial_size = 9705, name = "C")
demo.add_population(initial_size = 7306, name = "Rrrll")
demo.add_population(initial_size = 7817, name = "admixx")
demo.add_population(initial_size = 6759, name = "Rrrh")
demo.add_population(initial_size = 3429, name = "Rrrhm")
demo.add_population(initial_size = 5179, name = "admixv")
demo.add_population(initial_size = 7793, name = "Rrrllz")
demo.add_population(initial_size = 9277, name = "F")
demo.add_population(initial_size = 1074, name = "admixvl")
demo.add_population(initial_size = 2154, name = "A")
demo.add_population(initial_size = 4285, name = "Rrrz")
demo.add_population(initial_size = 4183, name = "admixb")
demo.add_population(initial_size = 6200, name = "Rrrf")
demo.add_population(initial_size = 3017, name = "admixz")
demo.add_population(initial_size = 3281, name = "Rllcq")
demo.add_population(initial_size = 5558, name = "admixl")
demo.add_population(initial_size = 6217, name = "H")
demo.add_population(initial_size = 2712, name = "Rrrrx")
demo.add_population(initial_size = 4745, name = "admixr")
demo.add_population(initial_size = 6520, name = "I")
demo.add_population_split(time=2842, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=2419, ancestral="Rl", derived=["Rll","G"])
demo.add_population_split(time=2325, ancestral="Rll", derived=["J","Rllc"])
demo.add_population_split(time=1901, ancestral="Rllc", derived=["E","Rllcq"])
demo.add_population_split(time=2419, ancestral="Rr", derived=["Rrr"])
demo.add_population_split(time=2325, ancestral="Rrr", derived=["Rrrr","Rrrf"])
demo.add_population_split(time=1901, ancestral="Rrrf", derived=["Rrrz"])
demo.add_population_split(time=1418, ancestral="Rrrh", derived=["Rrrhm"])
demo.add_population_split(time=1316, ancestral="Rrrhm", derived=["Rrrl"])
demo.add_population_split(time=1028, ancestral="Rrrl", derived=["Rrrlr"])
demo.add_population_split(time=701, ancestral="Rrrll", derived=["Rrrllz"])
demo.add_population_split(time=326, ancestral="Rrrllz", derived=["F"])
demo.add_population_split(time=900, ancestral="Rrrlr", derived=["Rrrlrl","B"])
demo.add_population_split(time=701, ancestral="Rrrlrl", derived=["Rrrlrlg"])
demo.add_population_split(time=326, ancestral="Rrrlrlg", derived=["D"])
demo.add_population_split(time=1901, ancestral="Rrrr", derived=["Rrrrr","Rrrrx"])
demo.add_population_split(time=326, ancestral="Rrrrl", derived=["L"])
demo.add_population_split(time=1552, ancestral="Rrrrr", derived=["M"])
demo.add_population_split(time=1552, ancestral="Rrrrx", derived=["Rrrrv"])
demo.add_population_split(time=1552, ancestral="Rrrz", derived=["Rrrh"])
demo.add_population_split(time=701, ancestral="admix", derived=["Rrrrl"])
demo.add_population_split(time=900, ancestral="admixb", derived=["Rrrll"])
demo.add_population_split(time=326, ancestral="admixl", derived=["H"])
demo.add_population_split(time=169, ancestral="admixo", derived=["K"])
demo.add_population_split(time=1316, ancestral="admixp", derived=["C"])
demo.add_population_split(time=169, ancestral="admixr", derived=["I"])
demo.add_population_split(time=169, ancestral="admixvl", derived=["A"])
demo.add_admixture(time=701, derived="admix", ancestral=["admixv","Rrrrv"], proportions=[0.19,0.81])
demo.add_admixture(time=900, derived="admixb", ancestral=["Rrrl","Rrrz"], proportions=[0.47,0.53])
demo.add_admixture(time=326, derived="admixl", ancestral=["Rrrlrl","Rllcq"], proportions=[0.85,0.15])
demo.add_admixture(time=169, derived="admixo", ancestral=["Rrrlrlg","Rrrrr"], proportions=[0.27,0.73])
demo.add_admixture(time=1316, derived="admixp", ancestral=["Rrrrv","Rr"], proportions=[0.52,0.48])
demo.add_admixture(time=169, derived="admixr", ancestral=["Rrrrl","Rrrrx"], proportions=[0.27,0.73])
demo.add_admixture(time=900, derived="admixv", ancestral=["Rrrhm","admixz"], proportions=[0.63,0.37])
demo.add_admixture(time=169, derived="admixvl", ancestral=["admixx","Rrrllz"], proportions=[0.14,0.86])
demo.add_admixture(time=326, derived="admixx", ancestral=["Rrrll","Rrrh"], proportions=[0.79,0.21])
demo.add_admixture(time=1418, derived="admixz", ancestral=["Rllcq","Rrrf"], proportions=[0.46,0.54])
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

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_31/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_31/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_31/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

