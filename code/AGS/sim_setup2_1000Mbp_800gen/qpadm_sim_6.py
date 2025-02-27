import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("L", "A", "I", "C", "F", "K", "B", "M", "E", "H", "J", "D", "G"), (150, 183, 223, 150, 150, 534, 534, 321, 321, 53, 300, 0, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 2535, name = "R")
demo.add_population(initial_size = 3216, name = "Rl")
demo.add_population(initial_size = 1653, name = "Rlr")
demo.add_population(initial_size = 2067, name = "Rlrl")
demo.add_population(initial_size = 8325, name = "Rlrll")
demo.add_population(initial_size = 3763, name = "Rlrlll")
demo.add_population(initial_size = 5802, name = "L")
demo.add_population(initial_size = 2331, name = "A")
demo.add_population(initial_size = 2753, name = "Rlrr")
demo.add_population(initial_size = 9957, name = "Rlrrr")
demo.add_population(initial_size = 1600, name = "I")
demo.add_population(initial_size = 9800, name = "Rr")
demo.add_population(initial_size = 2435, name = "Rrl")
demo.add_population(initial_size = 2715, name = "Rrlr")
demo.add_population(initial_size = 6302, name = "C")
demo.add_population(initial_size = 3387, name = "F")
demo.add_population(initial_size = 4137, name = "Rrr")
demo.add_population(initial_size = 2338, name = "K")
demo.add_population(initial_size = 7971, name = "Rli")
demo.add_population(initial_size = 8872, name = "B")
demo.add_population(initial_size = 2508, name = "admix")
demo.add_population(initial_size = 7901, name = "M")
demo.add_population(initial_size = 4442, name = "Rrls")
demo.add_population(initial_size = 7252, name = "admixa")
demo.add_population(initial_size = 1637, name = "E")
demo.add_population(initial_size = 4841, name = "Rlrrt")
demo.add_population(initial_size = 5244, name = "admixe")
demo.add_population(initial_size = 5422, name = "Rrra")
demo.add_population(initial_size = 9766, name = "admixeo")
demo.add_population(initial_size = 7704, name = "H")
demo.add_population(initial_size = 3321, name = "Rrrah")
demo.add_population(initial_size = 6909, name = "admixv")
demo.add_population(initial_size = 6685, name = "Rlrrti")
demo.add_population(initial_size = 3167, name = "J")
demo.add_population(initial_size = 1552, name = "admixt")
demo.add_population(initial_size = 3077, name = "D")
demo.add_population(initial_size = 6922, name = "Rrraho")
demo.add_population(initial_size = 5117, name = "admixr")
demo.add_population(initial_size = 2742, name = "Rrrahoe")
demo.add_population(initial_size = 2606, name = "admixq")
demo.add_population(initial_size = 9783, name = "admixei")
demo.add_population(initial_size = 6766, name = "admixel")
demo.add_population(initial_size = 2015, name = "admixell")
demo.add_population(initial_size = 1614, name = "admixj")
demo.add_population(initial_size = 8100, name = "G")
demo.add_population_split(time=429, ancestral="admix", derived=["M"])
demo.add_population_split(time=429, ancestral="admixa", derived=["E"])
demo.add_population_split(time=321, ancestral="admixe", derived=["admixei"])
demo.add_population_split(time=300, ancestral="admixei", derived=["Rlrll"])
demo.add_population_split(time=150, ancestral="admixel", derived=["admixell"])
demo.add_population_split(time=125, ancestral="admixeo", derived=["H"])
demo.add_population_split(time=53, ancestral="admixj", derived=["G"])
demo.add_population_split(time=223, ancestral="admixr", derived=["Rrlr"])
demo.add_population_split(time=53, ancestral="admixt", derived=["D"])
demo.add_population_split(time=321, ancestral="admixv", derived=["Rlrrr"])
demo.add_population_split(time=652, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=607, ancestral="Rl", derived=["Rlr","Rli"])
demo.add_population_split(time=573, ancestral="Rli", derived=["B"])
demo.add_population_split(time=573, ancestral="Rlr", derived=["Rlrl","Rlrr"])
demo.add_population_split(time=223, ancestral="Rlrll", derived=["Rlrlll","A"])
demo.add_population_split(time=183, ancestral="Rlrlll", derived=["L"])
demo.add_population_split(time=534, ancestral="Rlrr", derived=["Rlrrt"])
demo.add_population_split(time=300, ancestral="Rlrrr", derived=["I"])
demo.add_population_split(time=429, ancestral="Rlrrt", derived=["Rlrrti"])
demo.add_population_split(time=321, ancestral="Rlrrti", derived=["J"])
demo.add_population_split(time=607, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=573, ancestral="Rrl", derived=["Rrls"])
demo.add_population_split(time=183, ancestral="Rrlr", derived=["C","F"])
demo.add_population_split(time=573, ancestral="Rrr", derived=["K","Rrra"])
demo.add_population_split(time=534, ancestral="Rrra", derived=["Rrrah"])
demo.add_population_split(time=429, ancestral="Rrrah", derived=["Rrraho"])
demo.add_population_split(time=321, ancestral="Rrraho", derived=["Rrrahoe"])
demo.add_admixture(time=429, derived="admix", ancestral=["Rrra","Rli"], proportions=[0.54,0.46])
demo.add_admixture(time=429, derived="admixa", ancestral=["Rlrl","Rrls"], proportions=[0.16,0.84])
demo.add_admixture(time=321, derived="admixe", ancestral=["Rlrrt","Rlrl"], proportions=[0.88,0.12])
demo.add_admixture(time=150, derived="admixel", ancestral=["Rlrlll","admixei"], proportions=[0.54,0.46])
demo.add_admixture(time=125, derived="admixeo", ancestral=["admixq","Rrraho"], proportions=[0.29,0.71])
demo.add_admixture(time=53, derived="admixj", ancestral=["admixell","Rlrrr"], proportions=[0.22,0.78])
demo.add_admixture(time=223, derived="admixq", ancestral=["Rrrahoe","Rrl"], proportions=[0.8,0.2])
demo.add_admixture(time=223, derived="admixr", ancestral=["Rrrahoe","Rrls"], proportions=[0.52,0.48])
demo.add_admixture(time=53, derived="admixt", ancestral=["admixell","Rlrrti"], proportions=[0.41,0.59])
demo.add_admixture(time=321, derived="admixv", ancestral=["Rrrah","Rlrr"], proportions=[0.55,0.45])
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

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_6/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_6/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_6/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

