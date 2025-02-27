import numpy
import math
import msprime
import multiprocess

indnam = ["L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("L", "I", "A", "E", "J", "C", "M", "H", "G", "B", "K", "F", "D"), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 4781, name = "R")
demo.add_population(initial_size = 3544, name = "Rl")
demo.add_population(initial_size = 2061, name = "Rll")
demo.add_population(initial_size = 7234, name = "L")
demo.add_population(initial_size = 3462, name = "Rlr")
demo.add_population(initial_size = 1691, name = "I")
demo.add_population(initial_size = 4081, name = "Rr")
demo.add_population(initial_size = 5117, name = "Rrr")
demo.add_population(initial_size = 2355, name = "Rrrl")
demo.add_population(initial_size = 3379, name = "Rrrll")
demo.add_population(initial_size = 9694, name = "A")
demo.add_population(initial_size = 7603, name = "E")
demo.add_population(initial_size = 8213, name = "Rrrlr")
demo.add_population(initial_size = 6278, name = "Rrrlrl")
demo.add_population(initial_size = 9259, name = "J")
demo.add_population(initial_size = 2222, name = "Rrrla")
demo.add_population(initial_size = 7929, name = "admix")
demo.add_population(initial_size = 4277, name = "Rrrr")
demo.add_population(initial_size = 8438, name = "Rlrw")
demo.add_population(initial_size = 8633, name = "admixs")
demo.add_population(initial_size = 8927, name = "Rrl")
demo.add_population(initial_size = 9758, name = "Rrrrf")
demo.add_population(initial_size = 6573, name = "C")
demo.add_population(initial_size = 9106, name = "admixc")
demo.add_population(initial_size = 7518, name = "M")
demo.add_population(initial_size = 1148, name = "Rrlw")
demo.add_population(initial_size = 7874, name = "H")
demo.add_population(initial_size = 9090, name = "admixj")
demo.add_population(initial_size = 6424, name = "G")
demo.add_population(initial_size = 8047, name = "Rrrlrld")
demo.add_population(initial_size = 5104, name = "B")
demo.add_population(initial_size = 2139, name = "admixr")
demo.add_population(initial_size = 7135, name = "Rru")
demo.add_population(initial_size = 6045, name = "admixn")
demo.add_population(initial_size = 2087, name = "Rlk")
demo.add_population(initial_size = 9240, name = "admixl")
demo.add_population(initial_size = 6847, name = "K")
demo.add_population(initial_size = 1115, name = "Rrrlax")
demo.add_population(initial_size = 3310, name = "admixf")
demo.add_population(initial_size = 6916, name = "admixfz")
demo.add_population(initial_size = 3013, name = "F")
demo.add_population(initial_size = 1681, name = "admixlh")
demo.add_population(initial_size = 7039, name = "D")
demo.add_population(initial_size = 8346, name = "Rlkj")
demo.add_population(initial_size = 3631, name = "admixw")
demo.add_population_split(time=353, ancestral="admix", derived=["Rrrr"])
demo.add_population_split(time=139, ancestral="admixc", derived=["M"])
demo.add_population_split(time=324, ancestral="admixf", derived=["admixfz"])
demo.add_population_split(time=212, ancestral="admixfz", derived=["F"])
demo.add_population_split(time=65, ancestral="admixj", derived=["G"])
demo.add_population_split(time=212, ancestral="admixl", derived=["K"])
demo.add_population_split(time=139, ancestral="admixlh", derived=["D"])
demo.add_population_split(time=324, ancestral="admixn", derived=["Rrrll"])
demo.add_population_split(time=430, ancestral="admixs", derived=["Rrl"])
demo.add_population_split(time=353, ancestral="admixw", derived=["Rrrlrl"])
demo.add_population_split(time=616, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=588, ancestral="Rl", derived=["Rlr","Rlk"])
demo.add_population_split(time=553, ancestral="Rlk", derived=["Rlkj"])
demo.add_population_split(time=454, ancestral="Rlkj", derived=["Rll"])
demo.add_population_split(time=430, ancestral="Rll", derived=["L"])
demo.add_population_split(time=553, ancestral="Rlr", derived=["I","Rlrw"])
demo.add_population_split(time=588, ancestral="Rr", derived=["Rrr","Rru"])
demo.add_population_split(time=353, ancestral="Rrl", derived=["Rrlw"])
demo.add_population_split(time=324, ancestral="Rrlw", derived=["H"])
demo.add_population_split(time=553, ancestral="Rrr", derived=["Rrrl"])
demo.add_population_split(time=454, ancestral="Rrrl", derived=["Rrrlr","Rrrla"])
demo.add_population_split(time=430, ancestral="Rrrla", derived=["Rrrlax"])
demo.add_population_split(time=212, ancestral="Rrrll", derived=["A","E"])
demo.add_population_split(time=324, ancestral="Rrrlrl", derived=["J","Rrrlrld"])
demo.add_population_split(time=212, ancestral="Rrrlrld", derived=["B"])
demo.add_population_split(time=324, ancestral="Rrrr", derived=["Rrrrf"])
demo.add_population_split(time=212, ancestral="Rrrrf", derived=["C"])
demo.add_admixture(time=353, derived="admix", ancestral=["Rrrla","Rrr"], proportions=[0.85,0.15])
demo.add_admixture(time=139, derived="admixc", ancestral=["Rrrrf","Rll"], proportions=[0.76,0.24])
demo.add_admixture(time=324, derived="admixf", ancestral=["Rrrlax","Rlrw"], proportions=[0.74,0.26])
demo.add_admixture(time=65, derived="admixj", ancestral=["admixr","Rrrlr"], proportions=[0.39,0.61])
demo.add_admixture(time=212, derived="admixl", ancestral=["Rrrr","Rlk"], proportions=[0.58,0.42])
demo.add_admixture(time=139, derived="admixlh", ancestral=["admixfz","Rrl"], proportions=[0.45,0.55])
demo.add_admixture(time=324, derived="admixn", ancestral=["Rrrlax","Rru"], proportions=[0.13,0.87])
demo.add_admixture(time=139, derived="admixr", ancestral=["Rrrlrld","Rrlw"], proportions=[0.66,0.34])
demo.add_admixture(time=430, derived="admixs", ancestral=["Rlrw","Rru"], proportions=[0.29,0.71])
demo.add_admixture(time=353, derived="admixw", ancestral=["Rrrlr","Rlkj"], proportions=[0.6,0.4])
demo.sort_events()

r_chrom = 2e-08
r_break = math.log(2)
chrom_positions = [0, 1e+08, 2e+08, 3e+08]
map_positions = [0, 1e+08, 100000001, 2e+08, 200000001, 3e+08]
rates = [r_chrom, r_break, r_chrom, r_break, r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)

ts = msprime.sim_ancestry(samples=samples, model=[msprime.DiscreteTimeWrightFisher(duration=25), msprime.StandardCoalescent()], demography=demo, recombination_rate=rate_map)
tree_sequence = msprime.sim_mutations(ts, rate=1.25e-08, model='binary')
hap_gt = tree_sequence.genotype_matrix()
gt = hap_gt[:, range(0,nhap,2)] + hap_gt[:, range(1,nhap,2)]
nsnps = gt.shape[0]
ts_chroms = numpy.searchsorted(numpy.array(chrom_positions[1:]), tree_sequence.tables.sites.position, side='right') + 1

numpy.savetxt('./projects/qpadm/results/simout/qpadm_sim_fixed_13.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_13.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_13.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

