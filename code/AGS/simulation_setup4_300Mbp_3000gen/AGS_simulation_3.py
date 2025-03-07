import numpy
import math
import msprime
import multiprocess

indnam = ["G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("G", "C", "I", "F", "M", "B", "D", "K", "E", "J", "L", "H", "A"), (1766, 810, 2032, 810, 1279, 390, 2430, 810, 2032, 1110, 810, 274, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 1506, name = "R")
demo.add_population(initial_size = 9425, name = "Rl")
demo.add_population(initial_size = 8556, name = "Rll")
demo.add_population(initial_size = 5652, name = "Rllr")
demo.add_population(initial_size = 1939, name = "Rllrl")
demo.add_population(initial_size = 6228, name = "G")
demo.add_population(initial_size = 3569, name = "Rllrlr")
demo.add_population(initial_size = 8535, name = "C")
demo.add_population(initial_size = 4180, name = "I")
demo.add_population(initial_size = 7356, name = "Rlr")
demo.add_population(initial_size = 8613, name = "Rlrl")
demo.add_population(initial_size = 2229, name = "Rlrr")
demo.add_population(initial_size = 8012, name = "Rlrrl")
demo.add_population(initial_size = 3484, name = "Rlrrlr")
demo.add_population(initial_size = 7543, name = "F")
demo.add_population(initial_size = 7699, name = "Rlrrr")
demo.add_population(initial_size = 7994, name = "M")
demo.add_population(initial_size = 1325, name = "Rllrlo")
demo.add_population(initial_size = 2674, name = "admix")
demo.add_population(initial_size = 5329, name = "B")
demo.add_population(initial_size = 7338, name = "Rlrrlh")
demo.add_population(initial_size = 8587, name = "admixi")
demo.add_population(initial_size = 2685, name = "Rlrc")
demo.add_population(initial_size = 6550, name = "admixc")
demo.add_population(initial_size = 3004, name = "Ro")
demo.add_population(initial_size = 6837, name = "D")
demo.add_population(initial_size = 8072, name = "admixa")
demo.add_population(initial_size = 3920, name = "Rlrlv")
demo.add_population(initial_size = 8297, name = "admixu")
demo.add_population(initial_size = 9551, name = "Rlrlh")
demo.add_population(initial_size = 3488, name = "admixx")
demo.add_population(initial_size = 1469, name = "K")
demo.add_population(initial_size = 6254, name = "Rlls")
demo.add_population(initial_size = 7968, name = "E")
demo.add_population(initial_size = 6188, name = "admixb")
demo.add_population(initial_size = 4233, name = "Rlrrry")
demo.add_population(initial_size = 4701, name = "J")
demo.add_population(initial_size = 8342, name = "admixh")
demo.add_population(initial_size = 5130, name = "L")
demo.add_population(initial_size = 8288, name = "admixio")
demo.add_population(initial_size = 4272, name = "H")
demo.add_population(initial_size = 3328, name = "admixs")
demo.add_population(initial_size = 2114, name = "A")
demo.add_population(initial_size = 2206, name = "Rlrlva")
demo.add_population(initial_size = 1306, name = "admixcd")
demo.add_population_split(time=810, ancestral="admix", derived=["B"])
demo.add_population_split(time=2430, ancestral="admixa", derived=["Rlr"])
demo.add_population_split(time=1639, ancestral="admixc", derived=["Rllrlo"])
demo.add_population_split(time=1110, ancestral="admixh", derived=["L"])
demo.add_population_split(time=810, ancestral="admixi", derived=["admixio"])
demo.add_population_split(time=390, ancestral="admixio", derived=["H"])
demo.add_population_split(time=274, ancestral="admixs", derived=["A"])
demo.add_population_split(time=1279, ancestral="admixu", derived=["Rlrrlr"])
demo.add_population_split(time=1110, ancestral="admixx", derived=["K"])
demo.add_population_split(time=2659, ancestral="R", derived=["Rl","Ro"])
demo.add_population_split(time=2509, ancestral="Rl", derived=["Rll"])
demo.add_population_split(time=2430, ancestral="Rll", derived=["Rllr","Rlls"])
demo.add_population_split(time=2336, ancestral="Rllr", derived=["Rllrl","I"])
demo.add_population_split(time=2032, ancestral="Rllrl", derived=["G"])
demo.add_population_split(time=1279, ancestral="Rllrlo", derived=["Rllrlr"])
demo.add_population_split(time=1110, ancestral="Rllrlr", derived=["C"])
demo.add_population_split(time=2336, ancestral="Rlls", derived=["E"])
demo.add_population_split(time=2336, ancestral="Rlr", derived=["Rlrl","Rlrc"])
demo.add_population_split(time=2032, ancestral="Rlrc", derived=["Rlrr"])
demo.add_population_split(time=2032, ancestral="Rlrl", derived=["Rlrlv","Rlrlh"])
demo.add_population_split(time=1766, ancestral="Rlrlv", derived=["Rlrlva"])
demo.add_population_split(time=1766, ancestral="Rlrr", derived=["Rlrrl","Rlrrr"])
demo.add_population_split(time=1639, ancestral="Rlrrl", derived=["Rlrrlh"])
demo.add_population_split(time=1110, ancestral="Rlrrlr", derived=["F"])
demo.add_population_split(time=1639, ancestral="Rlrrr", derived=["M","Rlrrry"])
demo.add_population_split(time=1279, ancestral="Rlrrry", derived=["J"])
demo.add_population_split(time=2509, ancestral="Ro", derived=["D"])
demo.add_admixture(time=810, derived="admix", ancestral=["Rllrlo","admixcd"], proportions=[0.72,0.28])
demo.add_admixture(time=2430, derived="admixa", ancestral=["Rl","Ro"], proportions=[0.29,0.71])
demo.add_admixture(time=1766, derived="admixb", ancestral=["Rlrc","Rlls"], proportions=[0.61,0.39])
demo.add_admixture(time=1639, derived="admixc", ancestral=["admixb","Rllrl"], proportions=[0.26,0.74])
demo.add_admixture(time=1279, derived="admixcd", ancestral=["Rlrlva","Rlrlh"], proportions=[0.59,0.41])
demo.add_admixture(time=1110, derived="admixh", ancestral=["Rlrrry","Rlrlv"], proportions=[0.68,0.32])
demo.add_admixture(time=810, derived="admixi", ancestral=["Rllrlr","Rlrrlh"], proportions=[0.14,0.86])
demo.add_admixture(time=274, derived="admixs", ancestral=["admixio","Rlrrlr"], proportions=[0.39,0.61])
demo.add_admixture(time=1279, derived="admixu", ancestral=["Rlrrl","Rlrlva"], proportions=[0.84,0.16])
demo.add_admixture(time=1110, derived="admixx", ancestral=["Rlrrlh","Rlrlh"], proportions=[0.76,0.24])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_3.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_3.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_3.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

