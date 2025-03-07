import numpy
import math
import msprime
import multiprocess

indnam = ["D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("D", "I", "M", "K", "L", "A", "F", "E", "B", "C", "J", "H", "G"), (1031, 1200, 1646, 2010, 364, 1200, 0, 664, 848, 664, 848, 664, 364) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 7362, name = "Rl")
demo.add_population(initial_size = 3016, name = "Rll")
demo.add_population(initial_size = 7725, name = "Rlll")
demo.add_population(initial_size = 7022, name = "D")
demo.add_population(initial_size = 3698, name = "I")
demo.add_population(initial_size = 5499, name = "M")
demo.add_population(initial_size = 9222, name = "R")
demo.add_population(initial_size = 1414, name = "Rr")
demo.add_population(initial_size = 9203, name = "Rrl")
demo.add_population(initial_size = 3846, name = "Rrll")
demo.add_population(initial_size = 6877, name = "Rrr")
demo.add_population(initial_size = 1745, name = "K")
demo.add_population(initial_size = 2779, name = "Rrrr")
demo.add_population(initial_size = 6911, name = "Rrrrl")
demo.add_population(initial_size = 3574, name = "Rrrrlr")
demo.add_population(initial_size = 7447, name = "Rrrrr")
demo.add_population(initial_size = 4543, name = "Rrrrlrc")
demo.add_population(initial_size = 3316, name = "L")
demo.add_population(initial_size = 5251, name = "admix")
demo.add_population(initial_size = 8485, name = "Rrllv")
demo.add_population(initial_size = 2804, name = "admixj")
demo.add_population(initial_size = 4838, name = "Rri")
demo.add_population(initial_size = 4355, name = "admixn")
demo.add_population(initial_size = 5191, name = "A")
demo.add_population(initial_size = 5309, name = "Rh")
demo.add_population(initial_size = 9066, name = "admixl")
demo.add_population(initial_size = 7757, name = "F")
demo.add_population(initial_size = 2937, name = "Ry")
demo.add_population(initial_size = 1147, name = "admixln")
demo.add_population(initial_size = 8949, name = "E")
demo.add_population(initial_size = 1053, name = "Ryt")
demo.add_population(initial_size = 5961, name = "admixw")
demo.add_population(initial_size = 1333, name = "B")
demo.add_population(initial_size = 9689, name = "Rrrrlrv")
demo.add_population(initial_size = 8791, name = "admixg")
demo.add_population(initial_size = 6329, name = "C")
demo.add_population(initial_size = 6411, name = "admixjk")
demo.add_population(initial_size = 4083, name = "J")
demo.add_population(initial_size = 2299, name = "admixi")
demo.add_population(initial_size = 1806, name = "H")
demo.add_population(initial_size = 5600, name = "Rrrrlrvh")
demo.add_population(initial_size = 6816, name = "admixx")
demo.add_population(initial_size = 2950, name = "G")
demo.add_population(initial_size = 1689, name = "Rre")
demo.add_population(initial_size = 6524, name = "admixa")
demo.add_population_split(time=848, ancestral="admixg", derived=["C"])
demo.add_population_split(time=848, ancestral="admixi", derived=["H"])
demo.add_population_split(time=1200, ancestral="admixj", derived=["admixjk"])
demo.add_population_split(time=1031, ancestral="admixjk", derived=["J"])
demo.add_population_split(time=199, ancestral="admixl", derived=["F"])
demo.add_population_split(time=848, ancestral="admixln", derived=["E"])
demo.add_population_split(time=1646, ancestral="admixn", derived=["A"])
demo.add_population_split(time=1031, ancestral="admixw", derived=["B"])
demo.add_population_split(time=664, ancestral="admixx", derived=["G"])
demo.add_population_split(time=2779, ancestral="R", derived=["Rr","Ry"])
demo.add_population_split(time=2254, ancestral="Rh", derived=["Rl"])
demo.add_population_split(time=2010, ancestral="Rl", derived=["Rll","M"])
demo.add_population_split(time=1646, ancestral="Rll", derived=["Rlll","I"])
demo.add_population_split(time=1200, ancestral="Rlll", derived=["D"])
demo.add_population_split(time=2704, ancestral="Rr", derived=["Rri","Rre"])
demo.add_population_split(time=2344, ancestral="Rre", derived=["Rrr"])
demo.add_population_split(time=2344, ancestral="Rri", derived=["Rrl"])
demo.add_population_split(time=2254, ancestral="Rrl", derived=["Rrll"])
demo.add_population_split(time=2010, ancestral="Rrll", derived=["Rrllv"])
demo.add_population_split(time=2254, ancestral="Rrr", derived=["K","Rrrr"])
demo.add_population_split(time=2010, ancestral="Rrrr", derived=["Rrrrl","Rrrrr"])
demo.add_population_split(time=1646, ancestral="Rrrrl", derived=["Rrrrlr"])
demo.add_population_split(time=1200, ancestral="Rrrrlr", derived=["Rrrrlrv"])
demo.add_population_split(time=664, ancestral="Rrrrlrc", derived=["L"])
demo.add_population_split(time=1031, ancestral="Rrrrlrv", derived=["Rrrrlrvh"])
demo.add_population_split(time=848, ancestral="Rrrrlrvh", derived=["Rrrrlrc"])
demo.add_population_split(time=2704, ancestral="Ry", derived=["Ryt"])
demo.add_population_split(time=2344, ancestral="Ryt", derived=["Rh"])
demo.add_admixture(time=364, derived="admix", ancestral=["Rrrrlrc","Rrrrr"], proportions=[0.19,0.81])
demo.add_admixture(time=2254, derived="admixa", ancestral=["Rre","Ry"], proportions=[0.26,0.74])
demo.add_admixture(time=848, derived="admixg", ancestral=["Rrrrlrv","Rrrrr"], proportions=[0.11,0.89])
demo.add_admixture(time=848, derived="admixi", ancestral=["admixjk","Rrllv"], proportions=[0.26,0.74])
demo.add_admixture(time=1200, derived="admixj", ancestral=["Rrrrl","Rrllv"], proportions=[0.76,0.24])
demo.add_admixture(time=199, derived="admixl", ancestral=["admix","Rh"], proportions=[0.24,0.76])
demo.add_admixture(time=848, derived="admixln", ancestral=["Rrrrlr","admixa"], proportions=[0.55,0.45])
demo.add_admixture(time=1646, derived="admixn", ancestral=["Rrll","Rri"], proportions=[0.78,0.22])
demo.add_admixture(time=1031, derived="admixw", ancestral=["Rlll","Ryt"], proportions=[0.13,0.87])
demo.add_admixture(time=664, derived="admixx", ancestral=["Rrrrlrvh","Rrl"], proportions=[0.56,0.44])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_35.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_35.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_35.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

