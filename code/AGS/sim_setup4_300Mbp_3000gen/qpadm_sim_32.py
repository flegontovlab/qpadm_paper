import numpy
import math
import msprime
import multiprocess

indnam = ["M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("M", "J", "H", "E", "K", "I", "F", "B", "D", "C", "A", "G", "L"), (1714, 1358, 2258, 1358, 1050, 1050, 0, 0, 679, 679, 345, 679, 345) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 4560, name = "R")
demo.add_population(initial_size = 4487, name = "Rl")
demo.add_population(initial_size = 9137, name = "Rll")
demo.add_population(initial_size = 1637, name = "M")
demo.add_population(initial_size = 8164, name = "Rllr")
demo.add_population(initial_size = 2264, name = "J")
demo.add_population(initial_size = 5358, name = "H")
demo.add_population(initial_size = 9291, name = "Rr")
demo.add_population(initial_size = 4010, name = "Rrl")
demo.add_population(initial_size = 7483, name = "Rrll")
demo.add_population(initial_size = 5802, name = "Rrllr")
demo.add_population(initial_size = 1019, name = "E")
demo.add_population(initial_size = 3214, name = "Rrr")
demo.add_population(initial_size = 5710, name = "Rrrl")
demo.add_population(initial_size = 5824, name = "Rrrll")
demo.add_population(initial_size = 2317, name = "K")
demo.add_population(initial_size = 3248, name = "I")
demo.add_population(initial_size = 1565, name = "Rrrr")
demo.add_population(initial_size = 2758, name = "F")
demo.add_population(initial_size = 3231, name = "B")
demo.add_population(initial_size = 7732, name = "Rllo")
demo.add_population(initial_size = 3963, name = "admix")
demo.add_population(initial_size = 2336, name = "Rrrly")
demo.add_population(initial_size = 1538, name = "admixe")
demo.add_population(initial_size = 5967, name = "D")
demo.add_population(initial_size = 3714, name = "Rlloc")
demo.add_population(initial_size = 9216, name = "admixx")
demo.add_population(initial_size = 3665, name = "Rrrlye")
demo.add_population(initial_size = 6383, name = "admixy")
demo.add_population(initial_size = 2548, name = "Rrt")
demo.add_population(initial_size = 8026, name = "admixt")
demo.add_population(initial_size = 8876, name = "C")
demo.add_population(initial_size = 9468, name = "admixex")
demo.add_population(initial_size = 9720, name = "admixc")
demo.add_population(initial_size = 1413, name = "Rllrx")
demo.add_population(initial_size = 3950, name = "admixyp")
demo.add_population(initial_size = 5203, name = "A")
demo.add_population(initial_size = 8938, name = "Rrru")
demo.add_population(initial_size = 9244, name = "admixj")
demo.add_population(initial_size = 2839, name = "G")
demo.add_population(initial_size = 7245, name = "admixce")
demo.add_population(initial_size = 1921, name = "L")
demo.add_population(initial_size = 1348, name = "admixw")
demo.add_population(initial_size = 5997, name = "Rn")
demo.add_population(initial_size = 1658, name = "admixxw")
demo.add_population_split(time=1358, ancestral="admix", derived=["admixex"])
demo.add_population_split(time=904, ancestral="admixc", derived=["admixce"])
demo.add_population_split(time=679, ancestral="admixce", derived=["L"])
demo.add_population_split(time=904, ancestral="admixe", derived=["D"])
demo.add_population_split(time=904, ancestral="admixj", derived=["G"])
demo.add_population_split(time=904, ancestral="admixt", derived=["C"])
demo.add_population_split(time=345, ancestral="admixw", derived=["Rrrr"])
demo.add_population_split(time=2258, ancestral="admixxw", derived=["Rll"])
demo.add_population_split(time=679, ancestral="admixyp", derived=["A"])
demo.add_population_split(time=2655, ancestral="R", derived=["Rl","Rn"])
demo.add_population_split(time=2482, ancestral="Rl", derived=["H"])
demo.add_population_split(time=2081, ancestral="Rll", derived=["M","Rllo"])
demo.add_population_split(time=1714, ancestral="Rllo", derived=["Rllr","Rlloc"])
demo.add_population_split(time=1616, ancestral="Rllr", derived=["J","Rllrx"])
demo.add_population_split(time=2482, ancestral="Rn", derived=["Rr"])
demo.add_population_split(time=2258, ancestral="Rr", derived=["Rrl","Rrt"])
demo.add_population_split(time=2081, ancestral="Rrl", derived=["Rrll"])
demo.add_population_split(time=1714, ancestral="Rrll", derived=["Rrllr"])
demo.add_population_split(time=1616, ancestral="Rrllr", derived=["E"])
demo.add_population_split(time=1714, ancestral="Rrr", derived=["Rrrl","Rrru"])
demo.add_population_split(time=1616, ancestral="Rrrl", derived=["Rrrll","Rrrly"])
demo.add_population_split(time=1358, ancestral="Rrrll", derived=["K","I"])
demo.add_population_split(time=1358, ancestral="Rrrly", derived=["Rrrlye"])
demo.add_population_split(time=75, ancestral="Rrrr", derived=["F","B"])
demo.add_population_split(time=2081, ancestral="Rrt", derived=["Rrr"])
demo.add_admixture(time=1358, derived="admix", ancestral=["Rlloc","Rrl"], proportions=[0.44,0.56])
demo.add_admixture(time=904, derived="admixc", ancestral=["admixex","Rllrx"], proportions=[0.39,0.61])
demo.add_admixture(time=904, derived="admixe", ancestral=["admixx","Rrllr"], proportions=[0.37,0.63])
demo.add_admixture(time=904, derived="admixj", ancestral=["Rrrlye","Rrru"], proportions=[0.29,0.71])
demo.add_admixture(time=904, derived="admixt", ancestral=["admixex","Rrt"], proportions=[0.22,0.78])
demo.add_admixture(time=345, derived="admixw", ancestral=["admixce","Rrru"], proportions=[0.5,0.5])
demo.add_admixture(time=1050, derived="admixx", ancestral=["Rrrly","Rlloc"], proportions=[0.74,0.26])
demo.add_admixture(time=2258, derived="admixxw", ancestral=["Rl","Rn"], proportions=[0.48,0.52])
demo.add_admixture(time=904, derived="admixy", ancestral=["Rrrlye","Rrll"], proportions=[0.56,0.44])
demo.add_admixture(time=679, derived="admixyp", ancestral=["admixy","Rllrx"], proportions=[0.2,0.8])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_32.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_32.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_32.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

