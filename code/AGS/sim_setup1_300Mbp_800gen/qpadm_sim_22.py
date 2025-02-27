import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("F", "L", "I", "G", "B", "M", "H", "D", "J", "A", "E", "C", "K"), (665, 372, 465, 552, 491, 223, 0, 223, 223, 119, 372, 223, 372) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 8255, name = "R")
demo.add_population(initial_size = 8175, name = "Rl")
demo.add_population(initial_size = 2401, name = "F")
demo.add_population(initial_size = 5675, name = "Rlr")
demo.add_population(initial_size = 7412, name = "Rlrr")
demo.add_population(initial_size = 3899, name = "Rlrrl")
demo.add_population(initial_size = 6001, name = "L")
demo.add_population(initial_size = 7719, name = "Rlrrlr")
demo.add_population(initial_size = 6044, name = "I")
demo.add_population(initial_size = 8653, name = "Rr")
demo.add_population(initial_size = 9873, name = "Rrl")
demo.add_population(initial_size = 4738, name = "G")
demo.add_population(initial_size = 7639, name = "Rrr")
demo.add_population(initial_size = 1416, name = "Rrrr")
demo.add_population(initial_size = 3697, name = "Rrrrr")
demo.add_population(initial_size = 7669, name = "B")
demo.add_population(initial_size = 2987, name = "Rlrrlrw")
demo.add_population(initial_size = 2797, name = "M")
demo.add_population(initial_size = 9091, name = "admix")
demo.add_population(initial_size = 3121, name = "H")
demo.add_population(initial_size = 4890, name = "Rrrn")
demo.add_population(initial_size = 3635, name = "admixc")
demo.add_population(initial_size = 8820, name = "D")
demo.add_population(initial_size = 8474, name = "Rlre")
demo.add_population(initial_size = 2646, name = "admixj")
demo.add_population(initial_size = 8399, name = "Rlrex")
demo.add_population(initial_size = 5484, name = "admixa")
demo.add_population(initial_size = 5345, name = "J")
demo.add_population(initial_size = 8589, name = "Rlrexp")
demo.add_population(initial_size = 2529, name = "admixt")
demo.add_population(initial_size = 7166, name = "A")
demo.add_population(initial_size = 7530, name = "Rlt")
demo.add_population(initial_size = 4294, name = "Rrrl")
demo.add_population(initial_size = 6085, name = "admixy")
demo.add_population(initial_size = 7460, name = "Rlc")
demo.add_population(initial_size = 5186, name = "admixb")
demo.add_population(initial_size = 3937, name = "E")
demo.add_population(initial_size = 6887, name = "Rrrno")
demo.add_population(initial_size = 8577, name = "admixh")
demo.add_population(initial_size = 5219, name = "admixyp")
demo.add_population(initial_size = 4942, name = "admixl")
demo.add_population(initial_size = 8192, name = "C")
demo.add_population(initial_size = 2608, name = "Rltn")
demo.add_population(initial_size = 5940, name = "admixx")
demo.add_population(initial_size = 1472, name = "K")
demo.add_population_split(time=27, ancestral="admix", derived=["H"])
demo.add_population_split(time=321, ancestral="admixa", derived=["J"])
demo.add_population_split(time=465, ancestral="admixb", derived=["E"])
demo.add_population_split(time=321, ancestral="admixc", derived=["D"])
demo.add_population_split(time=321, ancestral="admixl", derived=["C"])
demo.add_population_split(time=223, ancestral="admixt", derived=["A"])
demo.add_population_split(time=465, ancestral="admixx", derived=["K"])
demo.add_population_split(time=465, ancestral="admixy", derived=["admixyp"])
demo.add_population_split(time=783, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=693, ancestral="Rl", derived=["F","Rlc"])
demo.add_population_split(time=665, ancestral="Rlc", derived=["Rlt"])
demo.add_population_split(time=519, ancestral="Rlr", derived=["Rlrr","Rlre"])
demo.add_population_split(time=491, ancestral="Rlre", derived=["Rlrex"])
demo.add_population_split(time=465, ancestral="Rlrex", derived=["Rlrexp"])
demo.add_population_split(time=491, ancestral="Rlrr", derived=["Rlrrl","I"])
demo.add_population_split(time=465, ancestral="Rlrrl", derived=["L","Rlrrlr"])
demo.add_population_split(time=372, ancestral="Rlrrlr", derived=["Rlrrlrw"])
demo.add_population_split(time=321, ancestral="Rlrrlrw", derived=["M"])
demo.add_population_split(time=552, ancestral="Rlt", derived=["Rlr","Rltn"])
demo.add_population_split(time=693, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=665, ancestral="Rrl", derived=["G"])
demo.add_population_split(time=665, ancestral="Rrr", derived=["Rrrr","Rrrn"])
demo.add_population_split(time=552, ancestral="Rrrn", derived=["Rrrno"])
demo.add_population_split(time=519, ancestral="Rrrno", derived=["Rrrl"])
demo.add_population_split(time=552, ancestral="Rrrr", derived=["Rrrrr"])
demo.add_population_split(time=519, ancestral="Rrrrr", derived=["B"])
demo.add_admixture(time=27, derived="admix", ancestral=["admixj","admixyp"], proportions=[0.43,0.57])
demo.add_admixture(time=321, derived="admixa", ancestral=["Rlrexp","Rrrr"], proportions=[0.45,0.55])
demo.add_admixture(time=465, derived="admixb", ancestral=["Rrrl","Rlc"], proportions=[0.81,0.19])
demo.add_admixture(time=321, derived="admixc", ancestral=["Rlrrlr","Rrrn"], proportions=[0.26,0.74])
demo.add_admixture(time=321, derived="admixh", ancestral=["Rlrexp","Rrrno"], proportions=[0.46,0.54])
demo.add_admixture(time=223, derived="admixj", ancestral=["Rlrrlrw","Rlrex"], proportions=[0.54,0.46])
demo.add_admixture(time=321, derived="admixl", ancestral=["admixyp","Rrrrr"], proportions=[0.44,0.56])
demo.add_admixture(time=223, derived="admixt", ancestral=["admixh","Rrl"], proportions=[0.34,0.66])
demo.add_admixture(time=465, derived="admixx", ancestral=["Rlre","Rltn"], proportions=[0.43,0.57])
demo.add_admixture(time=465, derived="admixy", ancestral=["Rrrl","Rltn"], proportions=[0.72,0.28])
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

numpy.savetxt('./projects/qpadm/results/simout/3chr/sim_22/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/3chr/sim_22/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/3chr/sim_22/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

