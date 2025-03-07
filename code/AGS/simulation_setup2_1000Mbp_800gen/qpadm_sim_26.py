import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("F", "K", "C", "D", "M", "I", "G", "B", "J", "H", "L", "E", "A"), (366, 197, 576, 509, 509, 509, 0, 509, 273, 392, 219, 37, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 1316, name = "R")
demo.add_population(initial_size = 9107, name = "Rl")
demo.add_population(initial_size = 8076, name = "Rll")
demo.add_population(initial_size = 9946, name = "Rlll")
demo.add_population(initial_size = 9424, name = "Rllll")
demo.add_population(initial_size = 5750, name = "F")
demo.add_population(initial_size = 2333, name = "Rlllrr")
demo.add_population(initial_size = 5979, name = "K")
demo.add_population(initial_size = 2964, name = "Rlr")
demo.add_population(initial_size = 7172, name = "C")
demo.add_population(initial_size = 1786, name = "Rlrr")
demo.add_population(initial_size = 9401, name = "D")
demo.add_population(initial_size = 1862, name = "M")
demo.add_population(initial_size = 8274, name = "Rr")
demo.add_population(initial_size = 1339, name = "Rrr")
demo.add_population(initial_size = 5971, name = "Rrrl")
demo.add_population(initial_size = 2531, name = "I")
demo.add_population(initial_size = 9314, name = "admix")
demo.add_population(initial_size = 3883, name = "Rlllg")
demo.add_population(initial_size = 6512, name = "admixx")
demo.add_population(initial_size = 8564, name = "G")
demo.add_population(initial_size = 9635, name = "Rre")
demo.add_population(initial_size = 6750, name = "B")
demo.add_population(initial_size = 6002, name = "Rllf")
demo.add_population(initial_size = 1980, name = "admixa")
demo.add_population(initial_size = 8582, name = "Rrh")
demo.add_population(initial_size = 8627, name = "admixv")
demo.add_population(initial_size = 7880, name = "Rllm")
demo.add_population(initial_size = 7856, name = "admixj")
demo.add_population(initial_size = 7823, name = "J")
demo.add_population(initial_size = 1744, name = "Rllfn")
demo.add_population(initial_size = 3862, name = "admixvl")
demo.add_population(initial_size = 1544, name = "Rrrlz")
demo.add_population(initial_size = 1818, name = "H")
demo.add_population(initial_size = 7511, name = "admixjv")
demo.add_population(initial_size = 6901, name = "Rlllr")
demo.add_population(initial_size = 9763, name = "Rlllri")
demo.add_population(initial_size = 8681, name = "admixe")
demo.add_population(initial_size = 6140, name = "Rlllle")
demo.add_population(initial_size = 4277, name = "admixf")
demo.add_population(initial_size = 9453, name = "L")
demo.add_population(initial_size = 9921, name = "admixw")
demo.add_population(initial_size = 5238, name = "E")
demo.add_population(initial_size = 7545, name = "admixp")
demo.add_population(initial_size = 6216, name = "A")
demo.add_population_split(time=197, ancestral="admix", derived=["admixw"])
demo.add_population_split(time=273, ancestral="admixf", derived=["L"])
demo.add_population_split(time=366, ancestral="admixj", derived=["J"])
demo.add_population_split(time=392, ancestral="admixjv", derived=["Rlllr"])
demo.add_population_split(time=37, ancestral="admixp", derived=["A"])
demo.add_population_split(time=509, ancestral="admixv", derived=["Rllll"])
demo.add_population_split(time=96, ancestral="admixw", derived=["E"])
demo.add_population_split(time=37, ancestral="admixx", derived=["G"])
demo.add_population_split(time=781, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=702, ancestral="Rl", derived=["Rll","Rlr"])
demo.add_population_split(time=632, ancestral="Rll", derived=["Rlll","Rllm"])
demo.add_population_split(time=509, ancestral="Rllf", derived=["Rllfn"])
demo.add_population_split(time=576, ancestral="Rlll", derived=["Rlllg"])
demo.add_population_split(time=392, ancestral="Rllll", derived=["F","Rlllle"])
demo.add_population_split(time=366, ancestral="Rlllr", derived=["Rlllri"])
demo.add_population_split(time=273, ancestral="Rlllri", derived=["Rlllrr"])
demo.add_population_split(time=219, ancestral="Rlllrr", derived=["K"])
demo.add_population_split(time=576, ancestral="Rllm", derived=["Rllf"])
demo.add_population_split(time=632, ancestral="Rlr", derived=["C","Rlrr"])
demo.add_population_split(time=576, ancestral="Rlrr", derived=["D","M"])
demo.add_population_split(time=702, ancestral="Rr", derived=["Rrr","Rrh"])
demo.add_population_split(time=576, ancestral="Rre", derived=["B"])
demo.add_population_split(time=632, ancestral="Rrh", derived=["Rre"])
demo.add_population_split(time=632, ancestral="Rrr", derived=["Rrrl"])
demo.add_population_split(time=576, ancestral="Rrrl", derived=["I","Rrrlz"])
demo.add_population_split(time=509, ancestral="Rrrlz", derived=["H"])
demo.add_admixture(time=197, derived="admix", ancestral=["Rlllrr","Rllfn"], proportions=[0.53,0.47])
demo.add_admixture(time=392, derived="admixa", ancestral=["Rllf","Rre"], proportions=[0.41,0.59])
demo.add_admixture(time=219, derived="admixe", ancestral=["Rlllri","Rlllg"], proportions=[0.8,0.2])
demo.add_admixture(time=273, derived="admixf", ancestral=["Rlllr","Rlllle"], proportions=[0.87,0.13])
demo.add_admixture(time=366, derived="admixj", ancestral=["admixa","Rllm"], proportions=[0.29,0.71])
demo.add_admixture(time=392, derived="admixjv", ancestral=["Rlllg","Rrrlz"], proportions=[0.48,0.52])
demo.add_admixture(time=37, derived="admixp", ancestral=["admixw","Rrr"], proportions=[0.13,0.87])
demo.add_admixture(time=509, derived="admixv", ancestral=["Rlll","Rrh"], proportions=[0.38,0.62])
demo.add_admixture(time=197, derived="admixvl", ancestral=["admixe","Rllfn"], proportions=[0.37,0.63])
demo.add_admixture(time=37, derived="admixx", ancestral=["admixvl","Rlllle"], proportions=[0.46,0.54])
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

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_26/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_26/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_26/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

