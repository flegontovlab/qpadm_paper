import numpy
import math
import msprime
import multiprocess

indnam = ["J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("J", "C", "M", "F", "K", "G", "I", "L", "A", "B", "D", "H", "E"), (566, 566, 1954, 1294, 1294, 911, 0, 278, 0, 1166, 1376, 1166, 105) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 5441, name = "Rl")
demo.add_population(initial_size = 7290, name = "Rll")
demo.add_population(initial_size = 4539, name = "R")
demo.add_population(initial_size = 1675, name = "Rr")
demo.add_population(initial_size = 7138, name = "Rrl")
demo.add_population(initial_size = 4617, name = "Rrll")
demo.add_population(initial_size = 8759, name = "J")
demo.add_population(initial_size = 5128, name = "Rrlr")
demo.add_population(initial_size = 2473, name = "C")
demo.add_population(initial_size = 6648, name = "Rrr")
demo.add_population(initial_size = 6167, name = "M")
demo.add_population(initial_size = 3299, name = "Rrrr")
demo.add_population(initial_size = 4771, name = "Rrrrr")
demo.add_population(initial_size = 9317, name = "F")
demo.add_population(initial_size = 9220, name = "K")
demo.add_population(initial_size = 2539, name = "Rlr")
demo.add_population(initial_size = 8640, name = "admix")
demo.add_population(initial_size = 6040, name = "G")
demo.add_population(initial_size = 5776, name = "Rlld")
demo.add_population(initial_size = 7825, name = "Rrllx")
demo.add_population(initial_size = 3456, name = "Rllr")
demo.add_population(initial_size = 9496, name = "admixq")
demo.add_population(initial_size = 4756, name = "I")
demo.add_population(initial_size = 7018, name = "Rrlro")
demo.add_population(initial_size = 3708, name = "L")
demo.add_population(initial_size = 6108, name = "admixe")
demo.add_population(initial_size = 5240, name = "A")
demo.add_population(initial_size = 2138, name = "Rlb")
demo.add_population(initial_size = 6990, name = "admixx")
demo.add_population(initial_size = 3064, name = "Rrrrrn")
demo.add_population(initial_size = 4726, name = "admixqc")
demo.add_population(initial_size = 3207, name = "B")
demo.add_population(initial_size = 8361, name = "Rllry")
demo.add_population(initial_size = 7397, name = "admixf")
demo.add_population(initial_size = 6665, name = "Rrrru")
demo.add_population(initial_size = 1980, name = "D")
demo.add_population(initial_size = 5829, name = "admixy")
demo.add_population(initial_size = 6105, name = "Rllrs")
demo.add_population(initial_size = 4398, name = "H")
demo.add_population(initial_size = 7720, name = "admixl")
demo.add_population(initial_size = 8296, name = "Rllz")
demo.add_population(initial_size = 9867, name = "admixr")
demo.add_population(initial_size = 3732, name = "Rn")
demo.add_population(initial_size = 1492, name = "admixd")
demo.add_population(initial_size = 9662, name = "E")
demo.add_population_split(time=1166, ancestral="admix", derived=["G"])
demo.add_population_split(time=278, ancestral="admixd", derived=["E"])
demo.add_population_split(time=105, ancestral="admixe", derived=["A"])
demo.add_population_split(time=1166, ancestral="admixl", derived=["Rrl"])
demo.add_population_split(time=105, ancestral="admixq", derived=["I"])
demo.add_population_split(time=1294, ancestral="admixqc", derived=["B"])
demo.add_population_split(time=1376, ancestral="admixy", derived=["Rlr"])
demo.add_population_split(time=2441, ancestral="R", derived=["Rr","Rn"])
demo.add_population_split(time=2036, ancestral="Rl", derived=["Rll","Rlb"])
demo.add_population_split(time=1954, ancestral="Rll", derived=["Rllry","Rllz"])
demo.add_population_split(time=1376, ancestral="Rllr", derived=["Rllrs"])
demo.add_population_split(time=1294, ancestral="Rllrs", derived=["H"])
demo.add_population_split(time=1538, ancestral="Rllry", derived=["Rlld"])
demo.add_population_split(time=1538, ancestral="Rllz", derived=["Rllr"])
demo.add_population_split(time=2284, ancestral="Rn", derived=["Rl"])
demo.add_population_split(time=2284, ancestral="Rr", derived=["Rrr"])
demo.add_population_split(time=911, ancestral="Rrl", derived=["Rrll","Rrlr"])
demo.add_population_split(time=712, ancestral="Rrll", derived=["J","Rrllx"])
demo.add_population_split(time=712, ancestral="Rrlr", derived=["C","Rrlro"])
demo.add_population_split(time=566, ancestral="Rrlro", derived=["L"])
demo.add_population_split(time=2036, ancestral="Rrr", derived=["M","Rrrr"])
demo.add_population_split(time=1954, ancestral="Rrrr", derived=["Rrrrrn","Rrrru"])
demo.add_population_split(time=1376, ancestral="Rrrrr", derived=["F","K"])
demo.add_population_split(time=1538, ancestral="Rrrrrn", derived=["Rrrrr"])
demo.add_population_split(time=1538, ancestral="Rrrru", derived=["D"])
demo.add_admixture(time=1166, derived="admix", ancestral=["Rlr","Rlld"], proportions=[0.51,0.49])
demo.add_admixture(time=278, derived="admixd", ancestral=["Rrllx","Rn"], proportions=[0.65,0.35])
demo.add_admixture(time=105, derived="admixe", ancestral=["Rrlro","Rlr"], proportions=[0.74,0.26])
demo.add_admixture(time=1376, derived="admixf", ancestral=["admixx","Rllry"], proportions=[0.3,0.7])
demo.add_admixture(time=1166, derived="admixl", ancestral=["Rllrs","admixf"], proportions=[0.69,0.31])
demo.add_admixture(time=105, derived="admixq", ancestral=["admixr","Rllr"], proportions=[0.47,0.53])
demo.add_admixture(time=1294, derived="admixqc", ancestral=["Rlld","Rrrrrn"], proportions=[0.36,0.64])
demo.add_admixture(time=278, derived="admixr", ancestral=["Rrllx","Rllz"], proportions=[0.55,0.45])
demo.add_admixture(time=1538, derived="admixx", ancestral=["Rlb","Rr"], proportions=[0.45,0.55])
demo.add_admixture(time=1376, derived="admixy", ancestral=["Rrrru","Rlb"], proportions=[0.67,0.33])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_15.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_15.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_15.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

