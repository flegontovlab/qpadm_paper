import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("C", "I", "K", "G", "L", "A", "E", "H", "D", "J", "F", "B", "M"), (53, 0, 0, 217, 295, 295, 83, 106, 106, 0, 0, 53, 53) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 9141, name = "R")
demo.add_population(initial_size = 2570, name = "Rl")
demo.add_population(initial_size = 1041, name = "Rll")
demo.add_population(initial_size = 5916, name = "Rlll")
demo.add_population(initial_size = 4949, name = "Rllrl")
demo.add_population(initial_size = 4338, name = "C")
demo.add_population(initial_size = 8473, name = "Rllrlr")
demo.add_population(initial_size = 8732, name = "I")
demo.add_population(initial_size = 5225, name = "K")
demo.add_population(initial_size = 4423, name = "Rllr")
demo.add_population(initial_size = 8568, name = "G")
demo.add_population(initial_size = 7368, name = "Rlr")
demo.add_population(initial_size = 1188, name = "L")
demo.add_population(initial_size = 1631, name = "Rr")
demo.add_population(initial_size = 7610, name = "A")
demo.add_population(initial_size = 6633, name = "Rrrl")
demo.add_population(initial_size = 7433, name = "E")
demo.add_population(initial_size = 3375, name = "Rlry")
demo.add_population(initial_size = 3859, name = "Rlrr")
demo.add_population(initial_size = 1220, name = "admix")
demo.add_population(initial_size = 9765, name = "Rlllk")
demo.add_population(initial_size = 1989, name = "admixb")
demo.add_population(initial_size = 8463, name = "Rlllkt")
demo.add_population(initial_size = 5459, name = "H")
demo.add_population(initial_size = 9572, name = "admixl")
demo.add_population(initial_size = 7610, name = "Rrm")
demo.add_population(initial_size = 5404, name = "Rrr")
demo.add_population(initial_size = 1812, name = "admixo")
demo.add_population(initial_size = 6383, name = "Rrrm")
demo.add_population(initial_size = 5636, name = "D")
demo.add_population(initial_size = 7507, name = "admixy")
demo.add_population(initial_size = 8495, name = "J")
demo.add_population(initial_size = 8265, name = "Rrrmj")
demo.add_population(initial_size = 7359, name = "admixj")
demo.add_population(initial_size = 8204, name = "F")
demo.add_population(initial_size = 5255, name = "Rlli")
demo.add_population(initial_size = 8385, name = "admixyd")
demo.add_population(initial_size = 1887, name = "Rlllkm")
demo.add_population(initial_size = 3988, name = "admixe")
demo.add_population(initial_size = 8337, name = "Ro")
demo.add_population(initial_size = 1783, name = "admixr")
demo.add_population(initial_size = 2376, name = "B")
demo.add_population(initial_size = 1135, name = "Rrrw")
demo.add_population(initial_size = 8760, name = "admixh")
demo.add_population(initial_size = 5270, name = "M")
demo.add_population_split(time=197, ancestral="admixe", derived=["Rrrl"])
demo.add_population_split(time=83, ancestral="admixh", derived=["M"])
demo.add_population_split(time=53, ancestral="admixj", derived=["F"])
demo.add_population_split(time=106, ancestral="admixl", derived=["Rllrl"])
demo.add_population_split(time=217, ancestral="admixo", derived=["Rlry"])
demo.add_population_split(time=83, ancestral="admixr", derived=["B"])
demo.add_population_split(time=53, ancestral="admixy", derived=["J"])
demo.add_population_split(time=510, ancestral="R", derived=["Rl","Ro"])
demo.add_population_split(time=487, ancestral="Rl", derived=["Rll","Rlr"])
demo.add_population_split(time=376, ancestral="Rll", derived=["Rlll","Rlli"])
demo.add_population_split(time=295, ancestral="Rlli", derived=["Rllr"])
demo.add_population_split(time=295, ancestral="Rlll", derived=["Rlllk"])
demo.add_population_split(time=249, ancestral="Rlllk", derived=["Rlllkm"])
demo.add_population_split(time=217, ancestral="Rlllkm", derived=["Rlllkt"])
demo.add_population_split(time=197, ancestral="Rlllkt", derived=["H"])
demo.add_population_split(time=249, ancestral="Rllr", derived=["G"])
demo.add_population_split(time=83, ancestral="Rllrl", derived=["C","Rllrlr"])
demo.add_population_split(time=53, ancestral="Rllrlr", derived=["I","K"])
demo.add_population_split(time=376, ancestral="Rlr", derived=["L"])
demo.add_population_split(time=197, ancestral="Rlry", derived=["Rlrr"])
demo.add_population_split(time=487, ancestral="Ro", derived=["Rr"])
demo.add_population_split(time=376, ancestral="Rr", derived=["A","Rrm"])
demo.add_population_split(time=295, ancestral="Rrm", derived=["Rrr"])
demo.add_population_split(time=249, ancestral="Rrr", derived=["Rrrw"])
demo.add_population_split(time=106, ancestral="Rrrl", derived=["E"])
demo.add_population_split(time=197, ancestral="Rrrm", derived=["D","Rrrmj"])
demo.add_population_split(time=217, ancestral="Rrrw", derived=["Rrrm"])
demo.add_admixture(time=106, derived="admix", ancestral=["Rlry","Rlll"], proportions=[0.66,0.34])
demo.add_admixture(time=217, derived="admixb", ancestral=["Rllr","Rlllk"], proportions=[0.8,0.2])
demo.add_admixture(time=197, derived="admixe", ancestral=["Rlllkm","Rrr"], proportions=[0.42,0.58])
demo.add_admixture(time=83, derived="admixh", ancestral=["Rlrr","Rrrw"], proportions=[0.49,0.51])
demo.add_admixture(time=53, derived="admixj", ancestral=["Rlrr","Rrrmj"], proportions=[0.45,0.55])
demo.add_admixture(time=106, derived="admixl", ancestral=["Rlllkt","admixb"], proportions=[0.77,0.23])
demo.add_admixture(time=217, derived="admixo", ancestral=["admixyd","Rrm"], proportions=[0.1,0.9])
demo.add_admixture(time=83, derived="admixr", ancestral=["admix","Ro"], proportions=[0.43,0.57])
demo.add_admixture(time=53, derived="admixy", ancestral=["Rrrl","Rrrmj"], proportions=[0.24,0.76])
demo.add_admixture(time=249, derived="admixyd", ancestral=["Rlli","Rlr"], proportions=[0.47,0.53])
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

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_37/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_37/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_37/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

