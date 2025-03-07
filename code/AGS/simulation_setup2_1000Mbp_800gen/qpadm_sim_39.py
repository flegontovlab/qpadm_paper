import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("G", "F", "J", "D", "I", "C", "E", "L", "B", "M", "K", "A", "H"), (244, 181, 244, 476, 244, 549, 93, 93, 93, 93, 128, 244, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 9483, name = "R")
demo.add_population(initial_size = 4914, name = "Rl")
demo.add_population(initial_size = 6303, name = "Rll")
demo.add_population(initial_size = 7722, name = "Rllr")
demo.add_population(initial_size = 8185, name = "G")
demo.add_population(initial_size = 1340, name = "Rlr")
demo.add_population(initial_size = 5868, name = "Rlrl")
demo.add_population(initial_size = 7855, name = "Rlrlr")
demo.add_population(initial_size = 2493, name = "Rlrlrl")
demo.add_population(initial_size = 5249, name = "F")
demo.add_population(initial_size = 5248, name = "J")
demo.add_population(initial_size = 2377, name = "D")
demo.add_population(initial_size = 3964, name = "Rr")
demo.add_population(initial_size = 7293, name = "Rrl")
demo.add_population(initial_size = 6595, name = "Rrll")
demo.add_population(initial_size = 4496, name = "I")
demo.add_population(initial_size = 3617, name = "Rrlr")
demo.add_population(initial_size = 4237, name = "C")
demo.add_population(initial_size = 8928, name = "Rlw")
demo.add_population(initial_size = 6614, name = "admix")
demo.add_population(initial_size = 8582, name = "E")
demo.add_population(initial_size = 2886, name = "Rrlro")
demo.add_population(initial_size = 5887, name = "admixx")
demo.add_population(initial_size = 5626, name = "Rrr")
demo.add_population(initial_size = 3820, name = "admixi")
demo.add_population(initial_size = 8950, name = "L")
demo.add_population(initial_size = 6707, name = "Rlrly")
demo.add_population(initial_size = 8990, name = "admixt")
demo.add_population(initial_size = 3423, name = "B")
demo.add_population(initial_size = 8652, name = "Rrlry")
demo.add_population(initial_size = 4628, name = "admixq")
demo.add_population(initial_size = 1177, name = "M")
demo.add_population(initial_size = 2367, name = "Rlrlrla")
demo.add_population(initial_size = 1873, name = "K")
demo.add_population(initial_size = 6166, name = "admixn")
demo.add_population(initial_size = 9218, name = "Rln")
demo.add_population(initial_size = 9323, name = "admixf")
demo.add_population(initial_size = 6540, name = "Rllk")
demo.add_population(initial_size = 9878, name = "A")
demo.add_population(initial_size = 4147, name = "admixd")
demo.add_population(initial_size = 4132, name = "Rlwk")
demo.add_population(initial_size = 7238, name = "admixg")
demo.add_population(initial_size = 6640, name = "H")
demo.add_population(initial_size = 2669, name = "Rrrl")
demo.add_population(initial_size = 5322, name = "admixnb")
demo.add_population_split(time=128, ancestral="admix", derived=["E"])
demo.add_population_split(time=244, ancestral="admixd", derived=["Rrlry"])
demo.add_population_split(time=93, ancestral="admixg", derived=["H"])
demo.add_population_split(time=128, ancestral="admixi", derived=["L"])
demo.add_population_split(time=128, ancestral="admixq", derived=["M"])
demo.add_population_split(time=128, ancestral="admixt", derived=["B"])
demo.add_population_split(time=772, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=655, ancestral="Rl", derived=["Rlr","Rln"])
demo.add_population_split(time=361, ancestral="Rll", derived=["Rllr","Rllk"])
demo.add_population_split(time=298, ancestral="Rllk", derived=["A"])
demo.add_population_split(time=298, ancestral="Rllr", derived=["G"])
demo.add_population_split(time=549, ancestral="Rln", derived=["Rlw"])
demo.add_population_split(time=549, ancestral="Rlr", derived=["Rlrl","D"])
demo.add_population_split(time=476, ancestral="Rlrl", derived=["Rlrly"])
demo.add_population_split(time=298, ancestral="Rlrlr", derived=["Rlrlrl","J"])
demo.add_population_split(time=244, ancestral="Rlrlrl", derived=["F","Rlrlrla"])
demo.add_population_split(time=181, ancestral="Rlrlrla", derived=["K"])
demo.add_population_split(time=361, ancestral="Rlrly", derived=["Rlrlr"])
demo.add_population_split(time=476, ancestral="Rlw", derived=["Rll","Rlwk"])
demo.add_population_split(time=655, ancestral="Rr", derived=["C","Rrr"])
demo.add_population_split(time=361, ancestral="Rrl", derived=["Rrll","Rrlr"])
demo.add_population_split(time=298, ancestral="Rrll", derived=["I"])
demo.add_population_split(time=298, ancestral="Rrlr", derived=["Rrlro"])
demo.add_population_split(time=549, ancestral="Rrr", derived=["Rrrl"])
demo.add_population_split(time=476, ancestral="Rrrl", derived=["Rrl"])
demo.add_admixture(time=128, derived="admix", ancestral=["admixx","Rlwk"], proportions=[0.11,0.89])
demo.add_admixture(time=244, derived="admixd", ancestral=["Rrlr","Rllk"], proportions=[0.83,0.17])
demo.add_admixture(time=181, derived="admixf", ancestral=["admixnb","Rln"], proportions=[0.43,0.57])
demo.add_admixture(time=93, derived="admixg", ancestral=["admixn","Rlwk"], proportions=[0.25,0.75])
demo.add_admixture(time=128, derived="admixi", ancestral=["admixf","Rrr"], proportions=[0.19,0.81])
demo.add_admixture(time=128, derived="admixn", ancestral=["Rlrlrla","Rlrl"], proportions=[0.69,0.31])
demo.add_admixture(time=244, derived="admixnb", ancestral=["Rllr","Rrrl"], proportions=[0.43,0.57])
demo.add_admixture(time=128, derived="admixq", ancestral=["Rrlry","Rrlro"], proportions=[0.69,0.31])
demo.add_admixture(time=128, derived="admixt", ancestral=["Rrlry","Rlrly"], proportions=[0.72,0.28])
demo.add_admixture(time=181, derived="admixx", ancestral=["Rrlro","Rrll"], proportions=[0.55,0.45])
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

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_39/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_39/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_39/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

