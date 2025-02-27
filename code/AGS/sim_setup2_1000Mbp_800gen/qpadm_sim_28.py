import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("I", "D", "H", "M", "F", "L", "E", "A", "C", "J", "G", "B", "K"), (232, 458, 425, 188, 349, 392, 392, 0, 425, 163, 232, 349, 188) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 2865, name = "Rl")
demo.add_population(initial_size = 2094, name = "Rlr")
demo.add_population(initial_size = 7529, name = "Rlrl")
demo.add_population(initial_size = 7457, name = "Rlrlr")
demo.add_population(initial_size = 4961, name = "I")
demo.add_population(initial_size = 7141, name = "Rlrr")
demo.add_population(initial_size = 5791, name = "D")
demo.add_population(initial_size = 3976, name = "Rr")
demo.add_population(initial_size = 1262, name = "Rrl")
demo.add_population(initial_size = 3786, name = "H")
demo.add_population(initial_size = 7768, name = "Rrlrl")
demo.add_population(initial_size = 2669, name = "M")
demo.add_population(initial_size = 6518, name = "Rrlr")
demo.add_population(initial_size = 1514, name = "F")
demo.add_population(initial_size = 6480, name = "Rrr")
demo.add_population(initial_size = 8247, name = "Rrrr")
demo.add_population(initial_size = 2362, name = "L")
demo.add_population(initial_size = 9183, name = "Ry")
demo.add_population(initial_size = 1808, name = "admix")
demo.add_population(initial_size = 5331, name = "Rrrf")
demo.add_population(initial_size = 9304, name = "E")
demo.add_population(initial_size = 2141, name = "admixx")
demo.add_population(initial_size = 9585, name = "A")
demo.add_population(initial_size = 1044, name = "Rlh")
demo.add_population(initial_size = 3196, name = "admixz")
demo.add_population(initial_size = 5516, name = "C")
demo.add_population(initial_size = 4811, name = "R")
demo.add_population(initial_size = 1726, name = "Rv")
demo.add_population(initial_size = 6968, name = "admixo")
demo.add_population(initial_size = 3150, name = "Rp")
demo.add_population(initial_size = 2064, name = "admixi")
demo.add_population(initial_size = 2123, name = "J")
demo.add_population(initial_size = 5127, name = "Rlrm")
demo.add_population(initial_size = 6629, name = "admixu")
demo.add_population(initial_size = 5599, name = "G")
demo.add_population(initial_size = 4542, name = "Rlrw")
demo.add_population(initial_size = 3930, name = "admixd")
demo.add_population(initial_size = 9775, name = "Rrrrs")
demo.add_population(initial_size = 5770, name = "B")
demo.add_population(initial_size = 1535, name = "admixv")
demo.add_population(initial_size = 4844, name = "Rrlrn")
demo.add_population(initial_size = 9888, name = "admixih")
demo.add_population(initial_size = 2602, name = "K")
demo.add_population(initial_size = 2677, name = "Rvh")
demo.add_population(initial_size = 9807, name = "admixdj")
demo.add_population_split(time=425, ancestral="admix", derived=["Rlrl"])
demo.add_population_split(time=425, ancestral="admixd", derived=["Rrlr"])
demo.add_population_split(time=188, ancestral="admixi", derived=["J"])
demo.add_population_split(time=232, ancestral="admixih", derived=["K"])
demo.add_population_split(time=349, ancestral="admixu", derived=["G"])
demo.add_population_split(time=73, ancestral="admixx", derived=["A"])
demo.add_population_split(time=458, ancestral="admixz", derived=["C"])
demo.add_population_split(time=686, ancestral="R", derived=["Rv","Rp"])
demo.add_population_split(time=564, ancestral="Rl", derived=["Rlr","Rlh"])
demo.add_population_split(time=537, ancestral="Rlr", derived=["Rlrr","Rlrw"])
demo.add_population_split(time=392, ancestral="Rlrl", derived=["Rlrlr"])
demo.add_population_split(time=349, ancestral="Rlrlr", derived=["I"])
demo.add_population_split(time=510, ancestral="Rlrr", derived=["D"])
demo.add_population_split(time=510, ancestral="Rlrw", derived=["Rlrm"])
demo.add_population_split(time=649, ancestral="Rp", derived=["Rl"])
demo.add_population_split(time=510, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=458, ancestral="Rrl", derived=["H"])
demo.add_population_split(time=392, ancestral="Rrlr", derived=["F","Rrlrn"])
demo.add_population_split(time=232, ancestral="Rrlrl", derived=["M"])
demo.add_population_split(time=349, ancestral="Rrlrn", derived=["Rrlrl"])
demo.add_population_split(time=458, ancestral="Rrr", derived=["Rrrr","Rrrf"])
demo.add_population_split(time=425, ancestral="Rrrf", derived=["E"])
demo.add_population_split(time=425, ancestral="Rrrr", derived=["L","Rrrrs"])
demo.add_population_split(time=392, ancestral="Rrrrs", derived=["B"])
demo.add_population_split(time=649, ancestral="Rv", derived=["Rvh"])
demo.add_population_split(time=564, ancestral="Rvh", derived=["Ry"])
demo.add_population_split(time=537, ancestral="Ry", derived=["Rr"])
demo.add_admixture(time=425, derived="admix", ancestral=["Rlrm","Ry"], proportions=[0.18,0.82])
demo.add_admixture(time=425, derived="admixd", ancestral=["Rrl","Rlrw"], proportions=[0.79,0.21])
demo.add_admixture(time=232, derived="admixdj", ancestral=["Rlrlr","Rvh"], proportions=[0.18,0.82])
demo.add_admixture(time=188, derived="admixi", ancestral=["Rrlrl","Rp"], proportions=[0.68,0.32])
demo.add_admixture(time=232, derived="admixih", ancestral=["Rrlrn","admixo"], proportions=[0.87,0.13])
demo.add_admixture(time=510, derived="admixo", ancestral=["Rlh","Rv"], proportions=[0.46,0.54])
demo.add_admixture(time=349, derived="admixu", ancestral=["Rlrl","Rlrm"], proportions=[0.34,0.66])
demo.add_admixture(time=349, derived="admixv", ancestral=["Rrrrs","Rrrf"], proportions=[0.89,0.11])
demo.add_admixture(time=73, derived="admixx", ancestral=["admixdj","admixv"], proportions=[0.19,0.81])
demo.add_admixture(time=458, derived="admixz", ancestral=["Rlrr","Rlh"], proportions=[0.63,0.37])
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

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_28/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_28/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_28/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

