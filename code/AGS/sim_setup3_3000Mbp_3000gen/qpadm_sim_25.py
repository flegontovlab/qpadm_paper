import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("L", "G", "F", "D", "K", "C", "B", "E", "M", "I", "J", "A", "H"), (870, 1088, 1001, 1406, 870, 870, 161, 0, 1088, 510, 1088, 870, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 3201, name = "R")
demo.add_population(initial_size = 7525, name = "Rl")
demo.add_population(initial_size = 9401, name = "Rll")
demo.add_population(initial_size = 4664, name = "Rlll")
demo.add_population(initial_size = 2540, name = "Rlllr")
demo.add_population(initial_size = 4276, name = "Rlllrl")
demo.add_population(initial_size = 1682, name = "Rlllrlr")
demo.add_population(initial_size = 1430, name = "L")
demo.add_population(initial_size = 9974, name = "G")
demo.add_population(initial_size = 8893, name = "Rllr")
demo.add_population(initial_size = 3481, name = "Rllrl")
demo.add_population(initial_size = 5617, name = "F")
demo.add_population(initial_size = 5712, name = "D")
demo.add_population(initial_size = 2447, name = "Rlr")
demo.add_population(initial_size = 5470, name = "Rlrl")
demo.add_population(initial_size = 1877, name = "Rlrr")
demo.add_population(initial_size = 3469, name = "K")
demo.add_population(initial_size = 4200, name = "C")
demo.add_population(initial_size = 9107, name = "Rlz")
demo.add_population(initial_size = 7739, name = "admix")
demo.add_population(initial_size = 5714, name = "B")
demo.add_population(initial_size = 9083, name = "Rle")
demo.add_population(initial_size = 9033, name = "admixv")
demo.add_population(initial_size = 1141, name = "E")
demo.add_population(initial_size = 3939, name = "Rlzq")
demo.add_population(initial_size = 2407, name = "admixq")
demo.add_population(initial_size = 3990, name = "Rlrlk")
demo.add_population(initial_size = 1252, name = "M")
demo.add_population(initial_size = 1950, name = "admixk")
demo.add_population(initial_size = 2583, name = "Rlllrlrw")
demo.add_population(initial_size = 3391, name = "admixw")
demo.add_population(initial_size = 5968, name = "Rlem")
demo.add_population(initial_size = 2551, name = "admixwj")
demo.add_population(initial_size = 6330, name = "I")
demo.add_population(initial_size = 8656, name = "Rlemz")
demo.add_population(initial_size = 7996, name = "admixt")
demo.add_population(initial_size = 7609, name = "Rlemzx")
demo.add_population(initial_size = 3031, name = "admixb")
demo.add_population(initial_size = 5469, name = "J")
demo.add_population(initial_size = 4133, name = "Rllrli")
demo.add_population(initial_size = 5761, name = "A")
demo.add_population(initial_size = 9918, name = "admixqu")
demo.add_population(initial_size = 6793, name = "H")
demo.add_population(initial_size = 8096, name = "Rlllz")
demo.add_population(initial_size = 9087, name = "admixqm")
demo.add_population_split(time=2835, ancestral="R", derived=["Rl"])
demo.add_population_split(time=2689, ancestral="Rl", derived=["Rlz","Rle"])
demo.add_population_split(time=2606, ancestral="Rle", derived=["Rlem"])
demo.add_population_split(time=2168, ancestral="Rlem", derived=["Rll","Rlemz"])
demo.add_population_split(time=2074, ancestral="Rlemz", derived=["Rlemzx"])
demo.add_population_split(time=2074, ancestral="Rll", derived=["Rlll","Rllr"])
demo.add_population_split(time=1755, ancestral="Rlll", derived=["Rlllr","Rlllz"])
demo.add_population_split(time=1406, ancestral="Rlllr", derived=["Rlllrl","G"])
demo.add_population_split(time=1088, ancestral="Rlllrl", derived=["Rlllrlr","Rlllrlrw"])
demo.add_population_split(time=1001, ancestral="Rlllrlr", derived=["L"])
demo.add_population_split(time=1755, ancestral="Rllr", derived=["D"])
demo.add_population_split(time=1088, ancestral="Rllrl", derived=["F","Rllrli"])
demo.add_population_split(time=1001, ancestral="Rllrli", derived=["A"])
demo.add_population_split(time=2074, ancestral="Rlr", derived=["Rlrl"])
demo.add_population_split(time=1755, ancestral="Rlrl", derived=["Rlrlk"])
demo.add_population_split(time=1406, ancestral="Rlrlk", derived=["M"])
demo.add_population_split(time=1001, ancestral="Rlrr", derived=["K","C"])
demo.add_population_split(time=2606, ancestral="Rlz", derived=["Rlzq"])
demo.add_population_split(time=2168, ancestral="Rlzq", derived=["Rlr"])
demo.add_population_split(time=341, ancestral="admix", derived=["B"])
demo.add_population_split(time=1406, ancestral="admixb", derived=["J"])
demo.add_population_split(time=1406, ancestral="admixq", derived=["Rllrl"])
demo.add_population_split(time=1088, ancestral="admixqm", derived=["Rlrr"])
demo.add_population_split(time=161, ancestral="admixqu", derived=["H"])
demo.add_population_split(time=161, ancestral="admixv", derived=["E"])
demo.add_population_split(time=870, ancestral="admixwj", derived=["I"])
demo.add_admixture(time=341, derived="admix", ancestral=["admixk","Rlllz"], proportions=[0.27,0.73])
demo.add_admixture(time=1406, derived="admixb", ancestral=["Rlemzx","R"], proportions=[0.58,0.42])
demo.add_admixture(time=1088, derived="admixk", ancestral=["Rlrlk","Rlz"], proportions=[0.44,0.56])
demo.add_admixture(time=1406, derived="admixq", ancestral=["Rllr","Rlzq"], proportions=[0.49,0.51])
demo.add_admixture(time=1088, derived="admixqm", ancestral=["Rlllz","Rlr"], proportions=[0.34,0.66])
demo.add_admixture(time=161, derived="admixqu", ancestral=["admixt","Rllrli"], proportions=[0.43,0.57])
demo.add_admixture(time=870, derived="admixt", ancestral=["Rlllrlr","Rlemzx"], proportions=[0.85,0.15])
demo.add_admixture(time=161, derived="admixv", ancestral=["admixw","Rlrl"], proportions=[0.3,0.7])
demo.add_admixture(time=870, derived="admixw", ancestral=["Rlllrlrw","Rle"], proportions=[0.82,0.18])
demo.add_admixture(time=870, derived="admixwj", ancestral=["Rlllrlrw","Rlemz"], proportions=[0.69,0.31])
demo.sort_events()

r_chrom = 2e-08
r_break = math.log(2)
chrom_positions = [0, 1e+08, 2e+08, 3e+08, 4e+08, 5e+08, 6e+08, 7e+08, 8e+08, 9e+08, 1e+09, 1.1e+09, 1.2e+09, 1.3e+09, 1.4e+09, 1.5e+09, 1.6e+09, 1.7e+09, 1.8e+09, 1.9e+09, 2e+09, 2.1e+09, 2.2e+09, 2.3e+09, 2.4e+09, 2.5e+09, 2.6e+09, 2.7e+09, 2.8e+09, 2.9e+09, 3e+09]
map_positions = [0, 1e+08, 100000001, 2e+08, 200000001, 3e+08, 300000001, 4e+08, 400000001, 5e+08, 500000001, 6e+08, 600000001, 7e+08, 700000001, 8e+08, 800000001, 9e+08, 900000001, 1e+09, 1000000001, 1.1e+09, 1100000001, 1.2e+09, 1200000001, 1.3e+09, 1300000001, 1.4e+09, 1400000001, 1.5e+09, 1500000001, 1.6e+09, 1600000001, 1.7e+09, 1700000001, 1.8e+09, 1800000001, 1.9e+09, 1900000001, 2e+09, 2000000001, 2.1e+09, 2100000001, 2.2e+09, 2200000001, 2.3e+09, 2300000001, 2.4e+09, 2400000001, 2.5e+09, 2500000001, 2.6e+09, 2600000001, 2.7e+09, 2700000001, 2.8e+09, 2800000001, 2.9e+09, 2900000001, 3e+09]
rates = [r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)

ts = msprime.sim_ancestry(samples=samples, model=[msprime.DiscreteTimeWrightFisher(duration=25), msprime.StandardCoalescent()], demography=demo, recombination_rate=rate_map)
tree_sequence = msprime.sim_mutations(ts, rate=1.25e-08, model='binary')
hap_gt = tree_sequence.genotype_matrix()
gt = hap_gt[:, range(0,nhap,2)] + hap_gt[:, range(1,nhap,2)]
nsnps = gt.shape[0]
ts_chroms = numpy.searchsorted(numpy.array(chrom_positions[1:]), tree_sequence.tables.sites.position, side='right') + 1

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_25/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_25/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_25/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

