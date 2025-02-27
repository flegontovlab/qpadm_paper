import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "N_1", "N_2", "N_3", "N_4", "N_5", "N_6", "N_7", "N_8", "N_9", "N_10"] ###

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("D", "K", "J", "F", "H", "L", "M", "B", "I", "G", "C", "E", "A", "N"), (268, 268, 268, 443, 61, 0, 0, 443, 187, 61, 187, 312, 337, 443) )] ##
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 7373, name = "R")
demo.add_population(initial_size = 1625, name = "Rl")
demo.add_population(initial_size = 5223, name = "Rll")
demo.add_population(initial_size = 1680, name = "Rllll")
demo.add_population(initial_size = 4897, name = "D")
demo.add_population(initial_size = 9617, name = "Rlllr")
demo.add_population(initial_size = 3447, name = "K")
demo.add_population(initial_size = 1286, name = "J")
demo.add_population(initial_size = 7291, name = "Rlr")
demo.add_population(initial_size = 1982, name = "F")
demo.add_population(initial_size = 8855, name = "Rr")
demo.add_population(initial_size = 3009, name = "Rrl")
demo.add_population(initial_size = 2712, name = "H")
demo.add_population(initial_size = 2041, name = "Rrlr")
demo.add_population(initial_size = 9187, name = "L")
demo.add_population(initial_size = 7283, name = "M")
demo.add_population(initial_size = 9071, name = "Rrr")
demo.add_population(initial_size = 2302, name = "B")
demo.add_population(initial_size = 5300, name = "Rlllv")
demo.add_population(initial_size = 3375, name = "admix")
demo.add_population(initial_size = 8059, name = "I")
demo.add_population(initial_size = 8824, name = "Rllf")
demo.add_population(initial_size = 6370, name = "Rlll")
demo.add_population(initial_size = 7392, name = "admixc")
demo.add_population(initial_size = 4294, name = "G")
demo.add_population(initial_size = 6453, name = "Rllg")
demo.add_population(initial_size = 9836, name = "admixx")
demo.add_population(initial_size = 3951, name = "C")
demo.add_population(initial_size = 7872, name = "Rlli")
demo.add_population(initial_size = 4397, name = "admixxf")
demo.add_population(initial_size = 1663, name = "Rlllvh")
demo.add_population(initial_size = 9719, name = "admixh")
demo.add_population(initial_size = 7318, name = "Rrri")
demo.add_population(initial_size = 6790, name = "Rlrl")
demo.add_population(initial_size = 4836, name = "admixq")
demo.add_population(initial_size = 6529, name = "Rllgl")
demo.add_population(initial_size = 6076, name = "admixv")
demo.add_population(initial_size = 6040, name = "Rlllrh")
demo.add_population(initial_size = 8404, name = "admixci")
demo.add_population(initial_size = 4377, name = "E")
demo.add_population(initial_size = 3176, name = "Rlru")
demo.add_population(initial_size = 3769, name = "admixl")
demo.add_population(initial_size = 5625, name = "Rlrlh")
demo.add_population(initial_size = 9899, name = "A")
demo.add_population(initial_size = 9811, name = "admixu")
demo.add_population(initial_size = 7727, name = "Rr1")    #
demo.add_population(initial_size = 5187, name = "N")      #
demo.add_population_split(time=268, ancestral="admix", derived=["I"])
demo.add_population_split(time=149, ancestral="admixc", derived=["G"])
demo.add_population_split(time=337, ancestral="admixci", derived=["E"])
demo.add_population_split(time=187, ancestral="admixh", derived=["Rrl"])
demo.add_population_split(time=337, ancestral="admixu", derived=["Rlllr"])
demo.add_population_split(time=268, ancestral="admixx", derived=["C"])
demo.add_population_split(time=611, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=527, ancestral="Rl", derived=["Rll","Rlr"])
demo.add_population_split(time=466, ancestral="Rll", derived=["Rllf","Rlli"])
demo.add_population_split(time=443, ancestral="Rllf", derived=["Rlll"])
demo.add_population_split(time=389, ancestral="Rllg", derived=["Rllgl"])
demo.add_population_split(time=443, ancestral="Rlli", derived=["Rllg"])
demo.add_population_split(time=389, ancestral="Rlll", derived=["Rlllrh"])
demo.add_population_split(time=312, ancestral="Rllll", derived=["D"])
demo.add_population_split(time=312, ancestral="Rlllr", derived=["K","J"])
demo.add_population_split(time=368, ancestral="Rlllrh", derived=["Rlllv"])
demo.add_population_split(time=337, ancestral="Rlllv", derived=["Rllll","Rlllvh"])
demo.add_population_split(time=466, ancestral="Rlr", derived=["F","Rlru"])
demo.add_population_split(time=389, ancestral="Rlrl", derived=["Rlrlh"])
demo.add_population_split(time=368, ancestral="Rlrlh", derived=["A"])
demo.add_population_split(time=443, ancestral="Rlru", derived=["Rlrl"])
demo.add_population_split(time=527, ancestral="Rr", derived=["Rrr", "Rr1"]) #
demo.add_population_split(time=466, ancestral="Rr1", derived=["N"])     #
demo.add_population_split(time=149, ancestral="Rrl", derived=["H","Rrlr"])
demo.add_population_split(time=61, ancestral="Rrlr", derived=["L","M"])
demo.add_population_split(time=466, ancestral="Rrr", derived=["B","Rrri"])
demo.add_admixture(time=268, derived="admix", ancestral=["Rlllvh","Rllg"], proportions=[0.55,0.45])
demo.add_admixture(time=149, derived="admixc", ancestral=["Rllll","admixv"], proportions=[0.71,0.29])
demo.add_admixture(time=337, derived="admixci", ancestral=["Rlllrh","Rrri"], proportions=[0.16,0.84])
demo.add_admixture(time=187, derived="admixh", ancestral=["admixl","Rr1"], proportions=[0.22,0.78])   #
demo.add_admixture(time=268, derived="admixl", ancestral=["Rlllvh","Rlru"], proportions=[0.49,0.51])
demo.add_admixture(time=368, derived="admixq", ancestral=["Rlrl","Rrri"], proportions=[0.23,0.77])
demo.add_admixture(time=337, derived="admixu", ancestral=["Rlrlh","Rlll"], proportions=[0.53,0.47])
demo.add_admixture(time=337, derived="admixv", ancestral=["Rllgl","Rllf"], proportions=[0.74,0.26])
demo.add_admixture(time=268, derived="admixx", ancestral=["admixxf","Rllgl"], proportions=[0.25,0.75])
demo.add_admixture(time=337, derived="admixxf", ancestral=["admixq","Rlli"], proportions=[0.49,0.51])
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

numpy.savetxt('./projects/qpadm/results/simout/3chr/sim_9/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/3chr/sim_9/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/3chr/sim_9/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

