import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("A", "G", "F", "I", "L", "M", "E", "B", "H", "J", "K", "C", "D"), (485, 112, 330, 381, 381, 187, 213, 0, 422, 259, 112, 187, 381) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 9391, name = "Rl")
demo.add_population(initial_size = 5886, name = "Rll")
demo.add_population(initial_size = 1861, name = "Rlr")
demo.add_population(initial_size = 8509, name = "A")
demo.add_population(initial_size = 4223, name = "Rlrr")
demo.add_population(initial_size = 5986, name = "Rlrrl")
demo.add_population(initial_size = 9990, name = "G")
demo.add_population(initial_size = 8732, name = "Rlrrr")
demo.add_population(initial_size = 6071, name = "F")
demo.add_population(initial_size = 8052, name = "Rrl")
demo.add_population(initial_size = 9174, name = "I")
demo.add_population(initial_size = 8764, name = "L")
demo.add_population(initial_size = 9445, name = "Rrrl")
demo.add_population(initial_size = 5457, name = "M")
demo.add_population(initial_size = 1734, name = "Rrr")
demo.add_population(initial_size = 4125, name = "Rrrr")
demo.add_population(initial_size = 5632, name = "E")
demo.add_population(initial_size = 8598, name = "Rrj")
demo.add_population(initial_size = 6119, name = "admix")
demo.add_population(initial_size = 4159, name = "B")
demo.add_population(initial_size = 2190, name = "Rlrrrf")
demo.add_population(initial_size = 5691, name = "admixf")
demo.add_population(initial_size = 5427, name = "Rrrd")
demo.add_population(initial_size = 1215, name = "admixi")
demo.add_population(initial_size = 6034, name = "R")
demo.add_population(initial_size = 8109, name = "Rb")
demo.add_population(initial_size = 7676, name = "admixz")
demo.add_population(initial_size = 7222, name = "H")
demo.add_population(initial_size = 7654, name = "Rz")
demo.add_population(initial_size = 3674, name = "Rr")
demo.add_population(initial_size = 4260, name = "admixo")
demo.add_population(initial_size = 4714, name = "J")
demo.add_population(initial_size = 3931, name = "Rlrrt")
demo.add_population(initial_size = 2308, name = "admixw")
demo.add_population(initial_size = 6285, name = "K")
demo.add_population(initial_size = 6677, name = "Rbj")
demo.add_population(initial_size = 3595, name = "admixr")
demo.add_population(initial_size = 9260, name = "Rrm")
demo.add_population(initial_size = 4057, name = "admixg")
demo.add_population(initial_size = 4079, name = "C")
demo.add_population(initial_size = 2747, name = "Rrg")
demo.add_population(initial_size = 4629, name = "admixu")
demo.add_population(initial_size = 1802, name = "admixrk")
demo.add_population(initial_size = 5393, name = "D")
demo.add_population(initial_size = 2256, name = "admixm")
demo.add_population_split(time=112, ancestral="admix", derived=["B"])
demo.add_population_split(time=381, ancestral="admixf", derived=["Rrr"])
demo.add_population_split(time=213, ancestral="admixg", derived=["C"])
demo.add_population_split(time=187, ancestral="admixi", derived=["Rlrrl"])
demo.add_population_split(time=381, ancestral="admixm", derived=["Rlrrt"])
demo.add_population_split(time=330, ancestral="admixo", derived=["J"])
demo.add_population_split(time=485, ancestral="admixr", derived=["admixrk"])
demo.add_population_split(time=422, ancestral="admixrk", derived=["D"])
demo.add_population_split(time=164, ancestral="admixw", derived=["K"])
demo.add_population_split(time=485, ancestral="admixz", derived=["H"])
demo.add_population_split(time=752, ancestral="R", derived=["Rb","Rz"])
demo.add_population_split(time=703, ancestral="Rb", derived=["Rl","Rbj"])
demo.add_population_split(time=656, ancestral="Rl", derived=["Rll","Rlr"])
demo.add_population_split(time=571, ancestral="Rlr", derived=["A","Rlrr"])
demo.add_population_split(time=485, ancestral="Rlrr", derived=["Rlrrrf"])
demo.add_population_split(time=164, ancestral="Rlrrl", derived=["G"])
demo.add_population_split(time=381, ancestral="Rlrrr", derived=["F"])
demo.add_population_split(time=422, ancestral="Rlrrrf", derived=["Rlrrr"])
demo.add_population_split(time=656, ancestral="Rr", derived=["Rrm","Rrg"])
demo.add_population_split(time=571, ancestral="Rrg", derived=["Rrj"])
demo.add_population_split(time=485, ancestral="Rrj", derived=["Rrl"])
demo.add_population_split(time=422, ancestral="Rrl", derived=["I","L"])
demo.add_population_split(time=330, ancestral="Rrr", derived=["Rrrr","Rrrd"])
demo.add_population_split(time=259, ancestral="Rrrd", derived=["Rrrl"])
demo.add_population_split(time=213, ancestral="Rrrl", derived=["M"])
demo.add_population_split(time=259, ancestral="Rrrr", derived=["E"])
demo.add_population_split(time=703, ancestral="Rz", derived=["Rr"])
demo.add_admixture(time=112, derived="admix", ancestral=["Rlrrl","Rrj"], proportions=[0.18,0.82])
demo.add_admixture(time=381, derived="admixf", ancestral=["Rlrrrf","Rrm"], proportions=[0.41,0.59])
demo.add_admixture(time=213, derived="admixg", ancestral=["Rrrr","admixu"], proportions=[0.51,0.49])
demo.add_admixture(time=187, derived="admixi", ancestral=["Rrrd","Rlrrt"], proportions=[0.9,0.1])
demo.add_admixture(time=381, derived="admixm", ancestral=["admixrk","Rlrr"], proportions=[0.46,0.54])
demo.add_admixture(time=330, derived="admixo", ancestral=["Rlrrr","Rz"], proportions=[0.5,0.5])
demo.add_admixture(time=485, derived="admixr", ancestral=["Rll","Rbj"], proportions=[0.64,0.36])
demo.add_admixture(time=485, derived="admixu", ancestral=["Rrm","Rrg"], proportions=[0.58,0.42])
demo.add_admixture(time=164, derived="admixw", ancestral=["Rrrl","Rlrrt"], proportions=[0.73,0.27])
demo.add_admixture(time=485, derived="admixz", ancestral=["Rll","Rbj"], proportions=[0.57,0.43])
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

numpy.savetxt('./projects/qpadm/results/simout/3chr/sim_14/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/3chr/sim_14/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/3chr/sim_14/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

