import numpy
import math
import msprime
import multiprocess

indnam = ["K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("K", "C", "H", "J", "A", "I", "F", "E", "B", "M", "D", "L", "G"), (600, 600, 1016, 1350, 1350, 1725, 600, 724, 855, 855, 199, 0, 1016) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 1362, name = "R")
demo.add_population(initial_size = 4807, name = "Rl")
demo.add_population(initial_size = 4885, name = "Rllr")
demo.add_population(initial_size = 6022, name = "K")
demo.add_population(initial_size = 2796, name = "C")
demo.add_population(initial_size = 8134, name = "Rlr")
demo.add_population(initial_size = 9868, name = "Rlrl")
demo.add_population(initial_size = 6737, name = "Rlrll")
demo.add_population(initial_size = 3457, name = "H")
demo.add_population(initial_size = 3447, name = "Rlrr")
demo.add_population(initial_size = 2672, name = "Rlrrl")
demo.add_population(initial_size = 2107, name = "J")
demo.add_population(initial_size = 3913, name = "A")
demo.add_population(initial_size = 2482, name = "Rr")
demo.add_population(initial_size = 5755, name = "I")
demo.add_population(initial_size = 9452, name = "Rrr")
demo.add_population(initial_size = 1538, name = "Rrrl")
demo.add_population(initial_size = 9794, name = "F")
demo.add_population(initial_size = 5983, name = "E")
demo.add_population(initial_size = 2095, name = "Rlrlf")
demo.add_population(initial_size = 2669, name = "B")
demo.add_population(initial_size = 8445, name = "Rll")
demo.add_population(initial_size = 5274, name = "admix")
demo.add_population(initial_size = 3637, name = "Rlls")
demo.add_population(initial_size = 2133, name = "M")
demo.add_population(initial_size = 9147, name = "admixl")
demo.add_population(initial_size = 9436, name = "Rlh")
demo.add_population(initial_size = 9741, name = "admixo")
demo.add_population(initial_size = 2346, name = "D")
demo.add_population(initial_size = 7078, name = "Rra")
demo.add_population(initial_size = 1649, name = "admixg")
demo.add_population(initial_size = 9806, name = "Rlrt")
demo.add_population(initial_size = 3314, name = "admixod")
demo.add_population(initial_size = 5764, name = "Rlhc")
demo.add_population(initial_size = 8254, name = "admixla")
demo.add_population(initial_size = 4930, name = "Rlk")
demo.add_population(initial_size = 8485, name = "admixn")
demo.add_population(initial_size = 1216, name = "Ru")
demo.add_population(initial_size = 7187, name = "admixb")
demo.add_population(initial_size = 2561, name = "Rlrrj")
demo.add_population(initial_size = 6191, name = "admixm")
demo.add_population(initial_size = 1877, name = "L")
demo.add_population(initial_size = 9201, name = "Rlrrju")
demo.add_population(initial_size = 4341, name = "G")
demo.add_population(initial_size = 9746, name = "admixw")
demo.add_population_split(time=855, ancestral="admix", derived=["Rllr"])
demo.add_population_split(time=1350, ancestral="admixg", derived=["Rlls"])
demo.add_population_split(time=199, ancestral="admixm", derived=["L"])
demo.add_population_split(time=1350, ancestral="admixn", derived=["Rlrlf"])
demo.add_population_split(time=600, ancestral="admixo", derived=["D"])
demo.add_population_split(time=1016, ancestral="admixod", derived=["Rrr"])
demo.add_population_split(time=2430, ancestral="R", derived=["Rl","Ru"])
demo.add_population_split(time=2062, ancestral="Rl", derived=["Rlr","Rlk"])
demo.add_population_split(time=1725, ancestral="Rlh", derived=["Rll","Rlhc"])
demo.add_population_split(time=1838, ancestral="Rlk", derived=["Rlh"])
demo.add_population_split(time=724, ancestral="Rllr", derived=["K","C"])
demo.add_population_split(time=1016, ancestral="Rlls", derived=["M"])
demo.add_population_split(time=1838, ancestral="Rlr", derived=["Rlrr","Rlrt"])
demo.add_population_split(time=1642, ancestral="Rlrl", derived=["Rlrll"])
demo.add_population_split(time=1016, ancestral="Rlrlf", derived=["B"])
demo.add_population_split(time=1350, ancestral="Rlrll", derived=["H"])
demo.add_population_split(time=1725, ancestral="Rlrr", derived=["Rlrrl","Rlrrj"])
demo.add_population_split(time=1642, ancestral="Rlrrj", derived=["Rlrrju"])
demo.add_population_split(time=1350, ancestral="Rlrrju", derived=["G"])
demo.add_population_split(time=1642, ancestral="Rlrrl", derived=["J","A"])
demo.add_population_split(time=1725, ancestral="Rlrt", derived=["Rlrl"])
demo.add_population_split(time=1838, ancestral="Rr", derived=["I","Rra"])
demo.add_population_split(time=855, ancestral="Rrr", derived=["Rrrl","E"])
demo.add_population_split(time=724, ancestral="Rrrl", derived=["F"])
demo.add_population_split(time=2062, ancestral="Ru", derived=["Rr"])
demo.add_admixture(time=855, derived="admix", ancestral=["Rlrlf","Rll"], proportions=[0.86,0.14])
demo.add_admixture(time=1642, derived="admixb", ancestral=["Rra","Ru"], proportions=[0.39,0.61])
demo.add_admixture(time=1350, derived="admixg", ancestral=["Rll","Rra"], proportions=[0.26,0.74])
demo.add_admixture(time=600, derived="admixl", ancestral=["admixla","Rlls"], proportions=[0.46,0.54])
demo.add_admixture(time=855, derived="admixla", ancestral=["admixw","Rlhc"], proportions=[0.19,0.81])
demo.add_admixture(time=199, derived="admixm", ancestral=["admixl","Rlrrj"], proportions=[0.31,0.69])
demo.add_admixture(time=1350, derived="admixn", ancestral=["Rlrl","Rlk"], proportions=[0.63,0.37])
demo.add_admixture(time=600, derived="admixo", ancestral=["Rrrl","Rlhc"], proportions=[0.5,0.5])
demo.add_admixture(time=1016, derived="admixod", ancestral=["admixb","Rlrt"], proportions=[0.13,0.87])
demo.add_admixture(time=1016, derived="admixw", ancestral=["Rlrll","Rlrrju"], proportions=[0.18,0.82])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_11.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_11.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_11.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

