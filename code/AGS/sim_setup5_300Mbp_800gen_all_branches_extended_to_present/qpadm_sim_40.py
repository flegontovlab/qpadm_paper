import numpy
import math
import msprime
import multiprocess

indnam = ["I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("I", "D", "L", "M", "B", "F", "E", "C", "G", "A", "K", "J", "H"), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 6620, name = "Rl")
demo.add_population(initial_size = 1770, name = "Rll")
demo.add_population(initial_size = 4126, name = "Rlll")
demo.add_population(initial_size = 3720, name = "I")
demo.add_population(initial_size = 9557, name = "Rlllr")
demo.add_population(initial_size = 5608, name = "Rlr")
demo.add_population(initial_size = 7708, name = "D")
demo.add_population(initial_size = 9381, name = "Rlrr")
demo.add_population(initial_size = 2670, name = "Rlrrl")
demo.add_population(initial_size = 5636, name = "L")
demo.add_population(initial_size = 3911, name = "M")
demo.add_population(initial_size = 9669, name = "R")
demo.add_population(initial_size = 4722, name = "Rr")
demo.add_population(initial_size = 6284, name = "Rrl")
demo.add_population(initial_size = 3128, name = "B")
demo.add_population(initial_size = 5882, name = "Rrr")
demo.add_population(initial_size = 3937, name = "Rrrl")
demo.add_population(initial_size = 4068, name = "F")
demo.add_population(initial_size = 3478, name = "admix")
demo.add_population(initial_size = 5579, name = "Rlllrj")
demo.add_population(initial_size = 5909, name = "admixt")
demo.add_population(initial_size = 2071, name = "E")
demo.add_population(initial_size = 8972, name = "Rlllrt")
demo.add_population(initial_size = 4763, name = "Rrlh")
demo.add_population(initial_size = 7245, name = "C")
demo.add_population(initial_size = 1699, name = "admixl")
demo.add_population(initial_size = 3969, name = "Rc")
demo.add_population(initial_size = 6180, name = "admixo")
demo.add_population(initial_size = 2925, name = "G")
demo.add_population(initial_size = 2075, name = "Rrrlh")
demo.add_population(initial_size = 8500, name = "A")
demo.add_population(initial_size = 2841, name = "admixp")
demo.add_population(initial_size = 6541, name = "Rrz")
demo.add_population(initial_size = 2580, name = "admixr")
demo.add_population(initial_size = 1918, name = "Rro")
demo.add_population(initial_size = 9233, name = "admixv")
demo.add_population(initial_size = 2668, name = "K")
demo.add_population(initial_size = 8486, name = "Rrrlm")
demo.add_population(initial_size = 8609, name = "admixh")
demo.add_population(initial_size = 4785, name = "Rlrrd")
demo.add_population(initial_size = 6406, name = "admixk")
demo.add_population(initial_size = 1283, name = "J")
demo.add_population(initial_size = 8291, name = "Rlllrjo")
demo.add_population(initial_size = 1195, name = "admixht")
demo.add_population(initial_size = 9777, name = "H")
demo.add_population_split(time=203, ancestral="admixht", derived=["H"])
demo.add_population_split(time=54, ancestral="admixk", derived=["J"])
demo.add_population_split(time=379, ancestral="admixl", derived=["Rlllrt"])
demo.add_population_split(time=203, ancestral="admixo", derived=["G"])
demo.add_population_split(time=203, ancestral="admixt", derived=["E"])
demo.add_population_split(time=203, ancestral="admixv", derived=["K"])
demo.add_population_split(time=751, ancestral="R", derived=["Rr","Rc"])
demo.add_population_split(time=656, ancestral="Rc", derived=["Rl"])
demo.add_population_split(time=554, ancestral="Rl", derived=["Rll","Rlr"])
demo.add_population_split(time=532, ancestral="Rll", derived=["Rlll"])
demo.add_population_split(time=506, ancestral="Rlll", derived=["I","Rlllr"])
demo.add_population_split(time=419, ancestral="Rlllr", derived=["Rlllrjo"])
demo.add_population_split(time=379, ancestral="Rlllrjo", derived=["Rlllrj"])
demo.add_population_split(time=532, ancestral="Rlr", derived=["D","Rlrr"])
demo.add_population_split(time=506, ancestral="Rlrr", derived=["M","Rlrrd"])
demo.add_population_split(time=419, ancestral="Rlrrd", derived=["Rlrrl"])
demo.add_population_split(time=379, ancestral="Rlrrl", derived=["L"])
demo.add_population_split(time=656, ancestral="Rr", derived=["Rrr","Rro"])
demo.add_population_split(time=506, ancestral="Rrl", derived=["B","Rrlh"])
demo.add_population_split(time=419, ancestral="Rrlh", derived=["C"])
demo.add_population_split(time=554, ancestral="Rro", derived=["Rrz"])
demo.add_population_split(time=554, ancestral="Rrr", derived=["Rrrl","F"])
demo.add_population_split(time=532, ancestral="Rrrl", derived=["Rrrlm"])
demo.add_population_split(time=419, ancestral="Rrrlh", derived=["A"])
demo.add_population_split(time=506, ancestral="Rrrlm", derived=["Rrrlh"])
demo.add_population_split(time=532, ancestral="Rrz", derived=["Rrl"])
demo.add_admixture(time=132, derived="admix", ancestral=["admixh","Rll"], proportions=[0.47,0.53])
demo.add_admixture(time=203, derived="admixh", ancestral=["Rlllrj","Rrrlm"], proportions=[0.2,0.8])
demo.add_admixture(time=203, derived="admixht", ancestral=["Rlllrt","Rlllrjo"], proportions=[0.62,0.38])
demo.add_admixture(time=54, derived="admixk", ancestral=["admix","Rlrrd"], proportions=[0.34,0.66])
demo.add_admixture(time=379, derived="admixl", ancestral=["Rlllr","Rrlh"], proportions=[0.47,0.53])
demo.add_admixture(time=203, derived="admixo", ancestral=["admixp","Rrrl"], proportions=[0.46,0.54])
demo.add_admixture(time=379, derived="admixp", ancestral=["Rrrlh","Rc"], proportions=[0.64,0.36])
demo.add_admixture(time=263, derived="admixr", ancestral=["Rlrrl","Rrz"], proportions=[0.23,0.77])
demo.add_admixture(time=203, derived="admixt", ancestral=["Rlllrj","Rlllrt"], proportions=[0.27,0.73])
demo.add_admixture(time=203, derived="admixv", ancestral=["admixr","Rro"], proportions=[0.41,0.59])
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

numpy.savetxt('./projects/qpadm/results/simout/qpadm_sim_fixed_40.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_40.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_40.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

