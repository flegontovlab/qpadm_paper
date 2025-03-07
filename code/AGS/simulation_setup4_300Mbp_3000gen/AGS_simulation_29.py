import numpy
import math
import msprime
import multiprocess

indnam = ["M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("M", "E", "K", "I", "J", "D", "A", "G", "B", "C", "F", "H", "L"), (2438, 2044, 2044, 2228, 1706, 1102, 1331, 1102, 0, 1612, 821, 1102, 821) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 6181, name = "R")
demo.add_population(initial_size = 5540, name = "Rl")
demo.add_population(initial_size = 6043, name = "M")
demo.add_population(initial_size = 9231, name = "Rlr")
demo.add_population(initial_size = 4118, name = "E")
demo.add_population(initial_size = 6182, name = "Rlrr")
demo.add_population(initial_size = 6395, name = "Rr")
demo.add_population(initial_size = 8043, name = "Rrl")
demo.add_population(initial_size = 4676, name = "Rrll")
demo.add_population(initial_size = 1930, name = "K")
demo.add_population(initial_size = 7927, name = "Rrllr")
demo.add_population(initial_size = 1959, name = "I")
demo.add_population(initial_size = 6590, name = "Rrr")
demo.add_population(initial_size = 4882, name = "Rrrl")
demo.add_population(initial_size = 8019, name = "Rrrr")
demo.add_population(initial_size = 3350, name = "J")
demo.add_population(initial_size = 9344, name = "Rrrrr")
demo.add_population(initial_size = 6991, name = "D")
demo.add_population(initial_size = 7580, name = "Rlrrm")
demo.add_population(initial_size = 3364, name = "admix")
demo.add_population(initial_size = 2109, name = "A")
demo.add_population(initial_size = 1941, name = "Rlrrk")
demo.add_population(initial_size = 1131, name = "G")
demo.add_population(initial_size = 6403, name = "admixg")
demo.add_population(initial_size = 7063, name = "B")
demo.add_population(initial_size = 7304, name = "Rlh")
demo.add_population(initial_size = 1038, name = "admixe")
demo.add_population(initial_size = 1629, name = "C")
demo.add_population(initial_size = 5017, name = "Rrrrg")
demo.add_population(initial_size = 6176, name = "admixet")
demo.add_population(initial_size = 5694, name = "Rrrrrd")
demo.add_population(initial_size = 2068, name = "F")
demo.add_population(initial_size = 9000, name = "admixb")
demo.add_population(initial_size = 3947, name = "Rlhr")
demo.add_population(initial_size = 8263, name = "admixr")
demo.add_population(initial_size = 5166, name = "Rrrrgh")
demo.add_population(initial_size = 6719, name = "admixs")
demo.add_population(initial_size = 4701, name = "H")
demo.add_population(initial_size = 3930, name = "Rlhrf")
demo.add_population(initial_size = 1603, name = "admixo")
demo.add_population(initial_size = 3560, name = "Rrb")
demo.add_population(initial_size = 8357, name = "admixn")
demo.add_population(initial_size = 6813, name = "Rlhrfl")
demo.add_population(initial_size = 9242, name = "admixc")
demo.add_population(initial_size = 2655, name = "L")
demo.add_population_split(time=1612, ancestral="admix", derived=["A"])
demo.add_population_split(time=1102, ancestral="admixc", derived=["L"])
demo.add_population_split(time=1706, ancestral="admixe", derived=["C"])
demo.add_population_split(time=338, ancestral="admixg", derived=["B"])
demo.add_population_split(time=1612, ancestral="admixr", derived=["Rlrrk"])
demo.add_population_split(time=1331, ancestral="admixs", derived=["H"])
demo.add_population_split(time=2872, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=2685, ancestral="Rl", derived=["M","Rlh"])
demo.add_population_split(time=2438, ancestral="Rlh", derived=["Rlr","Rlhr"])
demo.add_population_split(time=2228, ancestral="Rlhr", derived=["Rlhrf"])
demo.add_population_split(time=2044, ancestral="Rlhrf", derived=["Rlhrfl"])
demo.add_population_split(time=2228, ancestral="Rlr", derived=["E","Rlrr"])
demo.add_population_split(time=2044, ancestral="Rlrr", derived=["Rlrrm"])
demo.add_population_split(time=1331, ancestral="Rlrrk", derived=["G"])
demo.add_population_split(time=2685, ancestral="Rr", derived=["Rrl","Rrb"])
demo.add_population_split(time=2438, ancestral="Rrb", derived=["Rrr"])
demo.add_population_split(time=2438, ancestral="Rrl", derived=["Rrll","I"])
demo.add_population_split(time=2228, ancestral="Rrll", derived=["K","Rrllr"])
demo.add_population_split(time=2228, ancestral="Rrr", derived=["Rrrl","Rrrr"])
demo.add_population_split(time=2044, ancestral="Rrrr", derived=["J","Rrrrg"])
demo.add_population_split(time=1706, ancestral="Rrrrg", derived=["Rrrrgh"])
demo.add_population_split(time=1612, ancestral="Rrrrgh", derived=["Rrrrr"])
demo.add_population_split(time=1331, ancestral="Rrrrr", derived=["D","Rrrrrd"])
demo.add_population_split(time=1102, ancestral="Rrrrrd", derived=["F"])
demo.add_admixture(time=1612, derived="admix", ancestral=["Rlrrm","Rrrl"], proportions=[0.11,0.89])
demo.add_admixture(time=484, derived="admixb", ancestral=["admixet","Rrrrrd"], proportions=[0.19,0.81])
demo.add_admixture(time=1102, derived="admixc", ancestral=["admixo","Rlhrfl"], proportions=[0.49,0.51])
demo.add_admixture(time=1706, derived="admixe", ancestral=["Rrllr","Rlhr"], proportions=[0.54,0.46])
demo.add_admixture(time=1102, derived="admixet", ancestral=["Rlrrk","Rrrrg"], proportions=[0.49,0.51])
demo.add_admixture(time=338, derived="admixg", ancestral=["admixb","Rlrrm"], proportions=[0.31,0.69])
demo.add_admixture(time=1706, derived="admixn", ancestral=["Rrrl","Rrb"], proportions=[0.77,0.23])
demo.add_admixture(time=1331, derived="admixo", ancestral=["admixn","Rlhrf"], proportions=[0.18,0.82])
demo.add_admixture(time=1612, derived="admixr", ancestral=["Rlhrfl","Rlrr"], proportions=[0.29,0.71])
demo.add_admixture(time=1331, derived="admixs", ancestral=["Rrrrgh","Rrllr"], proportions=[0.7,0.3])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_29.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_29.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_29.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

