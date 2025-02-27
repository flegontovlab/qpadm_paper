import numpy
import math
import msprime
import multiprocess

indnam = ["B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("B", "M", "F", "D", "J", "L", "E", "C", "G", "I", "H", "A", "K"), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 9197, name = "R")
demo.add_population(initial_size = 8861, name = "Rl")
demo.add_population(initial_size = 9189, name = "Rlr")
demo.add_population(initial_size = 4330, name = "B")
demo.add_population(initial_size = 9101, name = "Rr")
demo.add_population(initial_size = 8247, name = "Rrl")
demo.add_population(initial_size = 3827, name = "Rrlr")
demo.add_population(initial_size = 1765, name = "Rrlrl")
demo.add_population(initial_size = 4109, name = "Rrlrll")
demo.add_population(initial_size = 3281, name = "Rrlrlll")
demo.add_population(initial_size = 6949, name = "M")
demo.add_population(initial_size = 5860, name = "F")
demo.add_population(initial_size = 1393, name = "Rrr")
demo.add_population(initial_size = 3536, name = "Rrrl")
demo.add_population(initial_size = 1313, name = "Rrrll")
demo.add_population(initial_size = 4809, name = "D")
demo.add_population(initial_size = 4821, name = "J")
demo.add_population(initial_size = 7933, name = "Rrlrlh")
demo.add_population(initial_size = 7577, name = "L")
demo.add_population(initial_size = 8841, name = "admix")
demo.add_population(initial_size = 5189, name = "E")
demo.add_population(initial_size = 1151, name = "Rrro")
demo.add_population(initial_size = 9047, name = "admixo")
demo.add_population(initial_size = 4899, name = "C")
demo.add_population(initial_size = 8623, name = "Rlrc")
demo.add_population(initial_size = 4893, name = "G")
demo.add_population(initial_size = 8971, name = "admixc")
demo.add_population(initial_size = 5912, name = "Rrrllm")
demo.add_population(initial_size = 8443, name = "I")
demo.add_population(initial_size = 3290, name = "admixi")
demo.add_population(initial_size = 6254, name = "Rrlg")
demo.add_population(initial_size = 8609, name = "admixv")
demo.add_population(initial_size = 4478, name = "Rrlrb")
demo.add_population(initial_size = 5523, name = "admixj")
demo.add_population(initial_size = 6058, name = "H")
demo.add_population(initial_size = 1889, name = "Rd")
demo.add_population(initial_size = 6312, name = "admixk")
demo.add_population(initial_size = 3285, name = "A")
demo.add_population(initial_size = 3325, name = "Rrlrbx")
demo.add_population(initial_size = 9751, name = "admixm")
demo.add_population(initial_size = 9360, name = "Rrroj")
demo.add_population(initial_size = 8001, name = "admixw")
demo.add_population(initial_size = 7218, name = "K")
demo.add_population(initial_size = 2173, name = "Rrroh")
demo.add_population(initial_size = 1788, name = "admixic")
demo.add_population_split(time=366, ancestral="admix", derived=["E"])
demo.add_population_split(time=390, ancestral="admixj", derived=["H"])
demo.add_population_split(time=40, ancestral="admixk", derived=["A"])
demo.add_population_split(time=210, ancestral="admixo", derived=["C"])
demo.add_population_split(time=390, ancestral="admixw", derived=["K"])
demo.add_population_split(time=800, ancestral="R", derived=["Rl","Rd"])
demo.add_population_split(time=760, ancestral="Rd", derived=["Rr"])
demo.add_population_split(time=760, ancestral="Rl", derived=["Rlr","Rlrc"])
demo.add_population_split(time=714, ancestral="Rlr", derived=["B"])
demo.add_population_split(time=714, ancestral="Rlrc", derived=["G"])
demo.add_population_split(time=714, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=671, ancestral="Rrl", derived=["Rrlr","Rrlg"])
demo.add_population_split(time=558, ancestral="Rrlr", derived=["Rrlrl","Rrlrb"])
demo.add_population_split(time=462, ancestral="Rrlrb", derived=["Rrlrbx"])
demo.add_population_split(time=462, ancestral="Rrlrl", derived=["Rrlrll","Rrlrlh"])
demo.add_population_split(time=390, ancestral="Rrlrlh", derived=["L"])
demo.add_population_split(time=390, ancestral="Rrlrll", derived=["Rrlrlll","F"])
demo.add_population_split(time=366, ancestral="Rrlrlll", derived=["M"])
demo.add_population_split(time=671, ancestral="Rrr", derived=["Rrrl","Rrro"])
demo.add_population_split(time=558, ancestral="Rrrl", derived=["Rrrll","J"])
demo.add_population_split(time=462, ancestral="Rrrll", derived=["D","Rrrllm"])
demo.add_population_split(time=390, ancestral="Rrrllm", derived=["I"])
demo.add_population_split(time=558, ancestral="Rrro", derived=["Rrroj","Rrroh"])
demo.add_admixture(time=366, derived="admix", ancestral=["Rrlrlh","Rrlrbx"], proportions=[0.61,0.39])
demo.add_admixture(time=462, derived="admixc", ancestral=["Rrlg","Rlrc"], proportions=[0.64,0.36])
demo.add_admixture(time=305, derived="admixi", ancestral=["Rrlrlll","Rrrllm"], proportions=[0.32,0.68])
demo.add_admixture(time=366, derived="admixic", ancestral=["Rrlrbx","Rrroh"], proportions=[0.3,0.7])
demo.add_admixture(time=390, derived="admixj", ancestral=["Rrlrb","Rrroh"], proportions=[0.3,0.7])
demo.add_admixture(time=40, derived="admixk", ancestral=["admixm","admixc"], proportions=[0.89,0.11])
demo.add_admixture(time=305, derived="admixm", ancestral=["admixic","Rd"], proportions=[0.22,0.78])
demo.add_admixture(time=210, derived="admixo", ancestral=["admixv","Rrroj"], proportions=[0.21,0.79])
demo.add_admixture(time=255, derived="admixv", ancestral=["admixi","Rrlg"], proportions=[0.21,0.79])
demo.add_admixture(time=390, derived="admixw", ancestral=["Rrroj","Rlr"], proportions=[0.58,0.42])
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

numpy.savetxt('./projects/qpadm/results/simout/qpadm_sim_fixed_12.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_12.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_12.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

