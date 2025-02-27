import numpy
import math
import msprime
import multiprocess

indnam = ["G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("G", "A", "I", "K", "F", "E", "L", "J", "H", "D", "M", "C", "B"), (1898, 1898, 622, 622, 622, 1399, 1399, 622, 334, 0, 892, 622, 334) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 1361, name = "R")
demo.add_population(initial_size = 1626, name = "Rr")
demo.add_population(initial_size = 7194, name = "Rrr")
demo.add_population(initial_size = 3292, name = "Rrrl")
demo.add_population(initial_size = 5152, name = "G")
demo.add_population(initial_size = 9705, name = "A")
demo.add_population(initial_size = 7792, name = "Rrrr")
demo.add_population(initial_size = 4896, name = "Rrrrl")
demo.add_population(initial_size = 3064, name = "Rrrrll")
demo.add_population(initial_size = 5099, name = "Rrrrllll")
demo.add_population(initial_size = 7029, name = "I")
demo.add_population(initial_size = 6554, name = "Rrrrlllr")
demo.add_population(initial_size = 6019, name = "K")
demo.add_population(initial_size = 2792, name = "F")
demo.add_population(initial_size = 2392, name = "Rrrrr")
demo.add_population(initial_size = 6961, name = "Rrrrrb")
demo.add_population(initial_size = 5031, name = "E")
demo.add_population(initial_size = 2691, name = "Rrrrlll")
demo.add_population(initial_size = 6874, name = "admix")
demo.add_population(initial_size = 2709, name = "Rrrrle")
demo.add_population(initial_size = 3522, name = "L")
demo.add_population(initial_size = 1002, name = "admixs")
demo.add_population(initial_size = 7233, name = "Rrrrrw")
demo.add_population(initial_size = 6418, name = "admixh")
demo.add_population(initial_size = 8053, name = "Rrn")
demo.add_population(initial_size = 8373, name = "admixf")
demo.add_population(initial_size = 7877, name = "J")
demo.add_population(initial_size = 8666, name = "Rrrrlld")
demo.add_population(initial_size = 7674, name = "Rrrrllr")
demo.add_population(initial_size = 4909, name = "admixe")
demo.add_population(initial_size = 2202, name = "Rrrrllllp")
demo.add_population(initial_size = 4241, name = "H")
demo.add_population(initial_size = 6131, name = "admixeo")
demo.add_population(initial_size = 4595, name = "Rrrrllrb")
demo.add_population(initial_size = 6175, name = "admixp")
demo.add_population(initial_size = 1312, name = "D")
demo.add_population(initial_size = 8987, name = "Rrrrrbb")
demo.add_population(initial_size = 1897, name = "admixb")
demo.add_population(initial_size = 8260, name = "Rrrrrbbx")
demo.add_population(initial_size = 5734, name = "admixw")
demo.add_population(initial_size = 7755, name = "M")
demo.add_population(initial_size = 8550, name = "Rrrrllrba")
demo.add_population(initial_size = 4958, name = "C")
demo.add_population(initial_size = 7494, name = "admixd")
demo.add_population(initial_size = 5699, name = "B")
demo.add_population_split(time=975, ancestral="admix", derived=["Rrrrlllr"])
demo.add_population_split(time=622, ancestral="admixd", derived=["B"])
demo.add_population_split(time=892, ancestral="admixf", derived=["J"])
demo.add_population_split(time=1399, ancestral="admixh", derived=["Rrrrlll"])
demo.add_population_split(time=176, ancestral="admixp", derived=["D"])
demo.add_population_split(time=975, ancestral="admixs", derived=["Rrrrllll"])
demo.add_population_split(time=975, ancestral="admixw", derived=["M"])
demo.add_population_split(time=2996, ancestral="R", derived=["Rr"])
demo.add_population_split(time=2632, ancestral="Rr", derived=["Rrr","Rrn"])
demo.add_population_split(time=2288, ancestral="Rrr", derived=["Rrrl","Rrrr"])
demo.add_population_split(time=2198, ancestral="Rrrl", derived=["G","A"])
demo.add_population_split(time=2198, ancestral="Rrrr", derived=["Rrrrl","Rrrrr"])
demo.add_population_split(time=1898, ancestral="Rrrrl", derived=["Rrrrll","Rrrrle"])
demo.add_population_split(time=1714, ancestral="Rrrrle", derived=["L"])
demo.add_population_split(time=1714, ancestral="Rrrrll", derived=["Rrrrlld"])
demo.add_population_split(time=1399, ancestral="Rrrrlld", derived=["Rrrrllr"])
demo.add_population_split(time=892, ancestral="Rrrrllll", derived=["I","Rrrrllllp"])
demo.add_population_split(time=622, ancestral="Rrrrllllp", derived=["H"])
demo.add_population_split(time=892, ancestral="Rrrrlllr", derived=["K","F"])
demo.add_population_split(time=1166, ancestral="Rrrrllr", derived=["Rrrrllrb"])
demo.add_population_split(time=975, ancestral="Rrrrllrb", derived=["Rrrrllrba"])
demo.add_population_split(time=892, ancestral="Rrrrllrba", derived=["C"])
demo.add_population_split(time=1898, ancestral="Rrrrr", derived=["Rrrrrb","Rrrrrw"])
demo.add_population_split(time=1714, ancestral="Rrrrrb", derived=["E","Rrrrrbb"])
demo.add_population_split(time=1399, ancestral="Rrrrrbb", derived=["Rrrrrbbx"])
demo.add_admixture(time=975, derived="admix", ancestral=["Rrrrlll","Rrrrrbbx"], proportions=[0.26,0.74])
demo.add_admixture(time=1166, derived="admixb", ancestral=["Rrrrrbb","Rrn"], proportions=[0.79,0.21])
demo.add_admixture(time=622, derived="admixd", ancestral=["Rrrrllrba","Rrrrrw"], proportions=[0.16,0.84])
demo.add_admixture(time=1166, derived="admixe", ancestral=["Rrrrlld","R"], proportions=[0.57,0.43])
demo.add_admixture(time=334, derived="admixeo", ancestral=["Rrrrllllp","Rrrrllr"], proportions=[0.34,0.66])
demo.add_admixture(time=892, derived="admixf", ancestral=["admixe","Rrn"], proportions=[0.23,0.77])
demo.add_admixture(time=1399, derived="admixh", ancestral=["Rrrrll","Rrrrrw"], proportions=[0.38,0.62])
demo.add_admixture(time=176, derived="admixp", ancestral=["admixeo","Rrrrllrb"], proportions=[0.32,0.68])
demo.add_admixture(time=975, derived="admixs", ancestral=["Rrrrlll","Rrrrle"], proportions=[0.53,0.47])
demo.add_admixture(time=975, derived="admixw", ancestral=["admixb","Rrrrrbbx"], proportions=[0.16,0.84])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_33.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_33.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_33.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

