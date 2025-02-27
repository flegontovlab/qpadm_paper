import numpy
import math
import msprime
import multiprocess

indnam = ["A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("A", "G", "F", "H", "K", "M", "J", "D", "L", "I", "E", "C", "B"), (1864, 1961, 1586, 1740, 1740, 1586, 945, 1864, 945, 0, 442, 686, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 1366, name = "R")
demo.add_population(initial_size = 7692, name = "Rr")
demo.add_population(initial_size = 3019, name = "Rrl")
demo.add_population(initial_size = 2572, name = "Rrll")
demo.add_population(initial_size = 9492, name = "A")
demo.add_population(initial_size = 8895, name = "Rrr")
demo.add_population(initial_size = 8997, name = "Rrrl")
demo.add_population(initial_size = 7212, name = "G")
demo.add_population(initial_size = 2457, name = "Rrrlr")
demo.add_population(initial_size = 2699, name = "Rrrlrl")
demo.add_population(initial_size = 9962, name = "F")
demo.add_population(initial_size = 2239, name = "Rrrr")
demo.add_population(initial_size = 9324, name = "Rrrrl")
demo.add_population(initial_size = 1684, name = "Rrrrr")
demo.add_population(initial_size = 1737, name = "H")
demo.add_population(initial_size = 8218, name = "Rrllv")
demo.add_population(initial_size = 2632, name = "K")
demo.add_population(initial_size = 4982, name = "admix")
demo.add_population(initial_size = 2315, name = "M")
demo.add_population(initial_size = 8043, name = "Rrx")
demo.add_population(initial_size = 3771, name = "admixs")
demo.add_population(initial_size = 5330, name = "Rrrrly")
demo.add_population(initial_size = 6804, name = "Rrrrlr")
demo.add_population(initial_size = 2556, name = "admixb")
demo.add_population(initial_size = 3057, name = "J")
demo.add_population(initial_size = 5124, name = "Rrrn")
demo.add_population(initial_size = 7608, name = "admixe")
demo.add_population(initial_size = 1302, name = "Rrlo")
demo.add_population(initial_size = 9397, name = "D")
demo.add_population(initial_size = 9072, name = "admixsa")
demo.add_population(initial_size = 3615, name = "Rrrrlru")
demo.add_population(initial_size = 3846, name = "L")
demo.add_population(initial_size = 2429, name = "admixy")
demo.add_population(initial_size = 3550, name = "I")
demo.add_population(initial_size = 2270, name = "admixsc")
demo.add_population(initial_size = 2139, name = "admixr")
demo.add_population(initial_size = 1890, name = "E")
demo.add_population(initial_size = 3345, name = "Rrrrlrt")
demo.add_population(initial_size = 7767, name = "C")
demo.add_population(initial_size = 5962, name = "admixm")
demo.add_population(initial_size = 2207, name = "B")
demo.add_population(initial_size = 4007, name = "Rrrrlrti")
demo.add_population(initial_size = 7419, name = "admixp")
demo.add_population(initial_size = 3819, name = "Rrrrlrd")
demo.add_population(initial_size = 2934, name = "admixsp")
demo.add_population_split(time=1740, ancestral="admix", derived=["M"])
demo.add_population_split(time=1350, ancestral="admixb", derived=["J"])
demo.add_population_split(time=1864, ancestral="admixe", derived=["Rrrlrl"])
demo.add_population_split(time=356, ancestral="admixm", derived=["B"])
demo.add_population_split(time=686, ancestral="admixr", derived=["E"])
demo.add_population_split(time=1740, ancestral="admixs", derived=["admixsc"])
demo.add_population_split(time=356, ancestral="admixy", derived=["I"])
demo.add_population_split(time=2966, ancestral="R", derived=["Rr"])
demo.add_population_split(time=2528, ancestral="Rr", derived=["Rrr","Rrx"])
demo.add_population_split(time=2198, ancestral="Rrl", derived=["Rrll","Rrlo"])
demo.add_population_split(time=1961, ancestral="Rrll", derived=["A","Rrllv"])
demo.add_population_split(time=1864, ancestral="Rrllv", derived=["K"])
demo.add_population_split(time=1961, ancestral="Rrlo", derived=["D"])
demo.add_population_split(time=2332, ancestral="Rrr", derived=["Rrrl","Rrrn"])
demo.add_population_split(time=2198, ancestral="Rrrl", derived=["G","Rrrlr"])
demo.add_population_split(time=1740, ancestral="Rrrlrl", derived=["F"])
demo.add_population_split(time=2198, ancestral="Rrrn", derived=["Rrrr"])
demo.add_population_split(time=1961, ancestral="Rrrr", derived=["Rrrrl","Rrrrr"])
demo.add_population_split(time=1864, ancestral="Rrrrl", derived=["Rrrrly"])
demo.add_population_split(time=1586, ancestral="Rrrrlr", derived=["Rrrrlru","Rrrrlrd"])
demo.add_population_split(time=1350, ancestral="Rrrrlrd", derived=["Rrrrlrt"])
demo.add_population_split(time=945, ancestral="Rrrrlrt", derived=["C","Rrrrlrti"])
demo.add_population_split(time=1350, ancestral="Rrrrlru", derived=["L"])
demo.add_population_split(time=1740, ancestral="Rrrrly", derived=["Rrrrlr"])
demo.add_population_split(time=1864, ancestral="Rrrrr", derived=["H"])
demo.add_population_split(time=2332, ancestral="Rrx", derived=["Rrl"])
demo.add_admixture(time=1740, derived="admix", ancestral=["Rrrrl","Rrllv"], proportions=[0.57,0.43])
demo.add_admixture(time=1350, derived="admixb", ancestral=["admixsc","Rrrrly"], proportions=[0.34,0.66])
demo.add_admixture(time=1864, derived="admixe", ancestral=["Rrrlr","Rrrn"], proportions=[0.25,0.75])
demo.add_admixture(time=356, derived="admixm", ancestral=["admixsp","Rrrlr"], proportions=[0.22,0.78])
demo.add_admixture(time=442, derived="admixp", ancestral=["Rrrrlrti","Rrrrlru"], proportions=[0.42,0.58])
demo.add_admixture(time=686, derived="admixr", ancestral=["admixsa","admixsc"], proportions=[0.38,0.62])
demo.add_admixture(time=1740, derived="admixs", ancestral=["Rrrrr","Rrx"], proportions=[0.49,0.51])
demo.add_admixture(time=1586, derived="admixsa", ancestral=["Rrrlrl","Rrlo"], proportions=[0.8,0.2])
demo.add_admixture(time=442, derived="admixsp", ancestral=["Rrrrlrti","Rrrrlrd"], proportions=[0.56,0.44])
demo.add_admixture(time=356, derived="admixy", ancestral=["admixp","R"], proportions=[0.48,0.52])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_38.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_38.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_38.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

