import numpy
import math
import msprime
import multiprocess
import sys

if len(sys.argv)  != 2:
  raise ValueError("Need to specify simulation no")
sim_id = str(sys.argv[1])

indnam = ["C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("C", "I", "A", "H", "K", "B", "D", "G", "F", "L", "J", "E", "M"), (2141, 1819, 506, 506, 412, 1538, 1538, 0, 412, 874, 1538, 506, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 2843, name = "R")
demo.add_population(initial_size = 7578, name = "Rr")
demo.add_population(initial_size = 9772, name = "Rrl")
demo.add_population(initial_size = 1047, name = "Rrll")
demo.add_population(initial_size = 4026, name = "C")
demo.add_population(initial_size = 4694, name = "Rrlr")
demo.add_population(initial_size = 4385, name = "Rrlrr")
demo.add_population(initial_size = 5124, name = "I")
demo.add_population(initial_size = 6677, name = "Rrrl")
demo.add_population(initial_size = 2822, name = "A")
demo.add_population(initial_size = 9059, name = "Rrrr")
demo.add_population(initial_size = 5682, name = "H")
demo.add_population(initial_size = 5108, name = "Rrrrr")
demo.add_population(initial_size = 1352, name = "K")
demo.add_population(initial_size = 3777, name = "Rrllf")
demo.add_population(initial_size = 3637, name = "admix")
demo.add_population(initial_size = 2279, name = "B")
demo.add_population(initial_size = 9064, name = "Rl")
demo.add_population(initial_size = 8887, name = "Rlx")
demo.add_population(initial_size = 6679, name = "D")
demo.add_population(initial_size = 9714, name = "admixx")
demo.add_population(initial_size = 2462, name = "Rlr")
demo.add_population(initial_size = 4307, name = "Rll")
demo.add_population(initial_size = 5375, name = "admixs")
demo.add_population(initial_size = 3001, name = "G")
demo.add_population(initial_size = 9024, name = "Rrg")
demo.add_population(initial_size = 7180, name = "Rrr")
demo.add_population(initial_size = 1986, name = "admixb")
demo.add_population(initial_size = 1261, name = "Rj")
demo.add_population(initial_size = 7996, name = "admixsv")
demo.add_population(initial_size = 3626, name = "F")
demo.add_population(initial_size = 7105, name = "Rllr")
demo.add_population(initial_size = 1347, name = "L")
demo.add_population(initial_size = 1969, name = "admixbm")
demo.add_population(initial_size = 2329, name = "Rrre")
demo.add_population(initial_size = 2569, name = "admixu")
demo.add_population(initial_size = 2877, name = "J")
demo.add_population(initial_size = 6349, name = "admixxi")
demo.add_population(initial_size = 9935, name = "E")
demo.add_population(initial_size = 9131, name = "admixj")
demo.add_population(initial_size = 8059, name = "Rln")
demo.add_population(initial_size = 2103, name = "admixr")
demo.add_population(initial_size = 2614, name = "M")
demo.add_population(initial_size = 7033, name = "Rlla")
demo.add_population(initial_size = 2691, name = "admixp")
demo.add_population_split(time=2865, ancestral="R", derived=["Rr","Rj"])
demo.add_population_split(time=2141, ancestral="Rl", derived=["Rlx","Rln"])
demo.add_population_split(time=1342, ancestral="Rll", derived=["Rllr","Rlla"])
demo.add_population_split(time=956, ancestral="Rllr", derived=["L"])
demo.add_population_split(time=1819, ancestral="Rln", derived=["Rlr"])
demo.add_population_split(time=1538, ancestral="Rlr", derived=["Rll"])
demo.add_population_split(time=1819, ancestral="Rlx", derived=["D"])
demo.add_population_split(time=2674, ancestral="Rr", derived=["Rrl","Rrg"])
demo.add_population_split(time=2542, ancestral="Rrg", derived=["Rrr"])
demo.add_population_split(time=2542, ancestral="Rrl", derived=["Rrll","Rrlr"])
demo.add_population_split(time=2378, ancestral="Rrll", derived=["C","Rrllf"])
demo.add_population_split(time=2378, ancestral="Rrlr", derived=["Rrlrr"])
demo.add_population_split(time=2141, ancestral="Rrlrr", derived=["I"])
demo.add_population_split(time=2378, ancestral="Rrr", derived=["Rrre"])
demo.add_population_split(time=638, ancestral="Rrrl", derived=["A"])
demo.add_population_split(time=638, ancestral="Rrrr", derived=["H","Rrrrr"])
demo.add_population_split(time=506, ancestral="Rrrrr", derived=["K"])
demo.add_population_split(time=1819, ancestral="admix", derived=["B"])
demo.add_population_split(time=2378, ancestral="admixb", derived=["Rl"])
demo.add_population_split(time=874, ancestral="admixbm", derived=["Rrrl"])
demo.add_population_split(time=874, ancestral="admixp", derived=["Rrrr"])
demo.add_population_split(time=412, ancestral="admixr", derived=["M"])
demo.add_population_split(time=412, ancestral="admixs", derived=["G"])
demo.add_population_split(time=506, ancestral="admixsv", derived=["F"])
demo.add_population_split(time=1819, ancestral="admixu", derived=["J"])
demo.add_population_split(time=874, ancestral="admixx", derived=["admixxi"])
demo.add_population_split(time=638, ancestral="admixxi", derived=["E"])
demo.add_admixture(time=1819, derived="admix", ancestral=["Rrlrr","Rrllf"], proportions=[0.72,0.28])
demo.add_admixture(time=2378, derived="admixb", ancestral=["Rrg","Rj"], proportions=[0.53,0.47])
demo.add_admixture(time=874, derived="admixbm", ancestral=["Rllr","Rrre"], proportions=[0.47,0.53])
demo.add_admixture(time=506, derived="admixj", ancestral=["admixxi","Rlr"], proportions=[0.34,0.66])
demo.add_admixture(time=874, derived="admixp", ancestral=["Rlla","Rrr"], proportions=[0.19,0.81])
demo.add_admixture(time=412, derived="admixr", ancestral=["Rrrrr","Rln"], proportions=[0.88,0.12])
demo.add_admixture(time=412, derived="admixs", ancestral=["admixj","Rrllf"], proportions=[0.12,0.88])
demo.add_admixture(time=506, derived="admixsv", ancestral=["Rrrl","Rj"], proportions=[0.54,0.46])
demo.add_admixture(time=1819, derived="admixu", ancestral=["Rrre","Rrlr"], proportions=[0.57,0.43])
demo.add_admixture(time=874, derived="admixx", ancestral=["Rlla","Rlx"], proportions=[0.21,0.79])
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

numpy.savetxt('./projects/qpadm/results/simout/30chr_deep/sim_17/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/30chr_deep/sim_17/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/30chr_deep/sim_17/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

