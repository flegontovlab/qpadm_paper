import numpy
import math
import msprime
import multiprocess

indnam = ["I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("I", "B", "M", "K", "A", "F", "H", "J", "L", "E", "C", "D", "G"), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 2913, name = "Rl")
demo.add_population(initial_size = 1868, name = "I")
demo.add_population(initial_size = 3379, name = "R")
demo.add_population(initial_size = 3419, name = "Rr")
demo.add_population(initial_size = 6401, name = "Rrl")
demo.add_population(initial_size = 3659, name = "Rrll")
demo.add_population(initial_size = 4680, name = "Rrlll")
demo.add_population(initial_size = 9752, name = "B")
demo.add_population(initial_size = 5999, name = "Rrllr")
demo.add_population(initial_size = 9165, name = "Rrllrl")
demo.add_population(initial_size = 7045, name = "M")
demo.add_population(initial_size = 4353, name = "K")
demo.add_population(initial_size = 2821, name = "Rrllrr")
demo.add_population(initial_size = 2615, name = "A")
demo.add_population(initial_size = 7330, name = "F")
demo.add_population(initial_size = 9233, name = "Rrlr")
demo.add_population(initial_size = 8041, name = "Rrlrl")
demo.add_population(initial_size = 9663, name = "Rrlrr")
demo.add_population(initial_size = 9620, name = "H")
demo.add_population(initial_size = 6839, name = "Rrlj")
demo.add_population(initial_size = 7160, name = "admix")
demo.add_population(initial_size = 7839, name = "J")
demo.add_population(initial_size = 5281, name = "admixz")
demo.add_population(initial_size = 5123, name = "Rrm")
demo.add_population(initial_size = 5085, name = "Rrlp")
demo.add_population(initial_size = 3509, name = "admixj")
demo.add_population(initial_size = 5908, name = "L")
demo.add_population(initial_size = 6023, name = "Rm")
demo.add_population(initial_size = 2944, name = "admixzy")
demo.add_population(initial_size = 8680, name = "Rlz")
demo.add_population(initial_size = 9234, name = "E")
demo.add_population(initial_size = 3765, name = "admixv")
demo.add_population(initial_size = 1132, name = "C")
demo.add_population(initial_size = 9897, name = "Rlb")
demo.add_population(initial_size = 3684, name = "admixzt")
demo.add_population(initial_size = 8039, name = "D")
demo.add_population(initial_size = 2023, name = "Rrlrz")
demo.add_population(initial_size = 9973, name = "admixd")
demo.add_population(initial_size = 1487, name = "G")
demo.add_population(initial_size = 6146, name = "Rlzn")
demo.add_population(initial_size = 5796, name = "admixx")
demo.add_population(initial_size = 7572, name = "Rrlw")
demo.add_population(initial_size = 3907, name = "admixjy")
demo.add_population(initial_size = 1733, name = "Rmk")
demo.add_population(initial_size = 9518, name = "admixa")
demo.add_population_split(time=130, ancestral="admix", derived=["J"])
demo.add_population_split(time=130, ancestral="admixa", derived=["Rrllrl"])
demo.add_population_split(time=219, ancestral="admixd", derived=["G"])
demo.add_population_split(time=28, ancestral="admixj", derived=["L"])
demo.add_population_split(time=219, ancestral="admixv", derived=["C"])
demo.add_population_split(time=28, ancestral="admixzt", derived=["D"])
demo.add_population_split(time=482, ancestral="admixzy", derived=["Rrm"])
demo.add_population_split(time=668, ancestral="R", derived=["Rr","Rm"])
demo.add_population_split(time=505, ancestral="Rl", derived=["I","Rlb"])
demo.add_population_split(time=482, ancestral="Rlb", derived=["Rlz"])
demo.add_population_split(time=366, ancestral="Rlz", derived=["E","Rlzn"])
demo.add_population_split(time=596, ancestral="Rm", derived=["Rl","Rmk"])
demo.add_population_split(time=596, ancestral="Rr", derived=["Rrl"])
demo.add_population_split(time=505, ancestral="Rrl", derived=["Rrlj","Rrlw"])
demo.add_population_split(time=482, ancestral="Rrlj", derived=["Rrlr"])
demo.add_population_split(time=325, ancestral="Rrll", derived=["Rrlll","Rrllr"])
demo.add_population_split(time=219, ancestral="Rrlll", derived=["B"])
demo.add_population_split(time=219, ancestral="Rrllr", derived=["Rrllrr"])
demo.add_population_split(time=28, ancestral="Rrllrl", derived=["M","K"])
demo.add_population_split(time=130, ancestral="Rrllrr", derived=["A","F"])
demo.add_population_split(time=366, ancestral="Rrlp", derived=["Rrll"])
demo.add_population_split(time=366, ancestral="Rrlr", derived=["Rrlrl","Rrlrz"])
demo.add_population_split(time=219, ancestral="Rrlrr", derived=["H"])
demo.add_population_split(time=325, ancestral="Rrlrz", derived=["Rrlrr"])
demo.add_population_split(time=482, ancestral="Rrlw", derived=["Rrlp"])
demo.add_admixture(time=130, derived="admix", ancestral=["Rrlll","admixjy"], proportions=[0.86,0.14])
demo.add_admixture(time=130, derived="admixa", ancestral=["Rrllr","Rmk"], proportions=[0.61,0.39])
demo.add_admixture(time=219, derived="admixd", ancestral=["Rrlrz","Rrm"], proportions=[0.43,0.57])
demo.add_admixture(time=28, derived="admixj", ancestral=["admixx","Rrlp"], proportions=[0.42,0.58])
demo.add_admixture(time=366, derived="admixjy", ancestral=["Rrlj","Rrlw"], proportions=[0.37,0.63])
demo.add_admixture(time=219, derived="admixv", ancestral=["Rrlrl","Rlzn"], proportions=[0.89,0.11])
demo.add_admixture(time=130, derived="admixx", ancestral=["Rrlrr","Rlzn"], proportions=[0.55,0.45])
demo.add_admixture(time=130, derived="admixz", ancestral=["Rrlrl","Rrm"], proportions=[0.15,0.85])
demo.add_admixture(time=28, derived="admixzt", ancestral=["admixz","Rlb"], proportions=[0.13,0.87])
demo.add_admixture(time=482, derived="admixzy", ancestral=["Rmk","Rr"], proportions=[0.49,0.51])
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

numpy.savetxt('./projects/qpadm/results/simout/qpadm_sim_fixed_2.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_2.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_2.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

