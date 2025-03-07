import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("L", "M", "I", "E", "H", "A", "D", "B", "F", "K", "C", "G", "J"), (375, 375, 423, 423, 646, 540, 375, 423, 423, 0, 423, 197, 197) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 8285, name = "Rlll")
demo.add_population(initial_size = 2511, name = "L")
demo.add_population(initial_size = 8408, name = "M")
demo.add_population(initial_size = 9113, name = "Rlr")
demo.add_population(initial_size = 5808, name = "Rlrl")
demo.add_population(initial_size = 8084, name = "I")
demo.add_population(initial_size = 1255, name = "E")
demo.add_population(initial_size = 2714, name = "R")
demo.add_population(initial_size = 9249, name = "Rr")
demo.add_population(initial_size = 3732, name = "Rrl")
demo.add_population(initial_size = 8967, name = "H")
demo.add_population(initial_size = 8993, name = "Rrr")
demo.add_population(initial_size = 9925, name = "Rrrl")
demo.add_population(initial_size = 7023, name = "A")
demo.add_population(initial_size = 3895, name = "Rrrlr")
demo.add_population(initial_size = 8892, name = "D")
demo.add_population(initial_size = 7878, name = "Rrrr")
demo.add_population(initial_size = 1156, name = "B")
demo.add_population(initial_size = 9704, name = "F")
demo.add_population(initial_size = 2197, name = "Rrrlz")
demo.add_population(initial_size = 1405, name = "admix")
demo.add_population(initial_size = 7508, name = "K")
demo.add_population(initial_size = 3226, name = "Rll")
demo.add_population(initial_size = 3271, name = "Rllj")
demo.add_population(initial_size = 7048, name = "C")
demo.add_population(initial_size = 3574, name = "admixx")
demo.add_population(initial_size = 7172, name = "Rrrlrg")
demo.add_population(initial_size = 7581, name = "admixq")
demo.add_population(initial_size = 7787, name = "Rlle")
demo.add_population(initial_size = 5509, name = "admixk")
demo.add_population(initial_size = 5230, name = "G")
demo.add_population(initial_size = 9550, name = "Rp")
demo.add_population(initial_size = 5590, name = "Rl")
demo.add_population(initial_size = 3067, name = "admixj")
demo.add_population(initial_size = 4403, name = "Rrc")
demo.add_population(initial_size = 5994, name = "admixc")
demo.add_population(initial_size = 2830, name = "J")
demo.add_population(initial_size = 8244, name = "Rllw")
demo.add_population(initial_size = 8223, name = "admixw")
demo.add_population(initial_size = 4453, name = "Rrrlzk")
demo.add_population(initial_size = 8387, name = "admixo")
demo.add_population(initial_size = 6900, name = "Rlc")
demo.add_population(initial_size = 8104, name = "admixn")
demo.add_population(initial_size = 2892, name = "Rrrx")
demo.add_population(initial_size = 4182, name = "admixkd")
demo.add_population_split(time=109, ancestral="admix", derived=["K"])
demo.add_population_split(time=241, ancestral="admixc", derived=["J"])
demo.add_population_split(time=241, ancestral="admixk", derived=["G"])
demo.add_population_split(time=791, ancestral="R", derived=["Rr","Rp"])
demo.add_population_split(time=683, ancestral="Rl", derived=["Rllw","Rlc"])
demo.add_population_split(time=646, ancestral="Rlc", derived=["Rlr"])
demo.add_population_split(time=592, ancestral="Rll", derived=["Rllj","Rlle"])
demo.add_population_split(time=540, ancestral="Rlle", derived=["Rlll"])
demo.add_population_split(time=540, ancestral="Rllj", derived=["C"])
demo.add_population_split(time=423, ancestral="Rlll", derived=["L","M"])
demo.add_population_split(time=646, ancestral="Rllw", derived=["Rll"])
demo.add_population_split(time=592, ancestral="Rlr", derived=["Rlrl"])
demo.add_population_split(time=540, ancestral="Rlrl", derived=["I","E"])
demo.add_population_split(time=751, ancestral="Rp", derived=["Rl"])
demo.add_population_split(time=751, ancestral="Rr", derived=["Rrl","Rrc"])
demo.add_population_split(time=683, ancestral="Rrc", derived=["Rrr"])
demo.add_population_split(time=683, ancestral="Rrl", derived=["H"])
demo.add_population_split(time=646, ancestral="Rrr", derived=["Rrrl","Rrrx"])
demo.add_population_split(time=592, ancestral="Rrrl", derived=["A","Rrrlz"])
demo.add_population_split(time=423, ancestral="Rrrlr", derived=["D","Rrrlrg"])
demo.add_population_split(time=540, ancestral="Rrrlz", derived=["Rrrlr","Rrrlzk"])
demo.add_population_split(time=540, ancestral="Rrrr", derived=["B","F"])
demo.add_population_split(time=592, ancestral="Rrrx", derived=["Rrrr"])
demo.add_admixture(time=109, derived="admix", ancestral=["admixx","Rrrlzk"], proportions=[0.42,0.58])
demo.add_admixture(time=241, derived="admixc", ancestral=["admixj","Rrc"], proportions=[0.28,0.72])
demo.add_admixture(time=338, derived="admixj", ancestral=["Rrrlrg","Rp"], proportions=[0.82,0.18])
demo.add_admixture(time=241, derived="admixk", ancestral=["admixn","Rlle"], proportions=[0.44,0.56])
demo.add_admixture(time=540, derived="admixkd", ancestral=["Rrrx","Rllw"], proportions=[0.88,0.12])
demo.add_admixture(time=375, derived="admixn", ancestral=["admixw","Rlc"], proportions=[0.42,0.58])
demo.add_admixture(time=241, derived="admixo", ancestral=["admixq","Rrrlzk"], proportions=[0.21,0.79])
demo.add_admixture(time=338, derived="admixq", ancestral=["Rrrlrg","Rllj"], proportions=[0.76,0.24])
demo.add_admixture(time=423, derived="admixw", ancestral=["admixkd","Rrl"], proportions=[0.3,0.7])
demo.add_admixture(time=197, derived="admixx", ancestral=["admixo","Rlr"], proportions=[0.17,0.83])
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

numpy.savetxt('./projects/qpadm/results/simout/3chr/sim_34/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/3chr/sim_34/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/3chr/sim_34/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

