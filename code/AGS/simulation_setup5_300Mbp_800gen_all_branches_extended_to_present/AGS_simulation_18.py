import numpy
import math
import msprime
import multiprocess

indnam = ["D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("D", "H", "G", "F", "M", "E", "L", "B", "A", "I", "J", "K", "C"), (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 7293, name = "R")
demo.add_population(initial_size = 2015, name = "Rl")
demo.add_population(initial_size = 7947, name = "Rll")
demo.add_population(initial_size = 4451, name = "Rlll")
demo.add_population(initial_size = 2152, name = "D")
demo.add_population(initial_size = 4165, name = "H")
demo.add_population(initial_size = 4731, name = "G")
demo.add_population(initial_size = 3977, name = "Rlr")
demo.add_population(initial_size = 5535, name = "Rlrl")
demo.add_population(initial_size = 9184, name = "Rlrll")
demo.add_population(initial_size = 4619, name = "F")
demo.add_population(initial_size = 6878, name = "Rlrr")
demo.add_population(initial_size = 8788, name = "Rlrrr")
demo.add_population(initial_size = 9161, name = "M")
demo.add_population(initial_size = 2972, name = "Rr")
demo.add_population(initial_size = 7526, name = "Rrl")
demo.add_population(initial_size = 8757, name = "E")
demo.add_population(initial_size = 5830, name = "Rrr")
demo.add_population(initial_size = 4795, name = "L")
demo.add_population(initial_size = 9208, name = "B")
demo.add_population(initial_size = 1784, name = "Rlp")
demo.add_population(initial_size = 8631, name = "admix")
demo.add_population(initial_size = 3408, name = "A")
demo.add_population(initial_size = 7882, name = "Rrlt")
demo.add_population(initial_size = 5716, name = "admixk")
demo.add_population(initial_size = 9535, name = "Rlla")
demo.add_population(initial_size = 7540, name = "admixv")
demo.add_population(initial_size = 9319, name = "Rlo")
demo.add_population(initial_size = 2321, name = "admixkm")
demo.add_population(initial_size = 3179, name = "I")
demo.add_population(initial_size = 3849, name = "admixkmf")
demo.add_population(initial_size = 8633, name = "admixz")
demo.add_population(initial_size = 6215, name = "Rrlti")
demo.add_population(initial_size = 6894, name = "J")
demo.add_population(initial_size = 4473, name = "admixo")
demo.add_population(initial_size = 2777, name = "K")
demo.add_population(initial_size = 5088, name = "Rrltiy")
demo.add_population(initial_size = 1116, name = "admixy")
demo.add_population(initial_size = 9074, name = "C")
demo.add_population(initial_size = 1961, name = "Rrltiv")
demo.add_population(initial_size = 2246, name = "admixvj")
demo.add_population(initial_size = 6958, name = "Rrltig")
demo.add_population(initial_size = 4077, name = "admixl")
demo.add_population(initial_size = 4778, name = "Rrltivf")
demo.add_population(initial_size = 3568, name = "admixvn")
demo.add_population_split(time=135, ancestral="admix", derived=["A"])
demo.add_population_split(time=207, ancestral="admixk", derived=["admixkmf"])
demo.add_population_split(time=42, ancestral="admixkm", derived=["I"])
demo.add_population_split(time=181, ancestral="admixo", derived=["K"])
demo.add_population_split(time=207, ancestral="admixvn", derived=["Rlll"])
demo.add_population_split(time=135, ancestral="admixy", derived=["C"])
demo.add_population_split(time=705, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=665, ancestral="Rl", derived=["Rll","Rlo"])
demo.add_population_split(time=583, ancestral="Rll", derived=["G"])
demo.add_population_split(time=557, ancestral="Rlla", derived=["Rlp"])
demo.add_population_split(time=181, ancestral="Rlll", derived=["D","H"])
demo.add_population_split(time=583, ancestral="Rlo", derived=["Rlla"])
demo.add_population_split(time=485, ancestral="Rlp", derived=["Rlr"])
demo.add_population_split(time=455, ancestral="Rlr", derived=["Rlrl","Rlrr"])
demo.add_population_split(time=350, ancestral="Rlrl", derived=["Rlrll","F"])
demo.add_population_split(time=350, ancestral="Rlrr", derived=["Rlrrr"])
demo.add_population_split(time=265, ancestral="Rlrrr", derived=["M"])
demo.add_population_split(time=665, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=583, ancestral="Rrl", derived=["E","Rrlt"])
demo.add_population_split(time=557, ancestral="Rrlt", derived=["Rrlti"])
demo.add_population_split(time=485, ancestral="Rrlti", derived=["J","Rrltig"])
demo.add_population_split(time=455, ancestral="Rrltig", derived=["Rrltiv"])
demo.add_population_split(time=350, ancestral="Rrltiv", derived=["Rrltiy","Rrltivf"])
demo.add_population_split(time=583, ancestral="Rrr", derived=["L","B"])
demo.add_admixture(time=135, derived="admix", ancestral=["admixkmf","Rlp"], proportions=[0.3,0.7])
demo.add_admixture(time=207, derived="admixk", ancestral=["admixv","Rrlt"], proportions=[0.31,0.69])
demo.add_admixture(time=42, derived="admixkm", ancestral=["admixz","Rlo"], proportions=[0.46,0.54])
demo.add_admixture(time=207, derived="admixl", ancestral=["Rrltiy","Rrltig"], proportions=[0.71,0.29])
demo.add_admixture(time=181, derived="admixo", ancestral=["admixvj","Rrltiy"], proportions=[0.39,0.61])
demo.add_admixture(time=265, derived="admixv", ancestral=["Rlrr","Rlla"], proportions=[0.52,0.48])
demo.add_admixture(time=207, derived="admixvj", ancestral=["Rlrll","Rrltivf"], proportions=[0.11,0.89])
demo.add_admixture(time=207, derived="admixvn", ancestral=["Rrltivf","Rll"], proportions=[0.73,0.27])
demo.add_admixture(time=135, derived="admixy", ancestral=["admixl","Rlrrr"], proportions=[0.35,0.65])
demo.add_admixture(time=135, derived="admixz", ancestral=["admixkmf","Rlrll"], proportions=[0.34,0.66])
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

numpy.savetxt('./projects/qpadm/results/simout/qpadm_sim_fixed_18.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_18.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/qpadm_sim_fixed_18.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

