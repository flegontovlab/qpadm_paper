import numpy
import math
import msprime
import multiprocess

indnam = ["K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("K", "E", "G", "F", "B", "J", "D", "L", "A", "M", "I", "C", "H"), (1612, 1196, 1954, 799, 296, 221, 720, 1196, 221, 799, 0, 221, 82) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 7084, name = "R")
demo.add_population(initial_size = 6332, name = "Rl")
demo.add_population(initial_size = 6101, name = "Rll")
demo.add_population(initial_size = 5027, name = "K")
demo.add_population(initial_size = 4380, name = "Rllr")
demo.add_population(initial_size = 8371, name = "E")
demo.add_population(initial_size = 6188, name = "Rr")
demo.add_population(initial_size = 1171, name = "Rrl")
demo.add_population(initial_size = 9085, name = "G")
demo.add_population(initial_size = 8080, name = "Rrr")
demo.add_population(initial_size = 4283, name = "Rrrl")
demo.add_population(initial_size = 4985, name = "F")
demo.add_population(initial_size = 6586, name = "Rrrr")
demo.add_population(initial_size = 5780, name = "Rrrrl")
demo.add_population(initial_size = 3793, name = "Rrrrll")
demo.add_population(initial_size = 7773, name = "B")
demo.add_population(initial_size = 1717, name = "Rrrrr")
demo.add_population(initial_size = 4764, name = "J")
demo.add_population(initial_size = 1650, name = "Rrrlf")
demo.add_population(initial_size = 3792, name = "D")
demo.add_population(initial_size = 5010, name = "admix")
demo.add_population(initial_size = 1834, name = "Rllw")
demo.add_population(initial_size = 7944, name = "admixo")
demo.add_population(initial_size = 7707, name = "admixx")
demo.add_population(initial_size = 3768, name = "L")
demo.add_population(initial_size = 1947, name = "Rrlg")
demo.add_population(initial_size = 1514, name = "Rllh")
demo.add_population(initial_size = 8866, name = "admixi")
demo.add_population(initial_size = 6460, name = "Rlg")
demo.add_population(initial_size = 7935, name = "admixq")
demo.add_population(initial_size = 1077, name = "A")
demo.add_population(initial_size = 7316, name = "Rrlgw")
demo.add_population(initial_size = 9085, name = "admixxs")
demo.add_population(initial_size = 6947, name = "Rrlgwd")
demo.add_population(initial_size = 3759, name = "admixd")
demo.add_population(initial_size = 1189, name = "M")
demo.add_population(initial_size = 4000, name = "Rrrrld")
demo.add_population(initial_size = 8837, name = "admixa")
demo.add_population(initial_size = 3422, name = "I")
demo.add_population(initial_size = 5736, name = "Rrrrlll")
demo.add_population(initial_size = 1613, name = "C")
demo.add_population(initial_size = 8420, name = "admixu")
demo.add_population(initial_size = 4607, name = "Rrls")
demo.add_population(initial_size = 3998, name = "admixg")
demo.add_population(initial_size = 1999, name = "H")
demo.add_population_split(time=82, ancestral="admixa", derived=["I"])
demo.add_population_split(time=1080, ancestral="admixd", derived=["M"])
demo.add_population_split(time=221, ancestral="admixg", derived=["H"])
demo.add_population_split(time=1466, ancestral="admixi", derived=["Rrr"])
demo.add_population_split(time=720, ancestral="admixo", derived=["Rrrrr"])
demo.add_population_split(time=296, ancestral="admixq", derived=["A"])
demo.add_population_split(time=1466, ancestral="admixx", derived=["L"])
demo.add_population_split(time=2925, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=2681, ancestral="Rl", derived=["Rllw","Rlg"])
demo.add_population_split(time=1954, ancestral="Rll", derived=["K","Rllh"])
demo.add_population_split(time=1612, ancestral="Rllh", derived=["Rllr"])
demo.add_population_split(time=1466, ancestral="Rllr", derived=["E"])
demo.add_population_split(time=2344, ancestral="Rllw", derived=["Rll"])
demo.add_population_split(time=2681, ancestral="Rr", derived=["Rrl"])
demo.add_population_split(time=2344, ancestral="Rrl", derived=["G","Rrls"])
demo.add_population_split(time=1612, ancestral="Rrlg", derived=["Rrlgw"])
demo.add_population_split(time=1466, ancestral="Rrlgw", derived=["Rrlgwd"])
demo.add_population_split(time=1954, ancestral="Rrls", derived=["Rrlg"])
demo.add_population_split(time=1196, ancestral="Rrr", derived=["Rrrl","Rrrr"])
demo.add_population_split(time=1080, ancestral="Rrrl", derived=["F","Rrrlf"])
demo.add_population_split(time=799, ancestral="Rrrlf", derived=["D"])
demo.add_population_split(time=1080, ancestral="Rrrr", derived=["Rrrrl"])
demo.add_population_split(time=799, ancestral="Rrrrl", derived=["Rrrrll","Rrrrld"])
demo.add_population_split(time=720, ancestral="Rrrrll", derived=["B","Rrrrlll"])
demo.add_population_split(time=296, ancestral="Rrrrlll", derived=["C"])
demo.add_population_split(time=296, ancestral="Rrrrr", derived=["J"])
demo.add_admixture(time=296, derived="admix", ancestral=["Rrrlf","Rrlgwd"], proportions=[0.8,0.2])
demo.add_admixture(time=82, derived="admixa", ancestral=["admixu","Rrrrld"], proportions=[0.48,0.52])
demo.add_admixture(time=1080, derived="admixd", ancestral=["Rrlgwd","Rllr"], proportions=[0.75,0.25])
demo.add_admixture(time=221, derived="admixg", ancestral=["admix","Rrls"], proportions=[0.27,0.73])
demo.add_admixture(time=1466, derived="admixi", ancestral=["Rllh","Rr"], proportions=[0.46,0.54])
demo.add_admixture(time=720, derived="admixo", ancestral=["admixxs","Rllw"], proportions=[0.37,0.63])
demo.add_admixture(time=296, derived="admixq", ancestral=["Rrrrld","Rlg"], proportions=[0.88,0.12])
demo.add_admixture(time=221, derived="admixu", ancestral=["Rrrrr","Rrrrlll"], proportions=[0.7,0.3])
demo.add_admixture(time=1466, derived="admixx", ancestral=["Rrlg","Rlg"], proportions=[0.68,0.32])
demo.add_admixture(time=799, derived="admixxs", ancestral=["Rrrr","Rrlgw"], proportions=[0.76,0.24])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_19.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_19.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_19.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

