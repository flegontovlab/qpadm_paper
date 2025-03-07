import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("B", "H", "C", "J", "F", "G", "L", "K", "E", "A", "M", "I", "D"), (356, 520, 356, 273, 543, 0, 0, 424, 273, 356, 197, 451, 197) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 6959, name = "R")
demo.add_population(initial_size = 9410, name = "Rl")
demo.add_population(initial_size = 8908, name = "Rll")
demo.add_population(initial_size = 8983, name = "Rllr")
demo.add_population(initial_size = 2930, name = "Rllrl")
demo.add_population(initial_size = 4573, name = "B")
demo.add_population(initial_size = 3313, name = "H")
demo.add_population(initial_size = 7088, name = "Rlr")
demo.add_population(initial_size = 7673, name = "Rlrl")
demo.add_population(initial_size = 6567, name = "Rlrll")
demo.add_population(initial_size = 6309, name = "C")
demo.add_population(initial_size = 8078, name = "Rlrllr")
demo.add_population(initial_size = 3576, name = "J")
demo.add_population(initial_size = 6816, name = "F")
demo.add_population(initial_size = 1915, name = "Rr")
demo.add_population(initial_size = 6352, name = "Rrl")
demo.add_population(initial_size = 1829, name = "G")
demo.add_population(initial_size = 1066, name = "L")
demo.add_population(initial_size = 6385, name = "Rrr")
demo.add_population(initial_size = 9843, name = "Rrrs")
demo.add_population(initial_size = 5419, name = "admix")
demo.add_population(initial_size = 1652, name = "Rllrld")
demo.add_population(initial_size = 3393, name = "admixn")
demo.add_population(initial_size = 2933, name = "Rlln")
demo.add_population(initial_size = 1485, name = "admixw")
demo.add_population(initial_size = 6268, name = "admixwi")
demo.add_population(initial_size = 4723, name = "K")
demo.add_population(initial_size = 6274, name = "admixb")
demo.add_population(initial_size = 3947, name = "E")
demo.add_population(initial_size = 1622, name = "Rllrd")
demo.add_population(initial_size = 7782, name = "admixu")
demo.add_population(initial_size = 8554, name = "A")
demo.add_population(initial_size = 5920, name = "Rrrst")
demo.add_population(initial_size = 5677, name = "admixwx")
demo.add_population(initial_size = 6825, name = "Rlrln")
demo.add_population(initial_size = 2884, name = "admixq")
demo.add_population(initial_size = 4814, name = "M")
demo.add_population(initial_size = 8401, name = "Rllrlc")
demo.add_population(initial_size = 7909, name = "admixe")
demo.add_population(initial_size = 8138, name = "Rllno")
demo.add_population(initial_size = 3926, name = "I")
demo.add_population(initial_size = 2573, name = "admixem")
demo.add_population(initial_size = 6402, name = "Rllrldd")
demo.add_population(initial_size = 1955, name = "D")
demo.add_population(initial_size = 2150, name = "admixqy")
demo.add_population_split(time=46, ancestral="admix", derived=["Rrl"])
demo.add_population_split(time=356, ancestral="admixb", derived=["E"])
demo.add_population_split(time=451, ancestral="admixe", derived=["Rlrll"])
demo.add_population_split(time=273, ancestral="admixq", derived=["M"])
demo.add_population_split(time=424, ancestral="admixu", derived=["A"])
demo.add_population_split(time=520, ancestral="admixw", derived=["admixwi"])
demo.add_population_split(time=451, ancestral="admixwi", derived=["K"])
demo.add_population_split(time=669, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=618, ancestral="Rl", derived=["Rll","Rlr"])
demo.add_population_split(time=566, ancestral="Rll", derived=["Rllr","Rlln"])
demo.add_population_split(time=543, ancestral="Rlln", derived=["Rllno"])
demo.add_population_split(time=520, ancestral="Rllno", derived=["I"])
demo.add_population_split(time=543, ancestral="Rllr", derived=["H","Rllrlc"])
demo.add_population_split(time=451, ancestral="Rllrd", derived=["Rllrl"])
demo.add_population_split(time=424, ancestral="Rllrl", derived=["B","Rllrld"])
demo.add_population_split(time=520, ancestral="Rllrlc", derived=["Rllrd"])
demo.add_population_split(time=356, ancestral="Rllrld", derived=["Rllrldd"])
demo.add_population_split(time=273, ancestral="Rllrldd", derived=["D"])
demo.add_population_split(time=566, ancestral="Rlr", derived=["Rlrl","F"])
demo.add_population_split(time=543, ancestral="Rlrl", derived=["Rlrln"])
demo.add_population_split(time=424, ancestral="Rlrll", derived=["C","Rlrllr"])
demo.add_population_split(time=356, ancestral="Rlrllr", derived=["J"])
demo.add_population_split(time=618, ancestral="Rr", derived=["Rrr"])
demo.add_population_split(time=21, ancestral="Rrl", derived=["G","L"])
demo.add_population_split(time=566, ancestral="Rrr", derived=["Rrrs"])
demo.add_population_split(time=543, ancestral="Rrrs", derived=["Rrrst"])
demo.add_admixture(time=46, derived="admix", ancestral=["admixqy","Rr"], proportions=[0.47,0.53])
demo.add_admixture(time=356, derived="admixb", ancestral=["admixwi","admixwx"], proportions=[0.41,0.59])
demo.add_admixture(time=451, derived="admixe", ancestral=["Rlrln","Rllrlc"], proportions=[0.61,0.39])
demo.add_admixture(time=273, derived="admixem", ancestral=["Rllrld","Rllno"], proportions=[0.58,0.42])
demo.add_admixture(time=197, derived="admixn", ancestral=["admixem","Rrrs"], proportions=[0.26,0.74])
demo.add_admixture(time=273, derived="admixq", ancestral=["Rlrllr","Rlrln"], proportions=[0.43,0.57])
demo.add_admixture(time=143, derived="admixqy", ancestral=["admixn","Rllrldd"], proportions=[0.48,0.52])
demo.add_admixture(time=424, derived="admixu", ancestral=["Rllrd","Rrrst"], proportions=[0.46,0.54])
demo.add_admixture(time=520, derived="admixw", ancestral=["Rlln","Rrr"], proportions=[0.18,0.82])
demo.add_admixture(time=451, derived="admixwx", ancestral=["Rrrst","Rlrl"], proportions=[0.64,0.36])
demo.sort_events()

r_chrom = 2e-08
r_break = math.log(2)
chrom_positions = [0, 1e+08, 2e+08, 3e+08, 4e+08, 5e+08, 6e+08, 7e+08, 8e+08, 9e+08, 1e+09]
map_positions = [0, 1e+08, 100000001, 2e+08, 200000001, 3e+08, 300000001, 4e+08, 400000001, 5e+08, 500000001, 6e+08, 600000001, 7e+08, 700000001, 8e+08, 800000001, 9e+08, 900000001, 1e+09]
rates = [r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom, r_break, r_chrom]
rate_map = msprime.RateMap(position=map_positions, rate=rates)

ts = msprime.sim_ancestry(samples=samples, model=[msprime.DiscreteTimeWrightFisher(duration=25), msprime.StandardCoalescent()], demography=demo, recombination_rate=rate_map)
tree_sequence = msprime.sim_mutations(ts, rate=1.25e-08, model='binary')
hap_gt = tree_sequence.genotype_matrix()
gt = hap_gt[:, range(0,nhap,2)] + hap_gt[:, range(1,nhap,2)]
nsnps = gt.shape[0]
ts_chroms = numpy.searchsorted(numpy.array(chrom_positions[1:]), tree_sequence.tables.sites.position, side='right') + 1

numpy.savetxt('./projects/qpadm/results/simout/10chr/sim_10/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/10chr/sim_10/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/10chr/sim_10/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

