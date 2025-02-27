import numpy
import math
import msprime
import multiprocess

indnam = ["J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("J", "K", "L", "F", "G", "D", "I", "E", "A", "C", "B", "M", "H"), (1669, 668, 1972, 2265, 904, 1166, 521, 0, 904, 521, 278, 278, 165) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 4621, name = "R")
demo.add_population(initial_size = 8185, name = "Rr")
demo.add_population(initial_size = 1550, name = "Rrl")
demo.add_population(initial_size = 1859, name = "Rrll")
demo.add_population(initial_size = 1128, name = "Rrlll")
demo.add_population(initial_size = 3224, name = "Rrllll")
demo.add_population(initial_size = 7686, name = "J")
demo.add_population(initial_size = 2525, name = "Rrlllr")
demo.add_population(initial_size = 8862, name = "Rrlllrl")
demo.add_population(initial_size = 7822, name = "Rrlllrll")
demo.add_population(initial_size = 1878, name = "Rrlllrlr")
demo.add_population(initial_size = 7496, name = "K")
demo.add_population(initial_size = 2566, name = "Rrlr")
demo.add_population(initial_size = 1811, name = "L")
demo.add_population(initial_size = 2403, name = "Rrr")
demo.add_population(initial_size = 2414, name = "F")
demo.add_population(initial_size = 7720, name = "Rrllly")
demo.add_population(initial_size = 9292, name = "admix")
demo.add_population(initial_size = 5483, name = "G")
demo.add_population(initial_size = 9318, name = "Rrllllp")
demo.add_population(initial_size = 3272, name = "D")
demo.add_population(initial_size = 6870, name = "admixn")
demo.add_population(initial_size = 7075, name = "Rrlllq")
demo.add_population(initial_size = 4944, name = "admixa")
demo.add_population(initial_size = 3168, name = "admixb")
demo.add_population(initial_size = 5024, name = "I")
demo.add_population(initial_size = 1867, name = "Rrlln")
demo.add_population(initial_size = 3270, name = "Rrlly")
demo.add_population(initial_size = 9745, name = "admixh")
demo.add_population(initial_size = 1826, name = "E")
demo.add_population(initial_size = 1044, name = "Rrllno")
demo.add_population(initial_size = 2036, name = "admixm")
demo.add_population(initial_size = 1114, name = "Rrlllrlv")
demo.add_population(initial_size = 2556, name = "admixi")
demo.add_population(initial_size = 4916, name = "A")
demo.add_population(initial_size = 3116, name = "Rrlllrlrx")
demo.add_population(initial_size = 2726, name = "C")
demo.add_population(initial_size = 7955, name = "admixp")
demo.add_population(initial_size = 1624, name = "B")
demo.add_population(initial_size = 6054, name = "Rrllnoe")
demo.add_population(initial_size = 9898, name = "M")
demo.add_population(initial_size = 5017, name = "admixl")
demo.add_population(initial_size = 9823, name = "H")
demo.add_population(initial_size = 8887, name = "Rrllllc")
demo.add_population(initial_size = 3921, name = "admixe")
demo.add_population_split(time=1166, ancestral="admix", derived=["G"])
demo.add_population_split(time=668, ancestral="admixb", derived=["I"])
demo.add_population_split(time=1571, ancestral="admixe", derived=["Rrlly"])
demo.add_population_split(time=165, ancestral="admixh", derived=["E"])
demo.add_population_split(time=1166, ancestral="admixi", derived=["A"])
demo.add_population_split(time=278, ancestral="admixl", derived=["H"])
demo.add_population_split(time=1166, ancestral="admixn", derived=["Rrlllrlr"])
demo.add_population_split(time=521, ancestral="admixp", derived=["B"])
demo.add_population_split(time=2636, ancestral="R", derived=["Rr"])
demo.add_population_split(time=2520, ancestral="Rr", derived=["Rrl","Rrr"])
demo.add_population_split(time=2362, ancestral="Rrl", derived=["Rrll","Rrlr"])
demo.add_population_split(time=2265, ancestral="Rrll", derived=["Rrlll"])
demo.add_population_split(time=1972, ancestral="Rrlll", derived=["Rrllly","Rrlllq"])
demo.add_population_split(time=1762, ancestral="Rrllll", derived=["J","Rrllllc"])
demo.add_population_split(time=1669, ancestral="Rrllllc", derived=["Rrllllp"])
demo.add_population_split(time=1571, ancestral="Rrllllp", derived=["D"])
demo.add_population_split(time=1886, ancestral="Rrlllq", derived=["Rrlllr"])
demo.add_population_split(time=1762, ancestral="Rrlllr", derived=["Rrlllrl","Rrlllrlv"])
demo.add_population_split(time=1669, ancestral="Rrlllrl", derived=["Rrlllrll"])
demo.add_population_split(time=904, ancestral="Rrlllrlr", derived=["K","Rrlllrlrx"])
demo.add_population_split(time=668, ancestral="Rrlllrlrx", derived=["C"])
demo.add_population_split(time=1886, ancestral="Rrllly", derived=["Rrllll"])
demo.add_population_split(time=904, ancestral="Rrlln", derived=["Rrllno"])
demo.add_population_split(time=668, ancestral="Rrllno", derived=["Rrllnoe"])
demo.add_population_split(time=521, ancestral="Rrllnoe", derived=["M"])
demo.add_population_split(time=1166, ancestral="Rrlly", derived=["Rrlln"])
demo.add_population_split(time=2265, ancestral="Rrlr", derived=["L"])
demo.add_population_split(time=2362, ancestral="Rrr", derived=["F"])
demo.add_admixture(time=1166, derived="admix", ancestral=["admixa","Rrllly"], proportions=[0.16,0.84])
demo.add_admixture(time=1762, derived="admixa", ancestral=["Rrlllq","Rrr"], proportions=[0.44,0.56])
demo.add_admixture(time=668, derived="admixb", ancestral=["Rrlln","Rrlr"], proportions=[0.41,0.59])
demo.add_admixture(time=1571, derived="admixe", ancestral=["Rrllllc","Rrll"], proportions=[0.67,0.33])
demo.add_admixture(time=165, derived="admixh", ancestral=["admixm","Rrlly"], proportions=[0.22,0.78])
demo.add_admixture(time=1166, derived="admixi", ancestral=["Rrlllrll","Rrlllrlv"], proportions=[0.37,0.63])
demo.add_admixture(time=278, derived="admixl", ancestral=["Rrllnoe","R"], proportions=[0.77,0.23])
demo.add_admixture(time=278, derived="admixm", ancestral=["Rrllno","Rrlllrll"], proportions=[0.83,0.17])
demo.add_admixture(time=1166, derived="admixn", ancestral=["Rrllllp","Rrlllrl"], proportions=[0.28,0.72])
demo.add_admixture(time=521, derived="admixp", ancestral=["Rrlllrlrx","Rrlllrlv"], proportions=[0.17,0.83])
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

numpy.savetxt('./projects/qpadm/results/simout/deep/qpadm_sim_20.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_20.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/deep/qpadm_sim_20.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

