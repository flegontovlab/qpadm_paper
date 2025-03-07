import numpy
import math
import sys
import msprime
import multiprocess

if len(sys.argv)  != 2:
	raise ValueError("Need to specify simulation no")

sim_id = str(sys.argv[1])

indnam = ["B_1", "B_2", "B_3", "B_4", "B_5", "B_6", "B_7", "B_8", "B_9", "B_10", "E_1", "E_2", "E_3", "E_4", "E_5", "E_6", "E_7", "E_8", "E_9", "E_10", "I_1", "I_2", "I_3", "I_4", "I_5", "I_6", "I_7", "I_8", "I_9", "I_10", "L_1", "L_2", "L_3", "L_4", "L_5", "L_6", "L_7", "L_8", "L_9", "L_10", "G_1", "G_2", "G_3", "G_4", "G_5", "G_6", "G_7", "G_8", "G_9", "G_10", "D_1", "D_2", "D_3", "D_4", "D_5", "D_6", "D_7", "D_8", "D_9", "D_10", "J_1", "J_2", "J_3", "J_4", "J_5", "J_6", "J_7", "J_8", "J_9", "J_10", "M_1", "M_2", "M_3", "M_4", "M_5", "M_6", "M_7", "M_8", "M_9", "M_10", "K_1", "K_2", "K_3", "K_4", "K_5", "K_6", "K_7", "K_8", "K_9", "K_10", "C_1", "C_2", "C_3", "C_4", "C_5", "C_6", "C_7", "C_8", "C_9", "C_10", "F_1", "F_2", "F_3", "F_4", "F_5", "F_6", "F_7", "F_8", "F_9", "F_10", "A_1", "A_2", "A_3", "A_4", "A_5", "A_6", "A_7", "A_8", "A_9", "A_10", "H_1", "H_2", "H_3", "H_4", "H_5", "H_6", "H_7", "H_8", "H_9", "H_10"]

samples = [msprime.SampleSet(num_samples=i, population=j, time=k) for i, j, k in zip( (10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10), ("B", "E", "I", "L", "G", "D", "J", "M", "K", "C", "F", "A", "H"), (716, 364, 364, 364, 0, 0, 561, 716, 219, 219, 443, 467, 257) )]
nhap = int(len(indnam) * 2)

demo = msprime.Demography()
demo.add_population(initial_size = 4673, name = "R")
demo.add_population(initial_size = 2807, name = "Rl")
demo.add_population(initial_size = 7615, name = "Rll")
demo.add_population(initial_size = 8227, name = "B")
demo.add_population(initial_size = 8654, name = "Rr")
demo.add_population(initial_size = 5119, name = "Rrl")
demo.add_population(initial_size = 5714, name = "Rrll")
demo.add_population(initial_size = 8888, name = "Rrlll")
demo.add_population(initial_size = 2427, name = "Rrllll")
demo.add_population(initial_size = 2154, name = "Rrlllll")
demo.add_population(initial_size = 2962, name = "E")
demo.add_population(initial_size = 3552, name = "I")
demo.add_population(initial_size = 1130, name = "Rrlllr")
demo.add_population(initial_size = 1487, name = "L")
demo.add_population(initial_size = 7989, name = "Rrllr")
demo.add_population(initial_size = 2046, name = "G")
demo.add_population(initial_size = 1962, name = "D")
demo.add_population(initial_size = 5868, name = "Rrlr")
demo.add_population(initial_size = 8722, name = "J")
demo.add_population(initial_size = 8626, name = "M")
demo.add_population(initial_size = 1390, name = "Rrlllc")
demo.add_population(initial_size = 3928, name = "admix")
demo.add_population(initial_size = 9020, name = "Rrls")
demo.add_population(initial_size = 6801, name = "admixq")
demo.add_population(initial_size = 9103, name = "K")
demo.add_population(initial_size = 8771, name = "Rrlllri")
demo.add_population(initial_size = 8269, name = "admixy")
demo.add_population(initial_size = 1393, name = "C")
demo.add_population(initial_size = 6681, name = "Rrlsu")
demo.add_population(initial_size = 3547, name = "admixn")
demo.add_population(initial_size = 9521, name = "admixu")
demo.add_population(initial_size = 5027, name = "Rrlrt")
demo.add_population(initial_size = 7230, name = "Rrlsm")
demo.add_population(initial_size = 6449, name = "admixym")
demo.add_population(initial_size = 8609, name = "F")
demo.add_population(initial_size = 4835, name = "Rrlsf")
demo.add_population(initial_size = 2850, name = "admixw")
demo.add_population(initial_size = 6460, name = "A")
demo.add_population(initial_size = 3715, name = "Rrlro")
demo.add_population(initial_size = 9125, name = "admixc")
demo.add_population(initial_size = 5964, name = "H")
demo.add_population(initial_size = 9392, name = "Rrlsug")
demo.add_population(initial_size = 8948, name = "admixm")
demo.add_population(initial_size = 7775, name = "Rrlsj")
demo.add_population(initial_size = 5109, name = "admixqr")
demo.add_population_split(time=364, ancestral="admixc", derived=["H"])
demo.add_population_split(time=103, ancestral="admixm", derived=["Rrllr"])
demo.add_population_split(time=257, ancestral="admixq", derived=["K"])
demo.add_population_split(time=509, ancestral="admixw", derived=["A"])
demo.add_population_split(time=257, ancestral="admixy", derived=["C"])
demo.add_population_split(time=467, ancestral="admixym", derived=["F"])
demo.add_population_split(time=783, ancestral="R", derived=["Rl","Rr"])
demo.add_population_split(time=759, ancestral="Rl", derived=["Rll","B"])
demo.add_population_split(time=759, ancestral="Rr", derived=["Rrl","M"])
demo.add_population_split(time=716, ancestral="Rrl", derived=["Rrlr","Rrls"])
demo.add_population_split(time=509, ancestral="Rrll", derived=["Rrlll"])
demo.add_population_split(time=467, ancestral="Rrlll", derived=["Rrllll","Rrlllc"])
demo.add_population_split(time=443, ancestral="Rrlllc", derived=["Rrlllr"])
demo.add_population_split(time=443, ancestral="Rrllll", derived=["Rrlllll"])
demo.add_population_split(time=411, ancestral="Rrlllll", derived=["E","I"])
demo.add_population_split(time=411, ancestral="Rrlllr", derived=["L","Rrlllri"])
demo.add_population_split(time=57, ancestral="Rrllr", derived=["G","D"])
demo.add_population_split(time=646, ancestral="Rrlr", derived=["J","Rrlro"])
demo.add_population_split(time=561, ancestral="Rrlro", derived=["Rrlrt"])
demo.add_population_split(time=646, ancestral="Rrls", derived=["Rrlsf","Rrlsj"])
demo.add_population_split(time=561, ancestral="Rrlsf", derived=["Rrlsm"])
demo.add_population_split(time=561, ancestral="Rrlsj", derived=["Rrll"])
demo.add_population_split(time=509, ancestral="Rrlsm", derived=["Rrlsu"])
demo.add_population_split(time=467, ancestral="Rrlsu", derived=["Rrlsug"])
demo.add_admixture(time=411, derived="admix", ancestral=["Rrlllc","admixu"], proportions=[0.64,0.36])
demo.add_admixture(time=364, derived="admixc", ancestral=["admixn","Rrlro"], proportions=[0.46,0.54])
demo.add_admixture(time=103, derived="admixm", ancestral=["admix","admixqr"], proportions=[0.12,0.88])
demo.add_admixture(time=411, derived="admixn", ancestral=["Rrlsug","Rll"], proportions=[0.65,0.35])
demo.add_admixture(time=257, derived="admixq", ancestral=["Rrlllri","Rrlsu"], proportions=[0.44,0.56])
demo.add_admixture(time=411, derived="admixqr", ancestral=["Rrlsug","Rrlsj"], proportions=[0.9,0.1])
demo.add_admixture(time=467, derived="admixu", ancestral=["Rrll","Rrlrt"], proportions=[0.27,0.73])
demo.add_admixture(time=509, derived="admixw", ancestral=["Rrlsf","Rll"], proportions=[0.64,0.36])
demo.add_admixture(time=257, derived="admixy", ancestral=["Rrlllri","Rrllll"], proportions=[0.53,0.47])
demo.add_admixture(time=467, derived="admixym", ancestral=["Rrlrt","Rrlsm"], proportions=[0.12,0.88])
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

numpy.savetxt('./projects/qpadm/results/simout/3chr/sim_7/qpadm_sim_' + sim_id + '.geno', gt, '%d', '')

with open('./projects/qpadm/results/simout/3chr/sim_7/qpadm_sim_' + sim_id + '.snp', 'w') as f:
  for i, j in enumerate(tree_sequence.tables.sites.position):
    chrm = int(ts_chroms[i])
    if chrm == 1:
      pos = int(j) 
    else:
      pos = int(j - chrom_positions[chrm-1])
    nsnps_chr = numpy.sum(ts_chroms == chrm)
    pos_chr = numpy.sum(ts_chroms[:i] == chrm)
    bytes = f.write('rs'+str(i+1)+'\t' + str(chrm) + '\t' + str(pos_chr*50/float(nsnps_chr-1)) + '\t' + str(pos) + '\tA\tC\n')

with open('./projects/qpadm/results/simout/3chr/sim_7/qpadm_sim_' + sim_id + '.ind', 'w') as f:
  for i in range(len(indnam)):
    bytes = f.write(indnam[i] + '\tU\t' + indnam[i].rstrip("1234567890").rstrip("_") + '\n')

