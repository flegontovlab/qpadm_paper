import random
import numpy as np


class SubSample(object):
    def __init__(self, inpref):
        self.inpref = inpref
        self.inds = np.loadtxt(inpref+'.ind', dtype=str, delimiter='\t')

        self.snps = np.loadtxt(inpref+'.snp', dtype=str, delimiter='\t')

        gens = []
        with open(inpref+'.geno') as f:
            lines = f.readlines()
            for ln in lines:
                gens.append(list(ln.rstrip()))
        self.genos = np.array(gens)

    @staticmethod
    def _get_pseudohaplotypes(genos):
        """Create pseudo haplotypes from genotype matrix"""
        het_pos = genos == '1'
        genos[het_pos] = random.choices(('0', '2'), k=np.sum(het_pos))

        return genos

    @staticmethod
    def _add_missing(genos):
        """Add missing data to chromosomes"""
        for col in range(genos.shape[1]):
            # 5- 95% of missing data for each individual
            missing_rate = random.randrange(5, 96)/100

            # Choose positions randomly to add missing data
            pos_to_missing = random.sample(range(genos.shape[0]), k=round(genos.shape[0] * missing_rate))

            # Introduce missing data
            genos[pos_to_missing, col] = '9'

        return genos

    @staticmethod
    def _subsample_inds(genos, inds):
        """Sample a number of individuals (n) for each population. Randomly select n between 1 to 10"""
        pops = np.unique(inds[:, 2])

        inds_include = []
        for pop in pops:
            num_inds = np.sum([p in pop for p in inds[:, 2]])
            num_to_get = random.choice(range(1, num_inds+1))
            inds_to_get = random.sample(range(num_inds), k=num_to_get)

            inds_include = inds_include + inds[inds[:, 2] == pop, 0][inds_to_get].tolist()

        selector = [ind in inds_include for ind in inds[:, 0]]

        filtered_inds = inds[selector, :]
        filtered_genos = genos[:, selector]

        return filtered_genos, filtered_inds

    def subsample(self):
        genos = self.genos.copy()

        genos_1 = self._get_pseudohaplotypes(genos)
        genos_2 = self._add_missing(genos_1)
        genos_3, inds = self._subsample_inds(genos_2, self.inds)

        return genos_3, inds, self.snps

    def save(self, genos, inds, snps, outpref=None):
        if outpref is None:
            outpref = self.inpref + ".sub"

        np.savetxt(outpref+'.ind', inds,  delimiter='\t', fmt="%s")
        np.savetxt(outpref+'.snp', snps,  delimiter='\t', fmt="%s")

        genos_to_write = [''.join(i) + '\n' for i in genos]
        with open(outpref+'.geno', 'w') as f:
            f.writelines(genos_to_write)


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description="Subsample msprime simulation output to create a more realistic dataset by:"
                                                 "(1) creating pseudohaplotypes, "
                                                 "(2) randomly adding missing data along the chromosomes (5 - 95 %), "
                                                 "(3) randomly selecting individuals for each population (1 - 10 individuals per population).")
    parser.add_argument('--in', dest='inpref', type=str, required=True, help='Prefix of input files')
    parser.add_argument('--out', dest='outpref', type=str, required=False, help='Prefix of output file (Optional). The default is <input_pref>.sub.*')
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    print("Reading files...")
    subsamp = SubSample(args.inpref)
    print("Subsampling...")
    fgeno, find, fsnp = subsamp.subsample()
    print("Saving subsampled files...")
    subsamp.save(genos=fgeno, inds=find, snps=fsnp, outpref=args.outpref)
