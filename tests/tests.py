import sys

sys.path.insert(0, '/Users/anthonywohns/Documents/mcvean_group/relative_allele_age/relative_age')

import tsinfer
import msprime

import simulations



test_ts = simulations.msprime_simulation("hi",10,10000,10000,2e-8,2e-8)

sample = simulations.generate_samples(test_ts,"hiya",0)


print(test_ts.genotype_matrix())
print("break")
print(sample.sites_genotypes[:])

test_ts.genotype_matrix() == sample.sites_genotypes[:]