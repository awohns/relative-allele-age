"""
Implements simulations with or without error.
"""

import numpy as np
import math
import msprime
import tsinfer
import logging
import warnings
warnings.filterwarnings("ignore",message="numpy.dtype size changed")

import geva_compare
import freq_compare



def msprime_simulation(name,sample_size, Ne, length, recombination_rate, mutation_rate):
    #Step one: simulate dataset with msprime
    simulation = msprime.simulate(sample_size=sample_size, Ne=Ne, length=length, recombination_rate=recombination_rate, mutation_rate=mutation_rate)
    with open("../tmp/"+name+".vcf", "w") as vcf_file:
        simulation.write_vcf(vcf_file, 2)
    return(simulation)

def make_errors(g, p):
    """
    Resample for a specific variant, whose genotypes g are passed in
    For each sample an error occurs with probability p. Errors are generated by
    sampling values from the stationary distribution, that is, if we have an
    allele frequency of f, a 1 is emitted with probability f and a
    0 with probability 1 - f. Thus, there is a possibility that an 'error'
    will in fact result in the same value.
    """
    w = np.copy(g)
    if p > 0:
        m = g.shape[0]
        frequency = np.sum(g) / m
        # Randomly choose samples with probability p
        samples = np.where(np.random.random(m) < p)[0]
        # Generate observations from the stationary distribution.
        errors = (np.random.random(samples.shape[0]) < frequency).astype(int)
        w[samples] = errors
    return w

def generate_samples(ts, filename, real_error_rate=0):
    """
    Generate a samples file from a simulated ts
    Samples may have bits flipped with a specified probability.
    (reject any variants that result in a fixed column)
    """
    if real_error_rate > 0:
        logging.debug("converting real error rate to an error param by multiplying by log(n)")
    error_param = real_error_rate * math.log(ts.num_samples)
    record_rate = logging.getLogger().isEnabledFor(logging.INFO)
    n_variants = bits_flipped = 0
    assert ts.num_sites != 0
    sample_data = tsinfer.SampleData(path="../tmp/"+filename + ".samples", sequence_length=ts.sequence_length)
    for v in ts.variants():
        n_variants += 1
        if error_param <=0:
            sample_data.add_site(
                position=v.site.position, alleles=v.alleles,
                genotypes=v.genotypes)
        else:
            #make new genotypes with error
            # Reject any columns that have no 1s or no zeros.
            # Unless the original also has them, as occasionally we have
            # some sims (e.g. under selection) where a variant is fixed
            while True:
                genotypes = make_errors(v.genotypes, error_param)
                s = np.sum(genotypes)
                if 0 < s < ts.sample_size:
                    break
                if s == np.sum(v.genotypes):
                    break
            if record_rate:
                bits_flipped += np.sum(np.logical_xor(genotypes, v.genotypes))
            sample_data.add_site(
                position=v.site.position, alleles=v.alleles,
                genotypes=genotypes)
    if real_error_rate>0:
        logging.info("Error of {} injected into {}".format(real_error_rate, os.path.basename(filename))
            + ": actual error rate = {} (error param = {})".format(
                bits_flipped/(n_variants*ts.sample_size), error_param) if record_rate else "")
    sample_data.finalise()
    return sample_data



