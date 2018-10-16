"""

Inject error (at empirically determined rates) into a simulated tree sequence.

Output is a ts-kit samples file. 

"""



import logging

import random

import numpy as np

import pandas as pd


import tsinfer



def make_errors_genotype_model(g, error_matrix):

    """

    Given an empirically estimated error matrix, resample for a particular

    variant. Given a true genotype of g0, g1, or g2, return observed genotype

    depending on the error_matrix. For example, given a variant with genotype

    g0, return g0 with probability 99.942%, g1 with probability 0.041%, and

    g2 with probability 0.017%. Treat each pair of alleles as a diploid 

    individual. 

    """

    w = np.copy(g)

    #Make diploid (iterate each pair of alleles)
    genos=[(w[i],w[i+1]) for i in range(0,w.shape[0],2)]
    
    #Record the true genotypes

    g0 = [i for i, x in enumerate(genos) if x == (0,0)]
    g1a = [i for i, x in enumerate(genos) if x == (1,0)]
    g1b = [i for i, x in enumerate(genos) if x == (0,1)]
    g2 = [i for i, x in enumerate(genos) if x == (1,1)]

    for idx in g0:

        result=random.choices(
            [(0,0),(1,0),(1,1)], weights=error_matrix[['p00','p01','p02']].values[0])

        if result == (1,0):
            genos[idx]=random.choices([(0,1),(1,0)])

        else:
            genos[idx] = result

    for idx in g1a:
        genos[idx]=random.choices(
            [(0,0),(1,0),(1,1)], weights=error_matrix[['p10','p11','p12']].values[0])

    for idx in g1b:
        genos[idx]=random.choices(
            [(0,0),(0,1),(1,1)], weights=error_matrix[['p10','p11','p12']].values[0])

    for idx in g2:
        result=random.choices(
            [(0,0),(1,0),(1,1)], weights=error_matrix[['p20','p21','p22']].values[0])

        if result == (1,0):
            genos[idx]=random.choices([(0,1),(1,0)])

        else:
            genos[idx] = result

    return(np.array(sum([chromo for tup in genos for chromo in tup], ())))



def generate_samples_empirical(ts, filename, error_matrix):

    """

    Generate a samples file from a simulated ts based on an empirically estimated error matrix

    Reject any variants that result in a fixed column. 

    """

    error_matrix=pd.read_csv(error_matrix)

    record_rate = logging.getLogger().isEnabledFor(logging.INFO)

    n_variants = bits_flipped = 0

    assert ts.num_sites != 0

    sample_data = tsinfer.SampleData(
        path="../data/"+filename + "/"+filename+".samples",
        sequence_length=ts.sequence_length)

    for v in ts.variants():
        n_variants += 1

        #Record the allele frequency

        m = v.genotypes.shape[0]

        frequency = np.sum(v.genotypes) / m

        

        #Find closest row in error matrix file

        closest_freq = error_matrix.iloc[(error_matrix['freq']-frequency).abs().argsort()[:1]]

       

        #make new genotypes with error

        # Reject any columns that have no 1s or no zeros.

        # Unless the original also has them, as occasionally we have

        # some sims (e.g. under selection) where a variant is fixed

        while True:

            genotypes = make_errors_genotype_model(v.genotypes,closest_freq)

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

  

    logging.info(": actual error rate = {}".format(

            bits_flipped/(n_variants*ts.sample_size)) if record_rate else "")

    sample_data.finalise()

    return sample_data

