"""
Script to compare frequency with GEVA
"""

import random
import itertools
import logging
import argparse
import multiprocessing
import subprocess
from argparse import ArgumentParser
import os
import random
import glob
import math

import pandas as pd
import numpy as np
import tqdm

import msprime
import tsinfer


error_matrix=pd.read_csv("EmpiricalErrorPlatinum1000G.csv")

def make_errors_genotype_model(g, error_probs):
    """
    Given an empirically estimated error probability matrix, resample for a particular
    variant. Determine variant frequency and true genotype (g0, g1, or g2),
    then return observed genotype based on row in error_probs with nearest
    frequency. Treat each pair of alleles as a diploid individual.
    """
    w = np.copy(g)

    # Make diploid (iterate each pair of alleles)
    genos = [(w[i], w[i+1]) for i in range(0, w.shape[0], 2)]

    # Record the true genotypes
    g0 = [i for i, x in enumerate(genos) if x == (0, 0)]
    g1a = [i for i, x in enumerate(genos) if x == (1, 0)]
    g1b = [i for i, x in enumerate(genos) if x == (0, 1)]
    g2 = [i for i, x in enumerate(genos) if x == (1, 1)]

    for idx in g0:
        result = [(0, 0), (1, 0), (1, 1)][
            np.random.choice(3, p=error_probs[['p00', 'p01', 'p02']].values[0])]
        if result == (1, 0):
            genos[idx] = [(0, 1), (1, 0)][np.random.choice(2)]
        else:
            genos[idx] = result
    for idx in g1a:
        genos[idx] = [(0, 0), (1, 0), (1, 1)][
            np.random.choice(3, p=error_probs[['p10', 'p11', 'p12']].values[0])]
    for idx in g1b:
        genos[idx] = [(0, 0), (0, 1), (1, 1)][
            np.random.choice(3, p=error_probs[['p10', 'p11', 'p12']].values[0])]
    for idx in g2:
        result = [(0, 0), (1, 0), (1, 1)][
            np.random.choice(3, p=error_probs[['p20', 'p21', 'p22']].values[0])]
        if result == (1, 0):
            genos[idx] = [(0, 1), (1, 0)][np.random.choice(2)]
        else:
            genos[idx] = result

    return np.array(sum(genos, ()))    
    

def generate_samples_empirical(ts):

    """
    Generate a samples file from a simulated ts based on an empirically estimated error matrix
    Reject any variants that result in a fixed column. 
    """
    assert ts.num_sites != 0
    sample_data = tsinfer.SampleData(sequence_length=ts.sequence_length)

    for v in ts.variants():
        #Record the allele frequency
        m = v.genotypes.shape[0]
        frequency = np.sum(v.genotypes) / m

        #Find closest row in error matrix file
        closest_freq = error_matrix.iloc[(error_matrix['freq']-frequency).abs().argsort()[:1]]

        #make new genotypes with error
        # Reject any columns that have no 1s or no zeros.
        # Unless the original also has them, as occasionally we have
        # some sims (e.g. under selection) where a variant is fixed
        genotypes = make_errors_genotype_model(v.genotypes,closest_freq)
        
        sample_data.add_site(
            position=v.site.position, alleles=v.alleles,
            genotypes=genotypes)

    sample_data.finalise()
    return sample_data

def generate_samples(ts):
    """
    Generate a samples file from a simulated ts
    Samples may have bits flipped with a specified probability.
    (reject any variants that result in a fixed column)
    """

    assert ts.num_sites != 0
    
    sample_data = tsinfer.SampleData(sequence_length=ts.sequence_length)
    for v in ts.variants():
        sample_data.add_site(
            position=v.site.position, alleles=v.alleles,
            genotypes=v.genotypes)

    sample_data.finalise()
    return sample_data

def geva_age_estimate(file_name,ne,rec_rate,mut_rate):
    subprocess.call(["/home/wilderwohns/relative_allele_age/bin/geva/geva_v1beta", '--vcf', file_name+".vcf","--rec", rec_rate, "--out",file_name],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(file_name+".positions.txt","wb") as out:
        subprocess.call(["awk", "NR>3 {print last} {last = $3}", file_name+".marker.txt"], stdout=out)
    with open(file_name+".positions.nodup.txt","wb") as out:
        subprocess.call(["awk", '!NF || !seen[$0]++', file_name+".positions.txt"], stdout=out)
    try:
        subprocess.call(["/home/wilderwohns/relative_allele_age/bin/geva/geva_v1beta", "-i", file_name+".bin", "--positions", file_name+".positions.nodup.txt","--hmm","/home/wilderwohns/relative_allele_age/bin/geva/hmm/hmm_initial_probs.txt ","/home/wilderwohns/relative_allele_age/bin/geva/hmm/hmm_emission_probs.txt","--Ne", ne, "--mut", mut_rate,"-o",file_name+"_estimation"],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as grepexc:
        print(grepexc.output)
    age_estimates = pd.read_csv(file_name+"_estimation.sites.txt", sep = " ")
    return(age_estimates)

def geva_with_error(name,sample_data,Ne,length,mut_rate):
    
    num_individuals = len(sample_data.individuals_metadata[:])
    ind_list = list()
    pos_geno_dict = {"POS":list()}

    for i in range(int(num_individuals/2)):
        pos_geno_dict["msp_"+str(i)] = list()
        ind_list.append("msp_"+str(i))

    #add all the sample positions and genotypes
    for i in sample_data.genotypes():
        pos_geno_dict["POS"].append(int(round(sample_data.sites_position[i[0]])))
        for j in range(0,len(i[1]),2):

            pos_geno_dict["msp_"+str(int(j/2))].append(str(i[1][j]) + "|" + str(i[1][j+1]))

    df = pd.DataFrame(pos_geno_dict)

    df["#CHROM"] = 1
    df["REF"] = "A"
    df["ALT"] = "T"
    df['ID'] = "."
    df['QUAL'] = "."
    df['FILTER'] = "PASS"
    df['INFO'] = "."
    df['FORMAT'] = "GT"

    cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT']+ind_list
    df = df[cols] 

    header = """##fileformat=VCFv4.2
##source=msprime 0.6.0
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length=""" + str(int(length)) +  """>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
    
    output_VCF = "tmp/"+name+".vcf"
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)

    df.to_csv(output_VCF, sep="\t", mode='a',index=False)
    
    return(geva_age_estimate("tmp/"+name,"10000","1e-8","1e-8"))

def recursive_search(combo,cur_node,origin_node,target_node,mutations,parent_table,child_table):
    for node in parent_table[child_table == cur_node]:
        if node == target_node:
            return(mutations[combo[0]].id,mutations[combo[1]].id)
        else:
            return(recursive_search(combo,node,origin_node,target_node,mutations,parent_table,child_table))


def run_comparisons(params):
    sample_size, Ne, length, rec_rate, mut_rate, seed = params
    n_bins = 20
    bin_size = length//n_bins
    assert bin_size == length/n_bins
    total = np.zeros(n_bins, dtype=np.int)
    total_error = np.zeros(n_bins, dtype=np.int)
    total_geva = np.zeros(n_bins, dtype=np.int)
    total_geva_error = np.zeros(n_bins, dtype=np.int)
    agree = np.zeros(n_bins, dtype=np.int)
    agree_error = np.zeros(n_bins, dtype=np.int)
    agree_geva = np.zeros(n_bins, dtype=np.int)
    agree_geva_error = np.zeros(n_bins, dtype=np.int)

    simulation = msprime.simulate(sample_size=sample_size, Ne=Ne, length=length, recombination_rate=rec_rate,mutation_rate=mut_rate, random_seed=seed)
        
    
    sample_data = generate_samples(simulation)
    sample_data_error = generate_samples_empirical(simulation)
    age_estimates = geva_with_error(str(seed)+"no_error",sample_data,10000,1e-8,1e-8)
    age_estimates_error = geva_with_error(str(seed)+"error",sample_data_error,10000,1e-8,1e-8)
    num_mutations = len(sample_data.sites_genotypes[:])
    mutations = list(simulation.mutations())
    variant_ages = [(mutation.id,simulation.node(mutation.node).time) for mutation in mutations]
    variant_frequencies = [len(geno[geno==1])/len(geno) for geno in sample_data.sites_genotypes[:]]
    variant_frequencies_error = [len(geno[geno==1])/len(geno) for geno in sample_data_error.sites_genotypes[:]]

    mutation_positions=[mutation.position for mutation in mutations]

    ancestrally_linked_pairs=list()
    
    tables=simulation.dump_tables()
    for combo in list(itertools.combinations(range(0,num_mutations),2)):
        result1=recursive_search(combo,mutations[combo[0]].node,mutations[combo[0]].node,mutations[combo[1]].node,mutations,tables.edges.parent,tables.edges.child)
        result2=recursive_search(combo,mutations[combo[1]].node,mutations[combo[1]].node,mutations[combo[0]].node,mutations,tables.edges.parent,tables.edges.child)
        if result1 is not None : ancestrally_linked_pairs.append(result1)
        if result2 is not None : ancestrally_linked_pairs.append(result2)
       


    for combo in ancestrally_linked_pairs:
        age1 = variant_ages[combo[0]]
        age2 = variant_ages[combo[1]]

        freq1 = variant_frequencies[combo[0]]
        freq2 = variant_frequencies[combo[1]]

        freq1_error = variant_frequencies_error[combo[0]]
        freq2_error = variant_frequencies_error[combo[1]]

        distance = int(abs(mutation_positions[combo[0]] - mutation_positions[combo[1]]))
        
        if (age1[1] != age2[1]) and (combo[0] != 0) and (combo[1] != 0) and (combo[0] != num_mutations - 1) and (combo[1] != num_mutations - 1):
            if (freq1 > 1/sample_size) and (freq2 > 1/sample_size) and (freq1 < (sample_size-1)/sample_size) and (freq2 < (sample_size-1)/sample_size):
        
                    
                total[distance//bin_size] += 1

                if (freq1 != freq2):
                    if ((age1[1] < age2[1]) == (freq1 < freq2)):
                        agree[distance//bin_size] += 1

                else:
                    if (random.choice([True, False])):
                        agree[distance//bin_size] += 1

                if (combo[0] in age_estimates.MarkerID.values) and (combo[1] in age_estimates.MarkerID.values):

                    total_geva[distance//5000] += 1
                    
                    geva_estimate_1 = age_estimates[(age_estimates.MarkerID == combo[0]) & (age_estimates.Filtered == 1) & (age_estimates.Clock == "J")].PostMode.values[0]
                    geva_estimate_2 = age_estimates[(age_estimates.MarkerID == combo[1]) & (age_estimates.Filtered == 1) & (age_estimates.Clock == "J")].PostMode.values[0]

                    if (geva_estimate_1 != geva_estimate_2):
                        if ((age1[1] < age2[1]) == (geva_estimate_1 < geva_estimate_2)):
                            agree_geva[distance//5000] += 1
                    else:
                        if (random.choice([True, False])):
                            agree_geva[distance//5000] += 1

            if (freq1_error > 1/sample_size) and (freq2_error > 1/sample_size) and (freq1_error < (sample_size-1)/sample_size) and (freq2_error < (sample_size-1)/sample_size):

                total_error[distance//5000] += 1
                if (freq1_error != freq2_error):
                    if ((age1[1] < age2[1]) == (freq1_error < freq2_error)):
                        agree_error[distance//bin_size] += 1

                else:
                    if (random.choice([True, False])):
                        agree_error[distance//bin_size] += 1

                if (combo[0] in age_estimates_error.MarkerID.values) and (combo[1] in age_estimates_error.MarkerID.values): 

                    geva_estimate_error_1 = age_estimates_error[(age_estimates_error.MarkerID == combo[0]) & (age_estimates_error.Filtered == 1) & (age_estimates_error.Clock == "J")].PostMode.values[0]
                    geva_estimate_error_2 = age_estimates_error[(age_estimates_error.MarkerID == combo[1]) & (age_estimates_error.Filtered == 1) & (age_estimates_error.Clock == "J")].PostMode.values[0]

                    total_geva_error[distance//5000] += 1

                    if (geva_estimate_error_1 != geva_estimate_error_2):
                        if ((age1[1] < age2[1]) == (geva_estimate_error_1 < geva_estimate_error_2)):
                            agree_geva_error[distance//5000] += 1
                    else:
                        if (random.choice([True, False])):
                            agree_geva_error[distance//5000] += 1

    for filename in glob.glob("tmp/"+str(seed)+"no_error*"):
        os.remove(filename) 
    for filename in glob.glob("tmp/"+str(seed)+"error*"):
        os.remove(filename) 
    return total, total_error, total_geva, total_geva_error, agree, agree_error, agree_geva, agree_geva_error


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
         '--replicates', '-r', type=int, default=4000, help="number of replicates")
    parser.add_argument(
         '--seed', '-s', type=int, default=123, help="use a non-default RNG seed")
    parser.add_argument(
        "--processes", '-p', type=int, default=1,
        help="number of worker processes, e.g. 40")
    parser.add_argument(
         '--progress',  "-P", action='store_true',
         help="Show a progress bar.", )


    args = parser.parse_args()

    np.random.seed(args.seed)
    random_seeds = np.random.randint(1 , 2**32, size = args.replicates)
    run_params = zip(
        itertools.repeat(50),     # sample_size
        itertools.repeat(10000),   # Ne 
        itertools.repeat(100000), # length
        itertools.repeat(1e-8), # recomb_rate
        itertools.repeat(1e-8), # mut_rate
        random_seeds # seed
        )

    # df = pd.DataFrame.from_dict({"Total":np.zeros(20), "TotalError":np.zeros(20), "TotalGeva":np.zeros(20), "TotalGevaError":np.zeros(20), "Agree":np.zeros(20), "AgreeError":np.zeros(20), "AgreeGeva":np.zeros(20), "AgreeGevaError":np.zeros(20)})
     
    df = pd.read_csv("data/ancestrally_linked_GEVA_frequency_distance_accuracy.csv")
    if args.processes > 1:
        logging.info("Setting up using multiprocessing ({} processes)".format(args.processes))
        with multiprocessing.Pool(processes=args.processes, maxtasksperchild=2) as pool:
            for total, total_error, total_geva, total_geva_error, agree, agree_error, agree_geva, agree_geva_error in tqdm.tqdm(
                pool.imap_unordered(run_comparisons, run_params), total=args.replicates, disable=not args.progress):
                df.Total += total
                df.TotalError += total_error
                df.TotalGeva += total_geva
                df.TotalGevaError += total_geva_error
                df.Agree += agree
                df.AgreeError += agree_error
                df.AgreeGeva += agree_geva
                df.AgreeGevaError += agree_geva_error
                df.to_csv("data/ancestrally_linked_GEVA_frequency_distance_accuracy.csv")
                
    else:
        # When we have only one process it's easier to keep everything in the same
        # process for debugging.
        logging.info("Setting up using a single process")
        for total, total_error, total_geva, total_geva_error, agree, agree_error, agree_geva, agree_geva_error in tqdm.tqdm(
            map(run_comparisons, run_params), total=args.replicates, disable=not args.progress):
                df.Total += total
                df.TotalError += total_error
                df.TotalGeva += total_geva
                df.TotalGevaError += total_geva_error
                df.Agree += agree
                df.AgreeError += agree_error
                df.AgreeGeva += agree_geva
                df.AgreeGevaError += agree_geva_error
                df.to_csv("data/ancestrally_linked_GEVA_frequency_distance_accuracy.csv")

    df.to_csv("data/ancestrally_linked_GEVA_frequency_distance_accuracy.csv")


