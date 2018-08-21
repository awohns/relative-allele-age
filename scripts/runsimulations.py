"""
Top-level module for running multiple simulations to 
test relative time ordering. The implementation of
age estimation, metrics on resulting trees, error
generation etc. are in other modules. 
"""
import sys
sys.path.insert(0, '../relative_age')
import simulations
import compare
import ordering

import msprime
import numpy as np
np.set_printoptions(threshold=np.nan)
import itertools
from itertools import combinations
from collections import defaultdict
import tsinfer
import pandas as pd
import time
from tqdm import tqdm
from scipy import stats
import csv
import warnings
warnings.filterwarnings("ignore")
from argparse import ArgumentParser
import os

 
#multiple replicates at given parameters
def multiple_replicates(replicates, samples, Ne, length, mut_rate, rec_rate,error_rate,delete_singletons,output,manual_check):
    replicates = int(replicates)
    summary_stats = np.zeros((replicates,2))
    geva_corrected = list()

    for replicate_num in tqdm(range(0,replicates)):
        
        msprime_ts = simulations.msprime_simulation(output, samples, Ne, length, mut_rate, rec_rate)
        error_sample = simulations.generate_samples(msprime_ts,output,float(error_rate))

        freq_matrix=compare.freq_relative_time(msprime_ts,error_sample,delete_singletons)
        geva_matrix,geva_age_estimates=compare.geva_all_time_orderings(output,msprime_ts,error_sample,Ne,length,mut_rate,delete_singletons)
        direct_matrix= compare.directly_comparable(msprime_ts,error_sample,delete_singletons)


        results = ordering.RelativeTimeResults(
            samples,Ne,length,mut_rate,rec_rate,msprime_ts,error_sample,geva_age_estimates,
            freq_matrix,geva_matrix,direct_matrix,error_rate,delete_singletons)

        

        results.write_summary_stats(output)
        results.update_accuracy_by_distance(output)
        
        if manual_check == True:
            results.manually_check_accuracy(output)

        results.get_frequencies(output)




def main():
    parser = ArgumentParser()
    
    parser.add_argument('replicates')
    parser.add_argument("samples")
    parser.add_argument("Ne")
    parser.add_argument("length")
    parser.add_argument("mut_rate")
    parser.add_argument("rec_rate")
    parser.add_argument("error_rate")
    parser.add_argument("delete_singletons")
    parser.add_argument("output")
    parser.add_argument('-manual_check', action='store_true')

    args = parser.parse_args()

    try:
        os.makedirs("../data/"+args.output,exist_ok=True)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    multiple_replicate_results = multiple_replicates(int(args.replicates),int(args.samples),
        int(args.Ne),int(args.length),float(args.mut_rate),float(args.rec_rate),args.error_rate,args.delete_singletons,
        str(args.output),args.manual_check)
    
   

    # with open(args.output, "a") as out:
    #     csv_out=csv.writer(out)
    #     csv_out.writerow(multiple_replicate_results)
        # for tup in multiple_replicate_results:
        #     csv_out.writerow(tup)


if __name__ == "__main__":
    main()
