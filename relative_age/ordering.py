import compare
import simulations
import distance
import frequency

import numpy as np
import csv
import os, errno


import warnings
warnings.filterwarnings("ignore")
from argparse import ArgumentParser




class RelativeTimeResults(object):
    """
    Matrices of relative time ordering results
    """
    def __init__(
            self,num_samples,Ne,length,mut_rate,rec_rate,msprime_ts,
            error_sample,geva_age_estimates,freq_matrix,freq_matrix_no_singletons,
            geva_matrix,geva_matrix_no_singletons,direct_matrix,direct_matrix_no_singletons,
            error_rate):
        self.samples = num_samples
        self.ne = Ne
        self.length = length
        self.mut_rate = mut_rate
        self.rec_rate = rec_rate

        self.msprime_ts=msprime_ts
        self.error_sample=error_sample
        self.geva_age_estimates=geva_age_estimates

        #matricies with singletons
        self.direct_matrix = direct_matrix
        self.freq_matrix = freq_matrix
        self.geva_matrix = geva_matrix

        #matricies without singletons
        self.direct_matrix_no_singletons = direct_matrix_no_singletons
        self.freq_matrix_no_singletons = freq_matrix_no_singletons
        self.geva_matrix_no_singletons = geva_matrix_no_singletons

        #if the matrices aren't the same size, something went wrong
        try:
            assert self.geva_matrix_no_singletons.size == self.freq_matrix_no_singletons.size == self.direct_matrix_no_singletons.size
        except AssertionError as e:
            e.args += ("ERROR: results matrices (no singletons) not same size")
            raise

        try:
            assert self.geva_matrix.size == self.freq_matrix.size == self.direct_matrix.size
        except AssertionError as e:
            e.args += ("ERROR: results matrices not same size")
            raise

        #take care of any nan values in the GEVA output
        where_are_NaNs = np.isnan(self.geva_matrix_no_singletons)
        self.geva_matrix_no_singletons[where_are_NaNs] = -99999

        #assign the size of the matrix as baseline
        self.baseline = self.freq_matrix_no_singletons.size
        self.baseline_comparable = len(self.direct_matrix_no_singletons[self.direct_matrix_no_singletons==1])                           
        
        direct_freq = np.add(self.freq_matrix_no_singletons,self.direct_matrix_no_singletons)                         
        direct_geva = np.add(self.geva_matrix_no_singletons,self.direct_matrix_no_singletons)
        
        #output distance separating directly comparable pairs of mutations
        #always outputs distances excluding singletons
        self.physical_distances_no_singletons = distance.get_pairwise_physical_distance(self.msprime_ts,self.error_sample)
        self.genetic_distances_no_singletons = distance.get_pairwise_genetic_distance(self.msprime_ts,self.error_sample)
        self.freq_direct_error_distance = distance.direct_comparison_pair_distance(msprime_ts,self.freq_matrix_no_singletons,
            self.direct_matrix_no_singletons,self.physical_distances_no_singletons,self.genetic_distances_no_singletons)
        self.geva_direct_error_distance = distance.direct_comparison_pair_distance(msprime_ts,self.geva_matrix_no_singletons,
            self.direct_matrix_no_singletons,self.physical_distances_no_singletons,self.genetic_distances_no_singletons)

        #output distance between all misordered pairs
        self.freq_all_error_distance = distance.all_misordered_pair_distance(msprime_ts,self.freq_matrix_no_singletons,
            self.physical_distances_no_singletons,self.genetic_distances_no_singletons)
        self.geva_all_error_distance = distance.all_misordered_pair_distance(msprime_ts,self.geva_matrix_no_singletons,
            self.physical_distances_no_singletons,self.genetic_distances_no_singletons)

        #get average distance between all pairs
        self.all_pair_average_distance = distance.all_pair_average_distance(msprime_ts,
           self.direct_matrix_no_singletons,self.physical_distances_no_singletons,self.genetic_distances_no_singletons)

        #get average distance between directly comparable pairs
        self.direct_comparison_distances = distance.direct_distances(msprime_ts,
           self.direct_matrix_no_singletons,self.physical_distances_no_singletons,self.genetic_distances_no_singletons)
        
        self.accuracy_by_distance = distance.accuracy_by_distance(self.msprime_ts,self.freq_matrix_no_singletons,self.geva_matrix_no_singletons,
            self.direct_matrix_no_singletons,self.physical_distances_no_singletons,self.genetic_distances_no_singletons)
       

        self.freq_metrics = {
            "correct" : len(self.freq_matrix_no_singletons[self.freq_matrix_no_singletons==0]),
            "incorrect" : len(self.freq_matrix_no_singletons[self.freq_matrix_no_singletons==1]),
            "comparable_incorrect" : len(direct_freq[direct_freq==2]),
            "average direct distance":np.nanmean(self.freq_direct_error_distance[:,1]),
            "average all distance":np.nanmean(self.freq_all_error_distance[:,1])
        }


        self.geva_metrics = {
            "correct" : len(self.geva_matrix_no_singletons[self.geva_matrix_no_singletons==0]),
            "incorrect" : len(self.geva_matrix_no_singletons[self.geva_matrix_no_singletons==1]),
            "comparable_incorrect" : len(direct_geva[direct_geva==2]),
            "average direct distance":np.nanmean(self.geva_direct_error_distance[:,1]),
            "average all distance":np.nanmean(self.geva_all_error_distance[:,1])
        }
        
        

 

    def write_summary_stats(self,output_filename):
        if self.baseline_comparable != 0:
            direct_freq_accuracy=self.freq_metrics["comparable_incorrect"]/self.baseline_comparable
            direct_geva_accuracy=self.geva_metrics["comparable_incorrect"]/self.baseline_comparable
        else:
            direct_freq_accuracy=0
            direct_geva_accuracy=0
        average_direct_physical_distance=np.nanmean(self.direct_comparison_distances[:,1])
        average_all_physical_distance=np.nanmean(self.all_pair_average_distance[:,1])
        with open("../data/"+output_filename+"/summary_stats", "a") as out:
            csv_out=csv.writer(out)
            row=[self.freq_metrics["correct"],self.freq_metrics["incorrect"],self.geva_metrics["correct"],
            self.geva_metrics["incorrect"],direct_freq_accuracy,direct_geva_accuracy,self.freq_metrics["average direct distance"],
                self.geva_metrics["average direct distance"],average_direct_physical_distance,self.freq_metrics["average all distance"],
                self.geva_metrics["average all distance"],average_all_physical_distance]
            csv_out.writerow(row)
        
    def update_accuracy_by_distance(self,output_filename):

        with open("../data/"+output_filename+"/freq_accuracy_by_distance", "a") as out:
            csv_out=csv.writer(out)
            csv_out.writerow(self.accuracy_by_distance[0,:])
        with open("../data/"+output_filename+"/geva_accuracy_by_distance", "a") as out:
            csv_out=csv.writer(out)
            csv_out.writerow(self.accuracy_by_distance[1,:])

    def manually_check_accuracy(self,output_filename): 
        np.savetxt("../data/"+output_filename+"/freq_matrix_no_singletons",
            self.freq_matrix_no_singletons,fmt='%i')
        np.savetxt("../data/"+output_filename+"/geva_matrix_no_singletons",
            self.geva_matrix_no_singletons,fmt='%i')
        np.savetxt("../data/"+output_filename+"/direct_matrix_no_singletons",
            self.direct_matrix_no_singletons,fmt='%i')
        np.savetxt("../data/"+output_filename+"/physical_distances_no_singletons",
            self.physical_distances_no_singletons)
        np.savetxt("../data/"+output_filename+"/freq error distances",
            self.freq_direct_error_distance[:,1])
        np.savetxt("../data/"+output_filename+"/geva error distances",
            self.geva_direct_error_distance[:,1])
        
        self.msprime_ts.dump("../data/"+output_filename+"/simulated_ts")
        self.geva_age_estimates.to_csv("../data/"+output_filename+"/geva_age_estimates")

    """
        function to ouput: 1. all the frequencies of simulated mutations 
        2. tuple of misordered freq pairs by frequency
        3. tuple of misordered freq pairs by GEVA
    """
    def get_frequencies(self,output_filename):
        all_frequencies=frequency.all_frequencies(self.error_sample)
        freq_misordered_pairs=frequency.misordered_frequencies(self.error_sample,self.freq_matrix)
        geva_misordered_pairs=frequency.misordered_frequencies(self.error_sample,self.geva_matrix)
        with open("../data/"+output_filename+"/all_frequencies", "a") as out:
            csv_out=csv.writer(out)
            csv_out.writerow(all_frequencies)

        with open("../data/"+output_filename+"/freq_misordered_frequencies", "a") as out:
            csv_out=csv.writer(out,escapechar=' ', quoting=csv.QUOTE_NONE)
            csv_out.writerow(freq_misordered_pairs)

        with open("../data/"+output_filename+"/geva_misordered_frequencies", "a") as out:
            csv_out=csv.writer(out)
            csv_out.writerow(geva_misordered_pairs)

    """
    this only uses the msprime tree, not the sample data simulated with error!
    """
    def tree_imbalance(self,output_filename):
        freq_colless_correlation=imbalance.colless_correlate_error(self.msprime_ts,
            self.freq_matrix,self.direct_matrix)
        freq_sackin_correlation=imbalance.sackin_correlate_error(self.msprime_ts,
            self.freq_matrix,self.direct_matrix)

    #CURRENTLY WORKING ON THIS
    # def overlap_geva_freq_misordering(self):
    #     np.equal(self.freq_matrix_no_singletons,self.geva_matrix_no_singletons)



# class SummaryStats(object):
#     """
#     Summary stats of relative time ordering results
#     """ 
  
