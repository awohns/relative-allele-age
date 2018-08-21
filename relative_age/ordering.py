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
            self, num_samples, Ne, length, mut_rate, rec_rate,msprime_ts,error_sample,geva_age_estimates,
             freq_matrix, geva_matrix, direct_matrix, error_rate,delete_singletons):
        self.samples = num_samples
        self.ne = Ne
        self.length = length
        self.mut_rate = mut_rate
        self.rec_rate = rec_rate

        self.msprime_ts=msprime_ts
        self.error_sample=error_sample
        self.geva_age_estimates=geva_age_estimates

        self.direct_matrix = direct_matrix
        self.freq_matrix = freq_matrix
        self.geva_matrix = geva_matrix
        self.delete_singletons = delete_singletons

        #if the matrices aren't the same size, something went wrong
        try:
            assert self.geva_matrix.size == self.freq_matrix.size == self.direct_matrix.size
        except AssertionError as e:
            e.args += ("ERROR: results matrices not same size")
            raise

        #take care of any nan values in the GEVA output
        where_are_NaNs = np.isnan(self.geva_matrix)
        self.geva_matrix[where_are_NaNs] = -99999

        #assign the size of the matrix as baseline
        self.baseline = self.freq_matrix.size
        self.baseline_comparable = len(self.direct_matrix[self.direct_matrix==1])                           
        
        direct_freq = np.add(self.freq_matrix,self.direct_matrix)                         
        direct_geva = np.add(self.geva_matrix,self.direct_matrix)
        
        #Now diagnose error where relative ordering is wrong
        self.physical_distances = distance.get_pairwise_physical_distance(self.msprime_ts,self.error_sample,delete_singletons)
        self.genetic_distances = distance.get_pairwise_genetic_distance(self.msprime_ts,self.error_sample,delete_singletons)
        self.freq_errors_diagnosis = distance.diagnose_error(msprime_ts,self.freq_matrix,
            self.direct_matrix,self.physical_distances,self.genetic_distances)
        self.geva_errors_diagnosis = distance.diagnose_error(msprime_ts,self.geva_matrix,
            self.direct_matrix,self.physical_distances,self.genetic_distances)

        #get average distance between directly comparable pairs
        self.direct_comparison_distances = distance.direct_distances(msprime_ts,
           self.direct_matrix,self.physical_distances,self.genetic_distances)
        
        self.accuracy_by_distance = distance.accuracy_by_distance(self.msprime_ts,self.freq_matrix,self.geva_matrix,
            self.direct_matrix,self.physical_distances,self.genetic_distances)
       

        self.freq_metrics = {
            "correct" : len(self.freq_matrix[self.freq_matrix==0]),
            "incorrect" : len(self.freq_matrix[self.freq_matrix==1]),
            "comparable_incorrect" : len(direct_freq[direct_freq==2]),
            "average distance":np.nanmean(self.freq_errors_diagnosis[:,1])
        }


        self.geva_metrics = {
            "correct" : len(self.geva_matrix[self.geva_matrix==0]),
            "incorrect" : len(self.geva_matrix[self.geva_matrix==1]),
            "comparable_incorrect" : len(direct_geva[direct_geva==2]),
            "average distance":np.nanmean(self.geva_errors_diagnosis[:,1])
        }
        
        

 

    def write_summary_stats(self,output_filename):
        if self.baseline_comparable != 0:
            direct_freq_accuracy=self.freq_metrics["comparable_incorrect"]/self.baseline_comparable
            direct_geva_accuracy=self.geva_metrics["comparable_incorrect"]/self.baseline_comparable
        else:
            direct_freq_accuracy=0
            direct_geva_accuracy=0
        average_direct_physical_distance=np.nanmean(self.direct_comparison_distances[:,1])
        with open("../data/"+output_filename+"/summary_stats", "a") as out:
            csv_out=csv.writer(out)
            row=[self.freq_metrics["correct"],self.freq_metrics["incorrect"],self.geva_metrics["correct"],
            self.geva_metrics["incorrect"],direct_freq_accuracy,direct_geva_accuracy,self.freq_metrics["average distance"],
                self.geva_metrics["average distance"],average_direct_physical_distance]
            csv_out.writerow(row)
        
    def update_accuracy_by_distance(self,output_filename):

        with open("../data/"+output_filename+"/freq_accuracy_by_distance", "a") as out:
            csv_out=csv.writer(out)
            csv_out.writerow(self.accuracy_by_distance[0,:])
        with open("../data/"+output_filename+"/geva_accuracy_by_distance", "a") as out:
            csv_out=csv.writer(out)
            csv_out.writerow(self.accuracy_by_distance[1,:])

    def manually_check_accuracy(self,output_filename): 
        np.savetxt("../data/"+output_filename+"/freq_matrix",self.freq_matrix,fmt='%i')
        np.savetxt("../data/"+output_filename+"/geva_matrix",self.geva_matrix,fmt='%i')
        np.savetxt("../data/"+output_filename+"/direct_matrix",self.direct_matrix,fmt='%i')
        np.savetxt("../data/"+output_filename+"/physical_distances",self.physical_distances)
        np.savetxt("../data/"+output_filename+"/freq error distances",self.freq_errors_diagnosis[:,1])
        np.savetxt("../data/"+output_filename+"/geva error distances",self.geva_errors_diagnosis[:,1])
     
        self.msprime_ts.dump("../data/"+output_filename+"/simulated_ts")
        self.geva_age_estimates.to_csv("../data/"+output_filename+"/geva_age_estimates")

    def get_frequencies(self,output_filename):
        """
        function to ouput: 1. all the frequencies of simulated mutations 
        2. frequencies of misordered younger mutations
        3. frequencies of misordered older mutations
        """

        all_frequencies=frequency.all_frequencies(self.msprime_ts,self.error_sample)
        # freq_misordered_pairs=frequency.misordered_frequencies(self.msprime_ts,self.error_sample)
        # geva_misordered_pairs=frequency.misordered_frequencies(self.msprime_ts,self.error_sample)
        with open("../data/"+output_filename+"/all_frequencies", "a") as out:
            csv_out=csv.writer(out)
            csv_out.writerow(all_frequencies)




# class SummaryStats(object):
#     """
#     Summary stats of relative time ordering results
#     """ 
  
