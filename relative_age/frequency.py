import compare

import numpy as np
import tsinfer

"""
Module for returning frequency values
For creation of site frequency spectrums
Only considers non-singleton sites
"""



def all_frequencies(msprime_ts, sample_file):

    all_freq=list()
    for variant in sample_file.sites_genotypes[:]:
        all_freq.append(len(variant[variant==1])/len(variant))
    
    #Delete singletons
    singletons=compare.find_singletons_sample(sample_file)

    all_freq_no_singletons=[x for i,x in enumerate(all_freq) if i not in singletons]

    return(all_freq_no_singletons)

# def freq_misordered_frequencies(msprime_ts,sample_file):
# 	freq_misordered=list()
# 	for variant in sample_


# def main():
# 	tsinfer.load("../data/")
# 	print(singletons=compare.find_singletons_sample(sample_file))


# if __name__ == "__main__":
#     main()