import compare

import numpy as np
import tsinfer

"""
Module for returning frequency values
For creation of site frequency spectrums
Only considers non-singleton sites
"""



def all_frequencies(sample_file):

    all_freq=list()
    num_vars=len(sample_file.sites_genotypes[:])

    for variant in sample_file.sites_genotypes[:]:
        for i in range(0,num_vars):
	        all_freq.append(len(variant[variant==1])/len(variant))
    
    #Delete singletons
    singletons=compare.find_singletons_sample(sample_file)

    all_freq_no_singletons=[x for i,x in enumerate(all_freq) if i not in singletons]

    return(all_freq_no_singletons)

def misordered_frequencies(sample_file,binary_matrix):
	freq_misordered=list()
	
	#only consider each pair once
	binary_matrix_triu=np.triu(binary_matrix)

	misordered_pairs=np.where(binary_matrix_triu==1)

	for var1,var2 in zip(misordered_pairs[0],misordered_pairs[1]):

		geno1=sample_file.sites_genotypes[var1]
		geno2=sample_file.sites_genotypes[var2]
		freq1=len(geno1[geno1==1])/len(geno1)
		freq2=len(geno2[geno2==1])/len(geno2)
		freq_misordered.append(freq1)
		freq_misordered.append(freq2)

	singletons=compare.find_singletons_sample(sample_file)
	freq_misordered_no_singletons=[x for i,x in enumerate(freq_misordered) if i not in singletons]
	return(freq_misordered)



# def main():
# 	sample_file=tsinfer.load("../data/new_test_dir_again/new_test_dir_again.samples")
# 	freq_matrix=np.loadtxt("../data/new_test_dir_again/freq_matrix")
# 	geva_matrix=np.loadtxt("../data/new_test_dir_again/geva_matrix")
# 	print(misordered_frequencies(sample_file,freq_matrix))
# 	print(misordered_frequencies(sample_file,geva_matrix))


# if __name__ == "__main__":
#     main()



