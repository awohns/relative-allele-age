import compare

import numpy as np
import itertools
import tsinfer

"""
Module for returning frequency values
For creation of site frequency spectrums
Only considers non-singleton sites
"""



def all_frequencies(sample_file,difference=False):
    
    all_freq=list()
    num_vars=len(sample_file.sites_genotypes[:])
    mutation_pairs=itertools.combinations(range(0,num_vars),2)
       
    singletons=compare.find_singletons_sample(sample_file)
    
    for var1,var2 in mutation_pairs:
        if (var1 not in singletons) and (var2 not in singletons):
            geno1=sample_file.sites_genotypes[var1]
            geno2=sample_file.sites_genotypes[var2]
            if difference == False:
                all_freq.append(len(geno1[geno1==1])/len(geno1))
                all_freq.append(len(geno2[geno2==1])/len(geno2))
            elif difference == True:
                all_freq.append(abs(len(geno1[geno1==1])/len(geno1) -len(geno2[geno2==1])/len(geno2)))

    return(all_freq)



def misordered_frequencies(sample_file,binary_matrix,second_binary=None,difference=False):
    freq_misordered=list()
    
    if second_binary is not None:
        binary_matrix_combo = np.add(binary_matrix,second_binary)
        binary_matrix_triu=np.triu(binary_matrix_combo)

        misordered_pairs=np.where(binary_matrix_triu==2)

    else:
        binary_matrix_triu=np.triu(binary_matrix)
        misordered_pairs=np.where(binary_matrix_triu==1)

    singletons=compare.find_singletons_sample(sample_file)

    for var1,var2 in zip(misordered_pairs[0],misordered_pairs[1]):
        if (var1 not in singletons) and (var2 not in singletons):
            geno1=sample_file.sites_genotypes[var1]
            geno2=sample_file.sites_genotypes[var2]
            freq1=len(geno1[geno1==1])/len(geno1)
            freq2=len(geno2[geno2==1])/len(geno2)
            if difference == False:
                freq_misordered.append(freq1)
                freq_misordered.append(freq2)
            elif difference == True:
                freq_misordered.append(abs(freq1-freq2))


    return(freq_misordered)





# def main():
# 	sample_file=tsinfer.load("../data/new_test_dir_again/new_test_dir_again.samples")
# 	freq_matrix=np.loadtxt("../data/new_test_dir_again/freq_matrix")
# 	geva_matrix=np.loadtxt("../data/new_test_dir_again/geva_matrix")
# 	print(misordered_frequencies(sample_file,freq_matrix))
# 	print(misordered_frequencies(sample_file,geva_matrix))


# if __name__ == "__main__":
#     main()



