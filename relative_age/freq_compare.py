"""
Frequency comparisons using samples object
"""

import msprime

#calculate frequency of variants
def freq_relative_time(msprime_ts,samples_file,delete_singletons):
    num_mutations = len(samples_file.sites_genotypes[:])
    pairwise_matrix_frequency = np.zeros((num_mutations,num_mutations))
    variants=list(msprime_ts.variants())
    
    for i,outer_geno in enumerate(samples_file.sites_genotypes[:]):
        outer_age = int(msprime_ts.node(variants[i].site.mutations[0].node).time)
        outer_freq = len(outer_geno[outer_geno ==1])/len(outer_geno)
        for j,inner_geno in enumerate(sample.sites_genotypes[:]):
            inner_age = int(msprime_ts.node(variants[j].site.mutations[0].node).time)
            inner_freq = len(inner_geno[inner_geno ==1])/len(inner_geno)
            
            #if the frequency is the same, can't compare ages and assign to 0
            if outer_age == inner_age or outer_freq == inner_freq:
                pairwise_matrix_frequency[i,j] = 0
            
            #if both the real age and the frequencies agree in relative magnitude, assign to 0
            elif ((outer_age < inner_age) == (outer_freq < inner_freq)):
                pairwise_matrix_frequency[i,j] = 0
            
            else:
                
                pairwise_matrix_frequency[i,j] = 1
    
    if delete_singletons:
        singletons = sample_find_singletons(samples_file)
        pairwise_matrix_frequency= np.delete(pairwise_matrix_frequency, singletons, 0)
        pairwise_matrix_frequency= np.delete(pairwise_matrix_frequency, singletons, 1)
        
    return(pairwise_matrix_frequency)


"""
Frequency comparisons old version (using tree sequence object)
"""

import msprime


#calculate frequency of variants
def old_freq_relative_time(msprime_ts,delete_singletons):
    num_mutations = len(list(msprime_ts.variants()))
    pairwise_matrix_frequency = np.zeros((num_mutations,num_mutations))
    for variant_outer in msprime_ts.variants():  
        outer_age = int(msprime_ts.node(variant_outer.site.mutations[0].node).time)
        outer_freq = len(variant_outer.genotypes[variant_outer.genotypes ==1])/len(variant_outer.genotypes)
        for variant_inner in msprime_ts.variants():
            inner_age = int(msprime_ts.node(variant_inner.site.mutations[0].node).time)
            inner_freq = len(variant_inner.genotypes[variant_inner.genotypes ==1])/len(variant_inner.genotypes)
            
            #if the frequency is the same, can't compare ages and assign to 0
            if outer_age == inner_age or outer_freq == inner_freq:
                pairwise_matrix_frequency[variant_outer.index,variant_inner.index] = 0
            
            #if both the real age and the frequencies agree in relative magnitude, assign to 0
            elif ((outer_age < inner_age) == (outer_freq < inner_freq)):
                pairwise_matrix_frequency[variant_outer.index,variant_inner.index] = 0
            
            
            else:

                pairwise_matrix_frequency[variant_outer.index,variant_inner.index] = 1
    if delete_singletons:
        singletons = ts_find_singletons(msprime_ts)
        pairwise_matrix_frequency= np.delete(pairwise_matrix_frequency, singletons, 0)
        pairwise_matrix_frequency= np.delete(pairwise_matrix_frequency, singletons, 1)
        
    return(pairwise_matrix_frequency)