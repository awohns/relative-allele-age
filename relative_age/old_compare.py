
def find_singletons(msprime_ts):
    singleton_list = list()
    for variant in msprime_ts.variants():
        if len(variant.genotypes[variant.genotypes ==1]) == 1:
            singleton_list.append(variant.index)
        if len(variant.genotypes[variant.genotypes == 1]) == (len(variant.genotypes)-1):
            singleton_list.append(variant.index)
        if 0 not in singleton_list:
            singleton_list.append(0)
        last_value = len(list(msprime_ts.variants())) - 1
        if last_value not in singleton_list:
            singleton_list.append(last_value)
    return(singleton_list)


#removes singletons
def old_geva_relative_time(msprime_ts, vcf_name, Ne, mutation_rate,directly_comparable_matrix,delete_singletons):
    num_mutations = len(list(msprime_ts.variants())) 
    pairwise_matrix_geva = np.zeros((num_mutations,num_mutations))

    #age estimation
    age_estimates = geva_age_estimate(vcf_name,Ne,mutation_rate)

    #Loop through all pairs of variants
    #for each pair, determine whether one should be older than other by msprime
    #then check whether GEVA: a. confirms, b. refutes, c. can't say
    for variant_outer in msprime_ts.variants():  
        outer_age = int(msprime_ts.node(variant_outer.site.mutations[0].node).time)

        outer_age_estimate = find_age(variant_outer.index,age_estimates)

        for variant_inner in msprime_ts.variants():


            inner_age = int(msprime_ts.node(variant_inner.site.mutations[0].node).time)
            inner_age_estimate = find_age(variant_inner.index,age_estimates)
            if outer_age_estimate is not None and inner_age_estimate is not None:
                if outer_age == inner_age:
                    pairwise_matrix_geva[variant_outer.index,variant_inner.index] = 0
                else:
                    if ((outer_age < inner_age) == (outer_age_estimate < inner_age_estimate)):
                        pairwise_matrix_geva[variant_outer.index,variant_inner.index] = 0
                    else:
                        pairwise_matrix_geva[variant_outer.index,variant_inner.index] = 1
            else:
                pairwise_matrix_geva[variant_outer.index,variant_inner.index] = np.nan
    
    if delete_singletons:
        singletons = find_singletons(msprime_ts)
        pairwise_matrix_geva= np.delete(pairwise_matrix_geva, singletons, 0)
        pairwise_matrix_geva= np.delete(pairwise_matrix_geva, singletons, 1)
        
        
    return(pairwise_matrix_geva)


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