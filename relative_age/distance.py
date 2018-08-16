import msprime
import numpy as np

import compare


def get_pairwise_genetic_distance(msprime_ts,samples_data,delete_singletons):
    #create matrix of pairwise genetic distances between mutations
    num_mutations = len(list(msprime_ts.variants())) 
    genetic_distances = np.zeros((num_mutations,num_mutations))

    for outer_tree in msprime_ts.trees():
        for outer_mutation in outer_tree.mutations():
            for inner_tree in msprime_ts.trees():
                for inner_mutation in inner_tree.mutations():
                    genetic_distances[outer_mutation.index,inner_mutation.index] = abs(outer_tree.index-inner_tree.index)
    
    if delete_singletons == "True":
        singletons = compare.find_singletons_sample(samples_data)
        genetic_distances= np.delete(genetic_distances, singletons, 0)
        genetic_distances= np.delete(genetic_distances, singletons, 1)
       
    return(genetic_distances)


def get_pairwise_physical_distance(msprime_ts,samples_data,delete_singletons):
    #create matrix of pairwise physical distances between mutations
    num_mutations = len(list(msprime_ts.variants())) 
    physical_distances = np.zeros((num_mutations,num_mutations))

    for outer_mutation in msprime_ts.mutations():
        for inner_mutation in msprime_ts.mutations():
            physical_distances[outer_mutation.index,inner_mutation.index] = abs(outer_mutation.position-inner_mutation.position)
    
    if delete_singletons == "True":
        singletons = compare.find_singletons_sample(samples_data)
        physical_distances= np.delete(physical_distances, singletons, 0)
        physical_distances= np.delete(physical_distances, singletons, 1)
       
    return(physical_distances)

def mutation_pair(mutation_1, mutation_2, msprime_ts,physical_distance,genetic_distance):
    #1. how far apart are the two mutations
    genetic_distance = genetic_distance[mutation_1,mutation_2]
    physical_distance = physical_distance[mutation_1,mutation_2]
    return(genetic_distance,physical_distance)
    #how many significant trees are they apart
    
    #how many nodes are they apart in the second

def diagnose_error(msprime_ts, pairwise_matrix_binary,direct_comparisons,physical_distance,genetic_distance):
    #look up which pairs of mutations are in the wrong order in pairwise array (either from freq or geva)
    direct_errors = np.add(pairwise_matrix_binary,direct_comparisons)
    wrong_order = np.argwhere(direct_errors == 2)
    results = np.zeros((len(wrong_order), 2))
    
    #wrong_order first array is first mutation's index, second array has second mutations's indices
    for idx,(mut1_index,mut2_index) in enumerate(wrong_order):
        cur_result = mutation_pair(mut1_index,mut2_index,msprime_ts,physical_distance,genetic_distance)
        results[idx,0] = cur_result[0]
        results[idx,1] = cur_result[1]
        #results[idx,2] = cur_result[2]
        
    #return: 1. array of genetic distances of pairwise error; 2. array of physical distances of pairwise error
    return(results)

def direct_distances(msprime_ts,direct_comparisons,physical_distance,genetic_distance):
    #look up which pairs of mutations are in the wrong order in pairwise array (either from freq or geva)
    all_direct = np.argwhere(direct_comparisons == 1)
    results = np.zeros((len(all_direct), 2))
    
    #all_direct first array is first mutation's index, second array has second mutations's indices
    for idx,(mut1_index,mut2_index) in enumerate(all_direct):
        cur_result = mutation_pair(mut1_index,mut2_index,msprime_ts,physical_distance,genetic_distance)
        results[idx,0] = cur_result[0]
        results[idx,1] = cur_result[1]
        #results[idx,2] = cur_result[2]
        
    #return: 1. array of genetic distances of pairwise error; 2. array of physical distances of pairwise error
    return(results)

#function to get frequency/geva accuracy by distance
def accuracy_by_distance(msprime_ts,freq_matrix,geva_matrix,direct_matrix,physical_distance,genetic_distance):
    #look up which pairs of mutations are in the wrong order in pairwise array (either from freq or geva)
    direct_errors = np.add(pairwise_matrix_binary,direct_comparisons)
    direct_baseline = np.argwhere(direct_matrix == 1)
    wrong_order = np.argwhere(direct_errors == 2)
    results = np.zeros((len(direct_baseline), 4))

    #Iterate through each mutation pair that is directly comparable, recording distance
    #Check if correct by frequency/GEVA. Make tuple of (idx,distance,freq_correct,geva_correct)
    for idx,(mut1_index,mut2_index) in enumerate(direct_baseline):
        cur_result = mutation_pair(mut1_index,mut2_index,msprime_ts,physical_distance,genetic_distance)


    #Using the resulting numpy array of tuples, make regularly spaced bins of 1000 bp and record accuracy
    #for each bin (number of 1's divded by 0's + 1's)





