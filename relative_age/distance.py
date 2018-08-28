import msprime
import numpy as np

import compare


def get_pairwise_genetic_distance(msprime_ts,sample_data):
    #create matrix of pairwise genetic distances between mutations
    num_mutations = len(list(msprime_ts.variants())) 
    genetic_distances = np.zeros((num_mutations,num_mutations))

    for outer_tree in msprime_ts.trees():
        for outer_mutation in outer_tree.mutations():
            for inner_tree in msprime_ts.trees():
                for inner_mutation in inner_tree.mutations():
                    genetic_distances[outer_mutation.index,inner_mutation.index] = abs(outer_tree.index-inner_tree.index)
    
    genetic_distances_no_singletons=compare.delete_singletons_matrix(sample_data,genetic_distances)
       
    return(genetic_distances_no_singletons)


def get_pairwise_physical_distance(msprime_ts,sample_data):
    #create matrix of pairwise physical distances between mutations
    num_mutations = len(list(msprime_ts.variants())) 
    physical_distances = np.zeros((num_mutations,num_mutations))

    for outer_mutation in msprime_ts.mutations():
        for inner_mutation in msprime_ts.mutations():
            physical_distances[outer_mutation.index,inner_mutation.index] = abs(outer_mutation.position-inner_mutation.position)
    

    physical_distances_no_singletons=compare.delete_singletons_matrix(sample_data,physical_distances)
       
    return(physical_distances_no_singletons)

def mutation_pair(mutation_1, mutation_2, msprime_ts,physical_distance,genetic_distance):
    #1. how far apart are the two mutations
    genetic_distance = genetic_distance[mutation_1,mutation_2]
    physical_distance = physical_distance[mutation_1,mutation_2]
    return(genetic_distance,physical_distance)
    #how many significant trees are they apart
    
    #how many nodes are they apart in the second

def direct_comparison_pair_distance(msprime_ts, pairwise_matrix_binary,direct_comparisons,physical_distance,genetic_distance):
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

def all_misordered_pair_distance(msprime_ts, pairwise_matrix_binary,physical_distance,genetic_distance):
    #look up which pairs of mutations are in the wrong order in pairwise array (either from freq or geva)
    wrong_order = np.argwhere(pairwise_matrix_binary == 1)
    results = np.zeros((len(wrong_order), 2))
    
    #wrong_order first array is first mutation's index, second array has second mutations's indices
    for idx,(mut1_index,mut2_index) in enumerate(wrong_order):
        cur_result = mutation_pair(mut1_index,mut2_index,msprime_ts,physical_distance,genetic_distance)
        results[idx,0] = cur_result[0]
        results[idx,1] = cur_result[1]
        #results[idx,2] = cur_result[2]
        
    #return: 1. array of genetic distances of pairwise error; 2. array of physical distances of pairwise error
    return(results)

def all_pair_average_distance(msprime_ts,direct_comparisons,physical_distance,genetic_distance):
    #look up which pairs of mutations are in the wrong order in pairwise array (either from freq or geva)
    all_direct = np.argwhere((direct_comparisons == 1) | (direct_comparisons == 0))
    results = np.zeros((len(all_direct), 2))
    
    #all_direct first array is first mutation's index, second array has second mutations's indices
    for idx,(mut1_index,mut2_index) in enumerate(all_direct):
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
    direct_baseline = np.argwhere(direct_matrix == 1)
    results = np.zeros((len(direct_baseline), 5))

    #Iterate through each mutation pair that is directly comparable, recording distance
    #Check if correct by frequency/GEVA. For index assign: genetic distance,physical distance,freq_correct,geva_correct
    for idx,(mut1_index,mut2_index) in enumerate(direct_baseline):
        cur_result = mutation_pair(mut1_index,mut2_index,msprime_ts,physical_distance,genetic_distance)
        results[idx,0] = cur_result[0]
        results[idx,1] = cur_result[1]
        results[idx,2] = freq_matrix[mut1_index,mut2_index]
        results[idx,3] = geva_matrix[mut1_index,mut2_index]
    
    #Using the resulting numpy array of tuples, make regularly spaced bins of 500 bp and record accuracy
    #for each bin (number of 1's divded by 0's + 1's)
    sequence_length=msprime_ts.get_sequence_length()

    bins=np.arange(0,sequence_length,500)
    results[:,4]=np.digitize(results[:,1],bins)


    accuracy_averages=np.zeros((2,len(bins)))
    for idx,cur_bin in enumerate(bins):
        num_in_bin=len(results[results[:,4]==idx,:])
        freq_mismatch=sum(results[results[:,4]==idx,:][:,2]==1)
        geva_mismatch=sum(results[results[:,4]==idx,:][:,3]==1)
        freq_accuracy=div_zero(freq_mismatch,num_in_bin)
        geva_accuracy=div_zero(geva_mismatch,num_in_bin)
        accuracy_averages[0,idx]=freq_accuracy
        accuracy_averages[1,idx]=geva_accuracy

    return(accuracy_averages)
 

#helper function for preventing division by zero
def div_zero(x,y):
    if y == 0:
        return 0
    return x / y


