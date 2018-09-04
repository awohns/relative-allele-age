import msprime
import numpy as np
from collections import defaultdict
import itertools
import compare

"""
Implement tree imbalance metrics: colless and sackin
Also output correlation with number of misordered mutation pairs
"""


def colless_index(tree):
    """
    Calculate colless statistic for given tree
    """
    colless_index=0
    for node in tree.nodes():
        children=tree.get_children(node)
        if children != ():
            colless_index=colless_index+abs(len(list(tree.get_leaves(children[0])))-len(list(tree.get_leaves(children[1]))))
    return(colless_index)


def sackin_index(tree):
    """
    Calculate sackin index for given tree
    """
    #add the number of internal nodes between each leaf of the tree and the root (inclusive)
    sackin_stat=list()
    for leaf in list(tree.leaves()):
        cur_leaf=0
        cur_node=tree.parent(leaf)
        while cur_node != -1:
            cur_node=tree.parent(cur_node)
            cur_leaf=cur_leaf+1
        sackin_stat.append(cur_leaf)
    return(np.mean(sackin_stat),np.var(sackin_stat))

def make_imbalance_dictionary(ts):
    """
    Create a dictionary of tree imbalance for each tree in ts
    """
    imbalance_dictionary = defaultdict(dict)
    for tree in ts.trees():
        imbalance_dictionary[tree.index] = {"colless":colless_index(tree),"sackin":sackin_index(tree)}

    return(imbalance_dictionary)

def all_avg_tree_imbalance(imbalance_dictionary, ts, sample_file):
    """
    Output average tree imbalance for all pairs 
    """
    imbalance_list = list()
    
    num_vars=len(sample_file.sites_genotypes[:])
    mutation_pairs=itertools.combinations(range(0,num_vars),2)
    
    singletons=compare.find_singletons_sample(sample_file)
   
    for var1,var2 in mutation_pairs:
        
        if (var1 not in singletons) and (var2 not in singletons):
            for tree in ts.trees():
                for tree_mut in tree.mutations():
                    if (tree_mut.index == var1) or (tree_mut.index == var2):
                        imbalance_list.append(imbalance_dictionary[tree.index])

    return(np.mean([stat["sackin"][0] for stat in imbalance_list]),
        np.mean([stat["sackin"][1] for stat in imbalance_list]),
        np.mean([stat["colless"] for stat in imbalance_list]))
  
def binary_matrix_avg_tree_imbalance(imbalance_dictionary, ts, sample_file, binary_matrix_with_singletons,second_binary_matrix=None):
    """
    Output average tree imbalance given a binary matrix with singletons
    """
    imbalance_list=list()
    
    if second_binary_matrix is not None:
        directly_comparable_errors=np.add(binary_matrix_with_singletons,second_binary_matrix)
        binary_matrix_triu=np.triu(directly_comparable_errors)
        direct_pairs=np.where(directly_comparable_errors==2)
    else:
        binary_matrix_triu=np.triu(binary_matrix_with_singletons)
        direct_pairs=np.where(binary_matrix_triu==1)
    
    singletons=compare.find_singletons_sample(sample_file)
 
    for var1,var2 in zip(direct_pairs[0],direct_pairs[1]):
    
        if (var1 not in singletons) and (var2 not in singletons):
            for tree in ts.trees():
                for tree_mut in tree.mutations():
                    if (tree_mut.index == var1) or (tree_mut.index == var2):
                            imbalance_list.append(imbalance_dictionary[tree.index])
        
        
  
    #singletons=compare.find_singletons_sample(sample_file)
    #freq_misordered_no_singletons=[x for i,x in enumerate(freq_misordered) if i not in singletons]
    return(np.mean([stat["sackin"][0] for stat in imbalance_list]),
        np.mean([stat["sackin"][1] for stat in imbalance_list]),
        np.mean([stat["colless"] for stat in imbalance_list]))

