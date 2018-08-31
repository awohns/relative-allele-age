import msprime
import numpy as np
from collections import defaultdict

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
    for tree in simulated_ts.trees():
        imbalance_dictionary[tree.index] = {"colless":colless_index(tree),"sackin":sackin_index(tree)}

    return(imbalance_dictionary)

def all_average_tree_imbalance(imbalance_dictionary, ts):
    """
    Output average tree imbalance for all pairs 
    """
    imbalance_list = list()
    for mutation in ts.variants():
        for tree in simulated_ts.trees():
            for tree_mut in tree.mutations():
                if (tree_mut.index == mutation.index):
                    imbalance_list.append(imbalance_dictionary[tree.index])
                    break
    return(np.mean([stat["sackin"][0] for stat in imbalance_list]),
        np.mean([stat["sackin"][1] for stat in imbalance_list]),
        np.mean([stat["colless"] for stat in imbalance_list]))

  
