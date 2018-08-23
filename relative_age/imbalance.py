import msprime
import numpy as np

"""
Implement tree imbalance metrics: colless and sackin
Also output correlation with number of misordered mutation pairs
"""

def colless(simulated_ts,freq_matrix_no_single,geva_matrix_no_single,direct_matrix_no_single):
	colless_indices=list()
	for tree in simulated_ts.trees():
	    colless_index=0
	    for node in tree.nodes():
	        children=tree.get_children(node)
	        if children != ():
	            colless_index=colless_index+abs(len(list(tree.get_leaves(children[0])))-len(list(tree.get_leaves(children[1]))))
	    colless_indices.append(colless_index)
	return(colless_indices)


def sackin(simulated_ts,freq_matrix_no_single,geva_matrix_no_single,direct_matrix_no_single):
	for tree in simulated_ts.trees():
    sackin_stat=list()
    for leaf in list(tree.leaves()):
        cur_leaf=0
        cur_node=tree.parent(leaf)
        while cur_node != -1:
            cur_node=tree.parent(cur_node)
            cur_leaf=cur_leaf+1
        sackin_stat.append(cur_leaf)
    return(np.mean(sackin_stat),np.var(sackin_stat))

