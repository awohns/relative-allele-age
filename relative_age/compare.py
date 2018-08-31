import subprocess
import pandas as pd
import numpy as np
import math
import msprime
import tsinfer
import logging
import warnings
warnings.filterwarnings("ignore")

"""
This module perform comparisons
based on frequency and GEVA
"""


def find_age(marker_id, age_estimates):
    if (age_estimates.iloc[:,0] == marker_id).any():
        return(float(age_estimates[(age_estimates.iloc[:,0] == marker_id) & (age_estimates.iloc[:,1] == "C") & (age_estimates.iloc[:,2] == 1)]["PostMode"]))
   


def find_singletons_sample(sample_data):
    singleton_list = list()
    for i,genos in enumerate(sample_data.sites_genotypes[:]):
        if len(genos[genos ==1]) == 1:
            singleton_list.append(i)
        if len(genos[genos == 1]) == (len(genos)-1):
            singleton_list.append(i)
        if 0 not in singleton_list:
            singleton_list.append(0)
        last_value = len(list(sample_data.sites_genotypes[:])) - 1
        if last_value not in singleton_list:
            singleton_list.append(last_value)
    return(singleton_list)


def delete_singletons_matrix(sample_data,mut_matrix):
    singletons = find_singletons_sample(sample_data)
    mut_matrix_no_singletons= np.delete(mut_matrix, singletons, 0)
    mut_matrix_no_singletons= np.delete(mut_matrix_no_singletons, singletons, 1)
    return(mut_matrix_no_singletons)

"""
Frequency comparisons using samples object
"""



#relative time orderings frequency
def freq_relative_time(msprime_ts,sample_data):
    num_mutations = len(sample_data.sites_genotypes[:])
    pairwise_matrix_frequency = np.zeros((num_mutations,num_mutations))
    variants=list(msprime_ts.variants())
    
    for i,outer_geno in enumerate(sample_data.sites_genotypes[:]):
        outer_age = int(msprime_ts.node(variants[i].site.mutations[0].node).time)
        outer_freq = len(outer_geno[outer_geno ==1])/len(outer_geno)
        for j,inner_geno in enumerate(sample_data.sites_genotypes[:]):
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
    
    
    pairwise_matrix_frequency_no_singletons=delete_singletons_matrix(sample_data,pairwise_matrix_frequency)
    
    return(pairwise_matrix_frequency,pairwise_matrix_frequency_no_singletons)




"""
GEVA comparisons
"""

def geva_age_estimate(file_name,ne,mut_rate):
    subprocess.check_output(["../bin/rvage_dev3/./rvage_dev3", '--vcf', file_name+".vcf", "-o","../tmp/"+file_name], cwd=r'/home/wilderwohns/relative_allele_age/scripts')
    with open("../tmp/"+file_name+".positions.txt","wb") as out:
        subprocess.call(["awk", "NR>3 {print last} {last = $3}", "../tmp/"+file_name+".marker.txt"], stdout=out)
    try:
        subprocess.check_output(["../bin/rvage_dev3/./rvage_dev3", "age", "-i", "../tmp/"+ file_name+".bin", "--positions", "../tmp/"+file_name+".positions.txt", "-m","hmm","--hmm","../bin/rvage_dev3/_initials.hhmm.tru","../bin/rvage_dev3/_emission.hhmm.tru.100","--Ne",'10000', "--mut", '2e-8',"--selectNN","1","--maxConcordant","100","--maxDiscordant", "100","-o","../tmp/"+file_name+"_estimation"])
    except subprocess.CalledProcessError as grepexc:
        print("error file:" + str(file_name))
        print(grepexc.output)
    age_estimates = pd.read_csv("../tmp/"+file_name+"_estimation.cle.txt", sep = " ")
    return(age_estimates)




#I compared with geva_age_estimate and it's pretty close
def geva_sample_estimate_age(file_name,sample_data,Ne,length,mut_rate):
    
    num_individuals = len(sample_data.individuals_metadata[:])
    ind_list = list()
    pos_geno_dict = {"POS":list()}

    for i in range(int(num_individuals/2)):
        pos_geno_dict["msp_"+str(i)] = list()
        ind_list.append("msp_"+str(i))

    #add all the sample positions and genotypes
    for i in sample_data.genotypes():
        cur_pos = int(round(sample_data.sites_position[i[0]]))
        if cur_pos in pos_geno_dict["POS"]:
            pos_geno_dict["POS"].append(cur_pos+1)
        else:
            pos_geno_dict["POS"].append(cur_pos)

        for j in range(0,len(i[1]),2):

            pos_geno_dict["msp_"+str(int(j/2))].append(str(i[1][j]) + "|" + str(i[1][j+1]))

    #Deal with duplicated position numbers
    df = pd.DataFrame(pos_geno_dict)
    while sum(df.POS.duplicated()) > 0:
        dup_row=df[df.POS.duplicated()].tail(1)
        dup_pos=dup_row.POS
        df.at[dup_row.index, 'POS'] = dup_pos+1

    df["#CHROM"] = 1
    df["REF"] = "A"
    df["ALT"] = "T"
    df['ID'] = "."
    df['QUAL'] = "."
    df['FILTER'] = "PASS"
    df['INFO'] = "."
    df['FORMAT'] = "GT"

    cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO','FORMAT']+ind_list
    df = df[cols] 

    header = """##fileformat=VCFv4.2
##source=msprime 0.6.0
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=1,length=""" + str(int(length)) +  """>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
"""
    
    output_VCF = "../tmp/"+file_name+"error.vcf"
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)

    df.to_csv(output_VCF, sep="\t", mode='a',index=False)
    age_estimates = geva_age_estimate("../tmp/"+file_name,10000,2e-8)
    return(age_estimates)



#removes singletons
def geva_all_time_orderings(file_name,msprime_ts,sample_data,Ne,length,mutation_rate):
    num_mutations = len(list(msprime_ts.variants())) 
    pairwise_matrix_geva = np.zeros((num_mutations,num_mutations))

    #age estimation
    age_estimates = geva_sample_estimate_age(file_name,sample_data,Ne,length,mutation_rate)

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
    
    pairwise_matrix_geva_no_singletons=delete_singletons_matrix(sample_data,pairwise_matrix_geva)    
        
    return(pairwise_matrix_geva,pairwise_matrix_geva_no_singletons,age_estimates)



def get_geva_corrected(geva,freq,direct_comparsion):
    geva_inverse=1-geva
    freq_direct_comparison = np.add(freq,direct_comparsion)
    geva_corrected = np.add(freq_direct_comparison,geva_inverse) 
    geva_corrected = np.add((geva_corrected == 3 ).astype(int),direct_comparsion)
    num_wrong_freq_direct = len(freq_direct_comparison[freq_direct_comparison==2])
    if num_wrong_freq_direct ==0:
        return(None)
    else:
        return(len(geva_corrected[geva_corrected == 2])/num_wrong_freq_direct)


#determine which mutations are "directly comparable"
#This means, which pair of mutations share a node as an ancestor or descendant
#if mutations have the same node as an ancestor or descedent, then they are directly comparable
#in other words, we know the relative time ordering (topology), locally
def directly_comparable(msprime_ts,sample_data):
    total_yes = 0
    #get the total number of variants and make a zero np array of that shape
    num_mutations = len(list(msprime_ts.variants()))
    pairwise_matrix_same_tree = np.zeros((num_mutations,num_mutations))

    #loop through all trees in the tree sequence
    for tree in msprime_ts.trees():
        #for each mutation in each tree
        for cur_mutation in tree.mutations():
            cur_node = cur_mutation.node
            #look up the tree
            while cur_node != tree.root:
                cur_node = tree.parent(cur_node)
                #check all the mutations in the tree sequence to see if they are directly above this node
                for mutation in msprime_ts.mutations():
                    if mutation.node == cur_node:
                        pairwise_matrix_same_tree[cur_mutation.index,mutation.index] = 1
                        pairwise_matrix_same_tree[mutation.index,cur_mutation.index] = 1
    
    pairwise_matrix_same_tree_no_singletons=delete_singletons_matrix(sample_data,pairwise_matrix_same_tree) 
   
    return(pairwise_matrix_same_tree,pairwise_matrix_same_tree_no_singletons)


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


def sackin_index(ts):
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
    return(np.var(sackin_stat))

def tree_imbalance_ts(ts, samples, freq, geva, direct):
    for tree in simulated_ts.trees():
        colless_index=0
        colless_index.append(colless_index(tree))


    for tree in simulated_ts.trees():
        sackin_stat=list()
        sackin_stat.append(sackin_index(tree))

