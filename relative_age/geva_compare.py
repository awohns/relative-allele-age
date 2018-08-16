import subprocess
import pandas as pd

"""
GEVA comparisons
"""

def geva_age_estimate(file_name,ne,mut_rate):
    subprocess.check_output(["../bin/rvage_dev3/./rvage_dev3", '--vcf', file_name+".vcf", "-o","../tmp/"+file_name], cwd=r'/Users/anthonywohns/Documents/mcvean_group/relative_allele_age/tests')
    with open("../tmp/"+file_name+".positions.txt","wb") as out:
        subprocess.call(["awk", "NR>3 {print last} {last = $3}", "../tmp/"+file_name+".marker.txt"], stdout=out)
    try:
        subprocess.check_output(["../bin/rvage_dev3/./rvage_dev3", "age", "-i", "../tmp/"+ file_name+".bin", "--positions", "../tmp/"+file_name+".positions.txt", "-m","hmm","--hmm","../bin/rvage_dev3/_initials.hhmm.tru","../bin/rvage_dev3/_emission.hhmm.tru.100","--Ne",'10000', "--mut", '2e-8',"--selectNN","1","--maxConcordant","100","--maxDiscordant", "100","-o","../tmp/"+file_name+"_estimation"])
    except subprocess.CalledProcessError as grepexc:
        print(grepexc.output)
    age_estimates = pd.read_csv("../tmp/"+file_name+"_estimation.cle.txt", sep = " ")
    return(age_estimates)


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




#I compared with geva_age_estimate and it's good
def geva_with_error(sample_data,Ne,length,mut_rate):
    
    num_individuals = len(sample_data.individuals_metadata[:])
    ind_list = list()
    pos_geno_dict = {"POS":list()}

    for i in range(int(num_individuals/2)):
        pos_geno_dict["msp_"+str(i)] = list()
        ind_list.append("msp_"+str(i))

    #add all the sample positions and genotypes
    for i in sample_data.genotypes():
        pos_geno_dict["POS"].append(int(round(sample_data.sites_position[i[0]])))
        for j in range(0,len(i[1]),2):

            pos_geno_dict["msp_"+str(int(j/2))].append(str(i[1][j]) + "|" + str(i[1][j+1]))

    df = pd.DataFrame(pos_geno_dict)

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
    
    output_VCF = "../tmp/error.vcf"
    with open(output_VCF, 'w') as vcf:
        vcf.write(header)

    df.to_csv(output_VCF, sep="\t", mode='a',index=False)
    
    return(geva_age_estimate("../tmp/error",10000,2e-8))



#removes singletons
def geva_all_time_orderings(msprime_ts,sample_data,Ne,length,mutation_rate,delete_singletons):
    num_mutations = len(list(msprime_ts.variants())) 
    pairwise_matrix_geva = np.zeros((num_mutations,num_mutations))

    #age estimation
    age_estimates = geva_with_error(sample_data,Ne,length,mutation_rate)

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
