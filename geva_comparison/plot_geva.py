import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
import os

# master_dataframe = pd.DataFrame(0,columns=["Agree", "Total", "Error Agree", "Error Total", "Agree Geva","Agree Error Geva"],index=[0,1,2,3,4])


# for filename in os.listdir("/home/wilderwohns/geva_comparison/data"):
#     master_dataframe=pd.DataFrame.add(pd.read_csv("/home/wilderwohns/geva_comparison/data/"+filename),master_dataframe)

master_dataframe=pd.read_csv("/home/wilderwohns/geva_comparison/data/GEVA_frequency_distance_accuracy.csv")

#Without SINGLETONS, don't consider any of the same age, flip coin for same frequency
plt.plot(master_dataframe.Agree/master_dataframe.Total,label="No Error")

plt.plot(master_dataframeAgreeError/master_dataframeTotalError,label="Error")
plt.plot(master_dataframe.AgreeGeva/master_dataframe.TotalGeva,label="Agree GEVA")
plt.plot(master_dataframe.AgreeGevaError/master_dataframe.TotalGevaError,label="Agree Error GEVA")
# plt.xlabel("Distance Separating Alleles (bp)")

plt.xlabel("Kb separating Alleles")
plt.ylabel("Proportion of Mutation Pairs Correctly Ordered")

plt.legend()
plt.savefig("geva_comparison_plot")