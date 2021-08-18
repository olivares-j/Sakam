import sys
import os
import numpy as np
import pandas as pd
import dill

#----- Global variables ----------
path_globals  = str(sys.argv[1]) 
dill.load_session(path_globals)
#---------------------------------

#------------- Load data ----------------------------------
df = pd.read_csv(file_abs,usecols=columns_abs)
df.replace(to_replace=nan_values,value=np.nan,inplace=True)
df.set_index(identifier,inplace=True)
n_sources = len(df)
#-----------------------------------------------------------

#---------------- Split data frame --------------------------------
assert n_sources => size, "There are more processes than sources!"
group_size = int(np.floor(n_sources/size))
reminder = n_sources % size
group_size = np.repeat(group_size,size)
group_size[-1] += reminder
groups = []
for g,gs in enumerate(group_size):
	groups.append(np.repeat(g+1,gs))

groups = np.concatenate(groups)

df.insert(loc=0,column="Groups",value=groups)
grouped_df = df.groupby("Groups")

#--- Write each chunk -----
for g in range(1,size+1):
	grouped_df.get_group(g).to_csv(dir_chunks + "data_{0}_of_{1}.csv".format(g,size))
