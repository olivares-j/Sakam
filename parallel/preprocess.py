import sys
import os
import numpy as np
import pandas as pd
from Globals import *

#-------- Create directories ------------
os.makedirs(dir_data,exist_ok=True)
os.makedirs(dir_plots,exist_ok=True)
os.makedirs(dir_chunks,exist_ok=True)

#------------- Load data ----------------------------------
df = pd.read_csv(file_data,usecols=columns_data)
df.replace(to_replace=nan_values,value=np.nan,inplace=True)
df.set_index(identifier,inplace=True)
n_init = len(df)
print("The data set contains {0} sources.".format(n_init))
#-----------------------------------------------------------

#+++++++++++++++++++ Filter data ++++++++++++++++++++++++++
#---- Set as NaN the BP values larger than limit_BP -------
if label_BP is not None:
	idx = np.where(df[label_BP] > limit_BP)[0]
	if len(idx) > 0:
		df.loc[df.iloc[idx].index,label_BP] = np.nan
#----------------------------------------------------------

#---- Set as missing if band or uncertainty are missing ---
for ob,un in zip(observables,uncertainties):
	mask = np.isnan(df.loc[:,ob]) | np.isnan(df.loc[:,un])
	df.loc[mask,[ob,un]] = np.nan
	# print("{0}: {1}".format(ob,sum(mask)))
#----------------------------------------------------------

#--- Remove objects with less than n_obs_min bands --------
df.dropna(thresh=n_obs_min,subset=observables,inplace=True)
#----------------------------------------------------------

#---- Minimum uncertainty --------------------------------
for un in uncertainties:
	df.loc[:,un] += add_unc
#----------------------------------------------------------

print("After filtering {0} sources were removed.".format(n_init - len(df)))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++ Split data frame +++++++++++++++++++++++++++++++++
n_sources = len(df)
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
	grouped_df.get_group(g).to_csv(dir_chunks + "/data_{0}_of_{1}.csv".format(g,size))
