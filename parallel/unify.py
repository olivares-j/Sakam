import sys
import os
import numpy as np
import pandas as pd
import h5py
from Globals import *


h5 = h5py.File(file_samples,'w')

#++++++++++++++++++ Loop over chunks +++++++++++++++++++++++++++
dfs = []
for i in range(1,size+1):
	#------------- Files ---------------------------------------------
	file_samp  = dir_chunks + "/samples_{0}_of_{1}.h5".format(i,size)
	file_stat  = dir_chunks + "/statistics_{0}_of_{1}.csv".format(i,size)
	#-----------------------------------------------------------------

	#------------- Load data frames into list --------------------------------
	dfi = pd.read_csv(file_stat)
	dfi.set_index(identifier,inplace=True)
	dfs.append(dfi)
	#---------------------------------------------------------

	#---------------- Read and save samples --------------------------
	with h5py.File(file_samp,'r') as hf:
		for ID in hf.keys():
			ogrp = hf.get(ID)
			ngrp = h5.create_group(ID)
			dset = ngrp.create_dataset("MAP",    data=ogrp.get("MAP"))
			dset = ngrp.create_dataset("Median", data=ogrp.get("Median"))
			dset = ngrp.create_dataset("SD",     data=ogrp.get("SD"))
			dset = ngrp.create_dataset("CI",     data=ogrp.get("CI"))
			dset = ngrp.create_dataset("sample", data=ogrp.get("sample"))
			h5.flush()
	#-----------------------------------------------------------------

h5.close()

#------- Concatenate data frames -------------------------
df = pd.concat(dfs)

#------------- Save statistics ---------------------
df.to_csv(file_statistics)
