#------------ Load libraries -------------------
from __future__ import absolute_import, unicode_literals, print_function
import sys
import dill
import numpy as np
import pandas as pd

#----- Global variables ----------
path_globals  = str(sys.argv[1]) 
dill.load_session(path_globals)
#---------------------------------

#-------------- Load distance --------------------------------
dst = pd.read_csv(file_distances,usecols=columns_distances)
dst.set_index(identifier,inplace=True)

if kalkayotl_dimension != 1:
	def cartesian_distance(x,y,z):
		return np.sqrt(x**2 + y**2 + z**2)

	dst["lower_distance"] = dst[lowers].apply(lambda x: cartesian_distance(*x),axis=1)
	dst["upper_distance"] = dst[uppers].apply(lambda x: cartesian_distance(*x),axis=1)

	dst.drop(columns=sum([lowers,uppers],[]),
				inplace=True)
else:
	dst.rename(columns={"lower":"lower_distance",
						"upper":"upper_distance",
						"mean":"mode_distance"},
				inplace=True)
#----------------------------------------------------------------------------------

#----------- Load photometry ----------------------------------
pht = pd.read_csv(file_photometry,usecols=columns_photometry)
pht.set_index(identifier,inplace=True)
#--------------------------------------------------------------

#--------- Join -----------------------
df = pht.join(dst,on=identifier)
n_init = len(df)
print("The data set contains {0} sources.".format(n_init))
#--------------------------------------

#+++++++++++++++++++ Filter data ++++++++++++++++++++++++++
#---- Set as NaN the BP values larger than limit_BP -------
if label_BP is not None:
	idx = np.where(df[label_BP] > limit_BP)[0]
	if len(idx) > 0:
		df.loc[df.iloc[idx].index,label_BP] = np.nan
#----------------------------------------------------------

#---- Set as missing if band or uncertainty are missing ---
for ob,un in zip(bands,errors):
	mask = np.isnan(df.loc[:,ob]) | np.isnan(df.loc[:,un])
	df.loc[mask,[ob,un]] = np.nan
	# print("{0}: {1}".format(ob,sum(mask)))
#----------------------------------------------------------

#--- Remove objects with less than n_obs_min bands --------
df.dropna(thresh=n_obs_min,subset=bands,inplace=True)
#----------------------------------------------------------

#---- Minimum uncertainty --------------------------------
for un in errors:
	df.loc[:,un] += add_unc
#----------------------------------------------------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#--------- Absolute photometry and Uncertainty ------------
for band,error in zip(bands,errors):
	idx = np.where(np.isfinite(df[band]) & np.isnan(df[error]))[0]

	if len(idx)>0:
		print("Band ",band," has conflicting measurement and uncertainties!")
		print("Uncertainties will be assigned as 0.1 mag")
		df[error].iloc[idx] = 0.1

	df[prefix + band]  = df[band]  - 5.0*np.log10(df["mode_distance"])+5.
	df[prefix + error] = df[error] + 5.0*np.abs(np.log10(df["upper_distance"])-np.log10(df["lower_distance"]))

df.drop(columns=["lower_distance","mode_distance","upper_distance"],inplace=True)
#-----------------------------------------------------------------------

print("The filtered data set contains {0} sources.".format(df.shape[0]))
#----------- Save ---------------------------
df.to_csv(file_abs)
