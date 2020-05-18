
#------------ Load libraries -------------------
from __future__ import absolute_import, unicode_literals, print_function
import sys
import os
import numpy as np
import pandas as pd


#============ Data and Directories =============================
#-------Main directory ---------------
dir_main  = "/home/javier/Cumulos/Taurus/"
#-------------------------------------

#------------------------- Files ----------------------------------------------------------
file_dst  = dir_main + "Distances/GMM/Sources_mean.csv"  # Kalkayotl distances
file_pht  = dir_main + "Luismatron/members.csv"           # File name of the data set
file_out  = dir_main + "Sakam/Absolute_magnitudes.csv"     # Output file with absolute magnitudes
#-------------------------------------------------------------------------------------------

#--------- Bands --------------------------------------------------------
gaia_bands      = ["g","rp","bp"]
gaia_errors     = ["g_error","rp_error","bp_error"]
twomass_bands   = ["Jmag","Hmag","Kmag"]
twomass_errors  = ["e_Jmag","e_Hmag","e_Kmag"]
panstars_bands  = ["gmag","rmag","imag","zmag","ymag"]
panstars_errors = ["e_gmag","e_rmag","e_imag","e_zmag","e_ymag"]

bands  = sum([gaia_bands,twomass_bands,panstars_bands],[])
errors = sum([gaia_errors,twomass_errors,panstars_errors],[])
bes    = sum([bands,errors],[])

#-------------- Load data --------------------------------
dst = pd.read_csv(file_dst,usecols=["source_id","mean","lower","upper"])
pht = pd.read_csv(file_pht,usecols=sum([["source_id"],bes],[]))
#-----------------------------------------------------------------------------------------------

#--------- Join -----------------
dst.set_index("source_id",inplace=True)
pht.set_index("source_id",inplace=True)
df = pht.join(dst,on="source_id")


#--------- Absolute photometry and Uncertainty ------------
for band,error in zip(bands,errors):
	idx = np.where(np.isfinite(df[band]) & np.isnan(df[error]))[0]

	if len(idx)>0:
		print("Band ",band," has conflicting measurement and uncertainties!")
		print("Uncertainties will be assigned as 0.1 mag")
		df[error].iloc[idx] = 0.1

	df["abs_"+band]  = df[band]  - 5.0*np.log10(df["mean"])+5.
	df["abs_"+error] = df[error] + 5.0*(np.log10(df["upper"])-np.log10(df["lower"]))


#----------- Save ---------------------------
df.to_csv(file_out)
