import sys
import numpy as np
import pandas as pd
from astropy import units as u
from astroquery.xmatch import XMatch
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy.table import Table

#=============== Data =======================================================
#------------ ComaBer ------------------------------------------------------------
# dir_base = "/home/jolivares/OCs/ComaBer/"
# file_in  = dir_base + "Mecayotl/iter_5_all/Classification/members_mecayotl.csv"
# file_out = dir_base + "Sakam/iter_5/data/members_GaiaEDR3+2MASS+PanSTARRS.csv"
#---------------------------------------------------------------------------------

#------------------- Taurus -----------------------------------------------------
dir_base = "/home/jolivares/OCs/Taurus/GEDR3/runs/d/"
file_in  = dir_base + "tables/members_bins.csv"
file_out = dir_base + "sakam/members_GaiaEDR3+2MASS+PanSTARRS.csv"
#--------------------------------------------------------------------------------

#-------- Catalogues --------------------------------------------------
catalogues = [
# Gaia EDR3
{"name":"GaiaEDR3","reference":"vizier:I/350/gaiaedr3","filters":{},
	"columns":['phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag',
	'phot_g_mean_mag_error','phot_bp_mean_mag_error', 'phot_rp_mean_mag_error'],
	"max_distance":5.0*u.arcsec},

# 2MASS
{"name":"2MASS","reference":"vizier:II/246/out","filters":{},
	"columns":['2MASS','Jmag','Hmag','Kmag','e_Jmag','e_Hmag','e_Kmag'],
	"max_distance":5.0*u.arcsec},

# PanSTARRS
{"name":"PanSTARRS","reference":"vizier:II/349/ps1","filters":{"Ns":0},
	"columns":['objID','gmag','rmag','imag','zmag','ymag',
				'e_gmag','e_rmag','e_imag','e_zmag','e_ymag'],
	"max_distance":5.0*u.arcsec}
]
#----------------------------------------------------------------------

#=======================================================================


#--------- Read members ------------------
df = pd.read_csv(file_in,
		usecols=["source_id","ra","dec"])
#--------------------------------------

#----- Rename to avoid duplicates --------
df.rename(columns={
		"source_id":"input_source_id",
		},inplace=True)
#----------------------------------------

print("Input members: {0}".format(df.shape[0]))

#---------- Coordinates J2000 -----------------------
coords_icrs = SkyCoord(
				ra=df["ra"]*u.degree,
				dec=df["dec"]*u.degree,
				frame='icrs')
coords_J2000 = coords_icrs.fk5
#----------------------------------------------

#------- Append new & drop old ----------------
df["input_RAJ2000"] = coords_J2000.ra.deg
df["input_DEJ2000"] = coords_J2000.dec.deg
df.drop(columns=["ra","dec"],inplace=True)
#-----------------------------------------------

#------ Transform to Table ---------
tab = Table.from_pandas(df)
#-----------------------------------

#------ Rename -----------------------
df.rename(columns={
	"input_source_id":"source_id",
	"input_RAJ2000":"RAJ2000",
	"input_DEJ2000":"DEJ2000"
	},inplace=True)
#-------------------------------------

#--------- set index -----------------
df.set_index("source_id",inplace=True)
#-------------------------------------

for cat in catalogues:
	#---------- Match ---------------------
	result = XMatch.query(
				cat1=tab,
                cat2=cat["reference"],
                max_distance=cat["max_distance"],
                colRA1='input_RAJ2000',
                colDec1='input_DEJ2000'
                )
	#--------------------------------------

	#-------- To pandas -----
	tmp = result.to_pandas()
	#------------------------

	#--------- Filter -----------------------
	for key,value in cat["filters"].items():
		tmp = tmp.loc[tmp[key] > value]
	#----------------------------------------
	
	#-------- Drop duplicates ------------------
	tmp.drop_duplicates(
				subset=["input_source_id"],
				keep="first",inplace=True)
	#-------------------------------------------

	#----- Gaia ------------------------------------
	if cat["name"] == "GaiaEDR3":
		tmp_bad = tmp.loc[tmp["input_source_id"] != tmp["source_id"]]
		if len(tmp_bad)>0:
			print("The following Gaia catalogue sources " + \
			       " do not match the input ones:\n")
			print(tmp_bad[["input_source_id","source_id","angDist"]])
			# sys.exit()

		tmp.drop(columns=["source_id"],inplace=True)
	#-----------------------------------------------

	print("Members in {0}: {1}".format(cat["name"],
										tmp.shape[0]))

	#------- Select columns ----------------
	columns = ["input_source_id","angDist"]
	columns.extend(cat["columns"])
	tmp = tmp[columns]
	#---------------------------------------

	#----- Rename angular distance -------------
	tmp.rename(columns={
			"angDist":"angDist_"+cat["name"]},
			inplace=True)
	#------------------------------------------

	#---- Set index ----------------
	tmp.set_index("input_source_id",
		inplace=True)
	#-------------------------------

	#-------- Merge --------------
	df = df.merge(tmp,
		left_index=True,
		right_index=True,
		how="left")
	#----------------------------


#-------- Save -----------------
df.to_csv(file_out,
	index_label="source_id")
#-------------------------------
