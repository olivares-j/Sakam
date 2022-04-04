import sys
import numpy as np
import pandas as pd

#---------- AMES-COND --------------------------------------------------
# dir_    = "/home/jolivares/OCs/Models/AMES-Cond/"

# instruments  = [
# 		{"name":"2MASS",
# 		"file":"model.AMES-Cond-2000.M-0.0.2MASS.Vega",
# 		"columns":["Age/Myr","M/Ms","Teff(K)","L/Ls","Li","J","H","K"]},
# 		{"name":"Gaia",
# 		"file":"model.AMES-Cond-2000.M-0.0.GAIA.Vega",
# 		"columns":["Age/Myr","M/Ms","G","G_BP","G_RP"]},
# 		{"name":"SDSS",
# 		"file":"model.AMES-Cond-2000.M-0.0.SDSS.AB",
# 		"columns":["Age/Myr","M/Ms","g","r","i","z"]},
# 		]
#-----------------------------------------------------------------------

#-------- BT-COND ------------------------------------------------------
dir_    = "/home/jolivares/OCs/Models/BT-Cond/"

instruments  = [
		{"name":"2MASS",
		"file":"model.BT-Cond.M-0.0.2MASS.Vega",
		"columns":["Age/Myr","M/Ms","Teff(K)","L/Ls","Li","J","H","K"]},
		{"name":"Gaia",
		"file":"model.BT-Cond.M-0.0.GAIA.Vega",
		"columns":["Age/Myr","M/Ms","G","G_BP","G_RP"]},
		{"name":"Subaru",
		"file":"model.BT-Cond.M-0.0.SUBARU.AB",
		"columns":["Age/Myr","M/Ms","g","r","i","z","Y"]},
		{"name":"UKIDSS",
		"file":"model.BT-Cond.M-0.0.UKIDSS.AB",
		"columns":["Age/Myr","M/Ms","y"]}
		]
#----------------------------------------------------------------------

file_out = dir_ + "+".join([instrument["name"] for instrument in instruments])+".csv"
#------------------------------------------------------------------------

#------------ Load isochrones --------------------
dfm = []
for instrument in instruments:
	dfs = []
	begin_file = True
	with open(dir_ + instrument["file"]) as file:
		for line in file.read().splitlines():
			#------- Skip usless lines --------------------
			if line.startswith("-") or line.strip() == '':
				continue
			#---------------------------------------------

			#------- Read ages -------------------------------
			if line.startswith("t",3):
				if not begin_file:
					data = np.asarray(data).astype(np.float32)
					df = pd.DataFrame(data=data,columns=hdr)
					df["Age"] = age
					dfs.append(df)

				begin_file = False
				data = []
				age  = float(line[13:])
			#-------------------------------------------------

			#------- Read header ----------------
			if line.startswith("M",4):
				line = line[:8]+str(" ")+line[8:]
				hdr = line.split()
			#------------------------------------

			#------ Read columns -----------
			if not line.startswith(" ",2):
				data.append(line.split())
			#-------------------------------

	#---------- Last block --------------------
	data = np.asarray(data).astype(np.float32)
	df = pd.DataFrame(data=data,columns=hdr)
	df["Age"] = age
	dfs.append(df)
	#------------------------------------------

	#----- Concatenate -------------------
	df = pd.concat(dfs,ignore_index=True)
	#-------------------------------------

	#------ Transform and extract -------
	df["Age/Myr"] = df["Age"] * 1.e3
	df = df.loc[:,instrument["columns"]]
	#------------------------------------

	#-- Append -------
	dfm.append(df)
	#----------------

#----------- Join -----------------------------------------
df = dfm[0]
for tmp in dfm[1:]:
	df = df.merge(tmp,on=["Age/Myr","M/Ms"],how="outer")
#----------------------------------------------------------

#-------------- Write ------------
df.to_csv(file_out,index=False)
#--------------------------------