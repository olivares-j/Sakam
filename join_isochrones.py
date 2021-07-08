import sys
import dill
import numpy as np
import pandas as pn
from scipy.interpolate import splrep,splev

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#---------- Input files --------------------------------------------------
dir_    = "/home/jolivares/Cumulos/Models/"
dir_out = "/home/jolivares/Cumulos/Perseus/Models/"

variate = "Mass" #"Luminosity"

models  = [
		{
		"name":"PARSEC",
		"file":dir_ + "PARSEC/PARSEC_1-10_Myr_Gaia+2MASS+PanSTARRS.csv",
		"columns":["age_Myr","Mass","logL",
					"Gmag","G_BPftmag","G_RPmag", #Gaia
					"gP1mag","rP1mag","iP1mag","zP1mag","yP1mag", #PanSTARRS
					"Jmag","Hmag","Ksmag" #2MASS
					],
		"mapper":{"Mass":"Mass", "logL":"Luminosity",
				"Gmag":"Gmag","G_BPftmag":"G_BPmag","G_RPmag":"G_RPmag", # Gaia
				"gP1mag":"gP1mag","rP1mag":"rP1mag","iP1mag":"iP1mag","zP1mag":"zP1mag","yP1mag":"yP1mag", #PanSTARRS
				"Jmag":"Jmag","Hmag":"Hmag","Ksmag":"Kmag" # 2MASS
				},
		"lower_mass":1.0,"upper_mass":10.0,
		},
		{
		"name":"MIST",
		"file":dir_ + "MIST/MIST_1-10_Myr_Gaia+2MASS+PanSTARRS.csv",
		"columns":["age_Myr","initial_mass","log_L",
					"Gaia_G_EDR3","Gaia_BP_EDR3","Gaia_RP_EDR3", #Gaia
					"PS_g","PS_r","PS_i","PS_z","PS_y", #PanSTARRS
					"2MASS_J","2MASS_H","2MASS_Ks" #2MASS
					],
		"mapper":{"initial_mass":"Mass", "log_L":"Luminosity",
				"Gaia_G_EDR3":"Gmag","Gaia_BP_EDR3":"G_BPmag","Gaia_RP_EDR3":"G_RPmag", # Gaia
				"PS_g":"gP1mag","PS_r":"rP1mag","PS_i":"iP1mag","PS_z":"zP1mag","PS_y":"yP1mag", #PanSTARRS
				"2MASS_J":"Jmag","2MASS_H":"Hmag","2MASS_Ks":"Kmag" # 2MASS
				},
		"lower_mass":1.0,"upper_mass":10.0,
		},
		{
		"name":"BT-Settl",
		"file":dir_ + "BT-Settl/BT-Settl_all_Myr_Gaia+2MASS+PanSTARRS.csv",
		"columns":["age_Myr","M/Ms","L/Ls",
					"G","G_BP","G_RP", #Gaia
					"g_p1","r_p1","i_p1","z_p1","y_p1", #PanSTARRS
					"J","H","K" #2MASS
					],
		"mapper":{"M/Ms":"Mass","L/Ls":"Luminosity",
				"G":"Gmag","G_BP":"G_BPmag","G_RP":"G_RPmag", # Gaia
				"g_p1":"gP1mag","r_p1":"rP1mag","i_p1":"iP1mag","z_p1":"zP1mag","y_p1":"yP1mag", #PanSTARRS
				"J":"Jmag","H":"Hmag","K":"Kmag" # 2MASS
				},
		"lower_mass":0.0,"upper_mass":1.0,
		}]
ages   = [3,5,7,10] # Myr
#------------------------------------------------------------------------

#---------- Output file ------------------------------
base_ = "PMB"
steps = 5000
file_plot = dir_out + "join.pdf"
#----------------------------------------------------------

#--------------------- Isochrone -----------------------------------------
stol_v2c = 0.5
stol_c2v = 1e-3
degspl = int(3)
steps_bt = 150
c2v_drop = [2.3,2.4]
covariates = [
			"G_BPmag","Gmag"
			# "Gmag","G_BPmag","G_RPmag", # Gaia
			# "gP1mag","rP1mag","iP1mag","zP1mag","yP1mag", #PanSTARRS
			# "Jmag","Hmag","Kmag" # 2MASS
			  ]
#----------------------------------------------------------------------

#------------ Load isochrones --------------------
dfm = []
for model in models:
	#------- Read input file ---------------------------------
	tmp = pn.read_csv(model["file"], usecols=model["columns"])
	#---------------------------------------------------------

	#-------------- Rename -----------------------------------
	tmp.rename(columns=model["mapper"],inplace=True)
	#---------------------------------------------------------

	#------- Ages to integers ------------------
	tmp["age_Myr"] = tmp["age_Myr"].astype(int)
	#-------------------------------------------

	#-------- Drop ages -------------------------------------------
	tmp.drop(tmp[~np.isin(tmp["age_Myr"],ages)].index,inplace=True)
	#--------------------------------------------------------------

	#-------- Drop masses outside limit ---------------------------------
	tmp.drop(tmp[(tmp["Mass"] > model["upper_mass"]) | 
				 (tmp["Mass"] < model["lower_mass"])].index,inplace=True)
	#--------------------------------------------------------------------

	#------- Add model and weight ----------------------------------------
	tmp["Model"] = model["name"]
	#--------------------------------------------------------------------

	dfm.append(tmp)

#----- Join data frames and sort them by mass --
dfm = pn.concat(dfm,sort=False)
dfm.sort_values(by=variate,inplace=True)
#-----------------------------------------------

#------ Open plot -------
pdf = PdfPages(file_plot)
#------------------------

#================ Loop over ages ================================================
#------ Group by ages ------------------
dfga = dfm.groupby(by="age_Myr")
#-----------------------------------------

for age in ages:
	#--- Get age data --------
	dfa = dfga.get_group(age)
	#-------------------------

	#---------- Min Max --------------
	max_variate  = dfa[variate].max()
	min_variate  = dfa[variate].min()
	#---------------------------------

	#------ Group by models ------------------
	dfgm = dfa.groupby(by="Model")
	#-----------------------------------------

	#--------- Create output grid ------------------------
	if variate == "Mass":
		dfo = pn.DataFrame(data={variate:np.logspace(start=np.log10(min_variate),
							stop=np.log10(max_variate),num=steps)})
	else:
		dfo = pn.DataFrame(data={variate:np.linspace(start=min_variate,
							stop=max_variate,num=steps)})
	#-----------------------------------------------------

	#========== Loop over magnitude =======================================
	tcks_v2c = []
	tcks_c2v = []
	for covariate in covariates:
		#----- Initialize plot ----
		plt.figure()
		plt.suptitle(str(age) + " Myr")
		#--------------------------

		#----- Plot original models -----------------
		for g in dfgm.groups:
			dfg = dfgm.get_group(g)
			plt.scatter(dfg[variate],dfg[covariate],
						s=1,zorder=0,label=g)
		#-------------------------------------------

		#======= Interpolate BT-Settl ==========
		dfb = dfgm.get_group("BT-Settl")
		#---------- Min Max --------------
		bmax = dfb[variate].max()
		bmin = dfb[variate].min()
		#---------------------------------

		#--------- Grid ------------------
		xb = np.logspace(start=np.log10(bmin),
						  stop=np.log10(bmax),
						  num=steps_bt)
		#----------------------------------

		#------- Interpolate ---------------
		yb = np.interp(xb,xp=dfb[variate],
						  fp=dfb[covariate],
						  left=np.nan,
						  right=np.nan)
		#-----------------------------------

		plt.scatter(xb,yb,s=1,zorder=0,
				label="Interp BT-Settl")
		#======================================

		#------ Join data frames---------------------------------------
		tmp = pn.DataFrame(
			data={"x":np.concatenate((xb,dfa[variate].to_numpy())),
			      "y":np.concatenate((yb,dfa[covariate].to_numpy()))})
		tmp.sort_values(by="x",inplace=True)
		#-------------------------------------------------------------- 


		#------ Fit spline ------------------------------------------
		tck_v2c = splrep(tmp["x"],tmp["y"],
					xb=min_variate,xe=max_variate,k=degspl,s=stol_v2c)
		#-----------------------------------------------------------

		#------ Join interpolated results  --------
		dfo[covariate] = splev(dfo[variate],tck_v2c)
		#------------------------------------------

		#========== Inverted relation ===================================
		dfc = dfo[[variate,covariate]].copy()
		dfc.drop(dfc[(dfc["Mass"] > c2v_drop[0]) &
					 (dfc["Mass"] < c2v_drop[1])].index,
					 inplace=True)

		dfc.sort_values(by=covariate,inplace=True)


		tck_c2v = splrep(dfc[covariate],dfc[variate],k=degspl,s=stol_c2v)
		#-------------------------------------------------------------

		#------- Covariate grid ---------------------------------
		min_covariate = dfc[covariate].min()
		max_covariate = dfc[covariate].max()
		grid = np.linspace(min_covariate,max_covariate,num=steps)
		#--------------------------------------------------------

		# ----- Plot interpolated points ---------------------------
		plt.scatter(dfo[variate],dfo[covariate],c="black",
						s=0.1,zorder=1,label="Joint forward")

		plt.scatter(splev(grid,tck_c2v),grid,c="red",
						s=0.1,zorder=1,label="Joint backward")
		# -------------------------------------------------------

		#----- Save -----
		plt.xlabel(variate)
		plt.ylabel(covariate)
		plt.xscale("log")
		plt.legend(loc="best")
		pdf.savefig()
		plt.close()
		#----------------

		#------- Append tcks -----
		tcks_v2c.append(tck_v2c)
		tcks_c2v.append(tck_c2v)
		#-------------------------


	#=====================================================================

	#---------- Write output file --------------------------
	base_out = "{0}{1}_{2}_Myr_{3}".format(dir_out,base_,age,variate)
	dfo.to_csv(base_out + ".csv")
	# #-------------------------------------------------------

	#------ Save splines -----------
	tck = {}
	tck["mass2mag"] = tcks_v2c
	tck["mag2mass"] = tcks_c2v

	with open(base_out + ".pkl", 'wb') as out_strm: 
		dill.dump(tck, out_strm,protocol=2) 
#==================================================================================

#- Close pdf ---
pdf.close()
#----------