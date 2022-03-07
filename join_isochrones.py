import sys
import dill
import numpy as np
import pandas as pn
from scipy.interpolate import splrep,splev

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#---------- Input files --------------------------------------------------
dir_    = "/home/jolivares/Cumulos/Models/"
dir_out = "/home/jolivares/Cumulos/ComaBer/Sakam/models/"

variate = "Mass" #"Luminosity"
ages = [
	# {"value":200,"max_variate":3.81},
	# {"value":300,"max_variate":3.27},
	# {"value":400,"max_variate":2.94},
	# {"value":600,"max_variate":2.48},
	# {"value":700,"max_variate":2.35},
	# {"value":800,"max_variate":2.23},
	# {"value":900,"max_variate":2.14},
	{"value":1000,"max_variate":2.06}
	]

models  = [
		{
		"name":"PARSEC",
		"file":dir_ + "PARSEC/Parsec_600-1000_Myr_Gaia+PanSTARRS+2MASS.csv",
		"columns":["age_Myr","Mini","logL",
					"Gmag","G_BPmag","G_RPmag", #Gaia
					"gP1mag","rP1mag","iP1mag","zP1mag","yP1mag", #PanSTARRS
					"Jmag","Hmag","Ksmag" #2MASS
					],
		"mapper":{"Mini":"Mass", "logL":"Luminosity",
				"Gmag":"Gmag","G_BPmag":"G_BPmag","G_RPmag":"G_RPmag", # Gaia
				"gP1mag":"gP1mag","rP1mag":"rP1mag","iP1mag":"iP1mag","zP1mag":"zP1mag","yP1mag":"yP1mag", #PanSTARRS
				"Jmag":"Jmag","Hmag":"Hmag","Ksmag":"Kmag" # 2MASS
				},
		"lower_variate":0.13,"upper_variate":10.0},
		# {
		# "name":"MIST",
		# "file":dir_ + "MIST/MIST_{0}_Myr_Gaia+2MASS+PanSTARRS+WISE.csv".format(age),
		# "columns":["age_Myr","initial_mass","log_L",
		# 			"Gaia_G_EDR3","Gaia_BP_EDR3","Gaia_RP_EDR3", #Gaia
		# 			"PS_g","PS_r","PS_i","PS_z","PS_y", #PanSTARRS
		# 			"2MASS_J","2MASS_H","2MASS_Ks" #2MASS
		# 			],
		# "mapper":{"initial_mass":"Mass", "log_L":"Luminosity",
		# 		"Gaia_G_EDR3":"Gmag","Gaia_BP_EDR3":"G_BPmag","Gaia_RP_EDR3":"G_RPmag", # Gaia
		# 		"PS_g":"gP1mag","PS_r":"rP1mag","PS_i":"iP1mag","PS_z":"zP1mag","PS_y":"yP1mag", #PanSTARRS
		# 		"2MASS_J":"Jmag","2MASS_H":"Hmag","2MASS_Ks":"Kmag" # 2MASS
		# 		},
		# "lower_mass":1.4,"upper_mass":2.3},
		# {
		# "name":"BT-Settl",
		# "file":dir_ + "BT-Settl/BT-Settl_all_Myr_Gaia+2MASS+PanSTARRS.csv",
		# "columns":["age_Myr","M/Ms","L/Ls",
		# 			"G","G_BP","G_RP", #Gaia
		# 			"g_p1","r_p1","i_p1","z_p1","y_p1", #PanSTARRS
		# 			"J","H","K" #2MASS
		# 			],
		# "mapper":{"M/Ms":"Mass","L/Ls":"Luminosity",
		# 		"G":"Gmag","G_BP":"G_BPmag","G_RP":"G_RPmag", # Gaia
		# 		"g_p1":"gP1mag","r_p1":"rP1mag","i_p1":"iP1mag","z_p1":"zP1mag","y_p1":"yP1mag", #PanSTARRS
		# 		"J":"Jmag","H":"Hmag","K":"Kmag" # 2MASS
		# 		},
		# "lower_variate":0.0,"upper_variate":0.14}
		]

cmds = [
	{"color":["Gmag","G_RPmag"],"magnitude":"Gmag"},
	{"color":["G_BPmag","G_RPmag"],"magnitude":"Gmag"},
	{"color":["Gmag","Kmag"],"magnitude":"Gmag"}
	]
#------------------------------------------------------------------------

#---------- Output file ------------------------------
base_ = "P"
steps_out = 1000
steps_int = 1000
file_v2c = dir_out + base_+"_{0}_Myr_relations.pdf"
file_cmd = dir_out + base_+"_{0}_Myr_cmds.pdf"
#----------------------------------------------------------

#--------------------- Isochrone -----------------------------------------
stol_v2c = 5e-2
stol_c2v = 1e-3
degspl = int(3)
c2v_drop = [2.3,2.4]
covariates = [
			"Gmag","G_BPmag","G_RPmag", # Gaia
			"gP1mag","rP1mag","iP1mag","zP1mag","yP1mag", #PanSTARRS
			"Jmag","Hmag","Kmag" # 2MASS
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
	tmp.drop(tmp[~np.isin(tmp["age_Myr"],
				[a["value"] for a in ages])].index,
				inplace=True)
	#--------------------------------------------------------------

	#-------- Drop masses outside limit ---------------------------------
	tmp.drop(tmp[(tmp[variate] > model["upper_variate"]) | 
				 (tmp[variate] < model["lower_variate"])].index,inplace=True)
	#--------------------------------------------------------------------

	#------- Add model and weight ----------------
	tmp["Model"] = model["name"]
	#------------------------------------

	dfm.append(tmp)

#----- Join data frames and sort them by mass --
dfm = pn.concat(dfm,sort=False)
dfm.sort_values(by=variate,inplace=True)
#-----------------------------------------------

#================ Loop over ages ================================================
#------ Group by ages ------------------
dfga = dfm.groupby(by="age_Myr")
#-----------------------------------------

for age in ages:
	#------ Open plot -------
	pdf_v2c = PdfPages(file_v2c.format(age["value"]))
	pdf_cmd = PdfPages(file_cmd.format(age["value"]))
	#------------------------

	#--- Get age data --------
	dfa = dfga.get_group(age["value"])
	#-------------------------

	#-------- Drop masses outside limit ---------------------------------
	dfa.drop(dfa[dfa[variate] > age["max_variate"]].index,inplace=True)
	#--------------------------------------------------------------------

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
							stop=np.log10(max_variate),num=steps_out)})
	else:
		dfo = pn.DataFrame(data={variate:np.linspace(start=min_variate,
							stop=max_variate,num=steps_out)})
	#-----------------------------------------------------

	#========== Loop over magnitude =======================================
	tcks_v2c = []
	tcks_c2v = []
	for covariate in covariates:
		#----- Initialize plot ----
		fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
		plt.subplots_adjust(hspace=0.0)
		fig.suptitle(str(age["value"]) + " Myr",y=0.95)
		#--------------------------

		#========= Loop over models ==================
		tmp = []
		for g in dfgm.groups:
			dfg = dfgm.get_group(g)

			#---------- Min Max --------------
			gmax = dfg[variate].max()
			gmin = dfg[variate].min()
			#---------------------------------

			#--------- Grid ------------------
			xg = np.logspace(start=np.log10(gmin),
							  stop=np.log10(gmax),
							  num=steps_int)
			#----------------------------------

			#------- Interpolate ---------------
			yg = np.interp(xg,xp=dfg[variate],
							  fp=dfg[covariate],
							  left=np.nan,
							  right=np.nan)
			#-----------------------------------

			#------------ Plot ---------------------
			ax1.scatter(dfg[variate],dfg[covariate],
						s=3,zorder=0,marker="s",
						label="{0} original".format(g))
			ax1.scatter(xg,yg,
						s=0.1,zorder=0,marker=".",
						label="{0} interpolated".format(g))
			#---------------------------------------
			
			#---------------------------------------------
			tmp.append(pn.DataFrame(data={"x":xg,"y":yg}))
			#---------------------------------------------

			#===========================================================

		tmp = pn.concat(tmp,ignore_index=True)
		tmp.drop_duplicates(subset="x",inplace=True)
		tmp.dropna(inplace=True)
		tmp.sort_values(by="x",inplace=True)
		#--------------------------------------------------------------

		#------ Fit spline -------------------------------------------
		tck_v2c = splrep(tmp["x"],tmp["y"],
					xb=min_variate,xe=max_variate,k=degspl,s=stol_v2c)
		#-------------------------------------------------------------

		#------ Join interpolated results  --------
		dfo[covariate] = splev(dfo[variate],tck_v2c)
		#------------------------------------------

		# ----- Plot Spline and its derivative ----------------
		ax1.plot(dfo[variate],dfo[covariate],c="black",
						lw=0.5,zorder=1,label="Spline fit")
		ax2.plot(dfo[variate],splev(dfo[variate],tck_v2c,der=1),
						c="black",lw=0.5,zorder=1)
		#--------------------------------------------------------


		#========== Inverted relation ===================================
		dfc = dfo[[variate,covariate]].copy()
		dfc.drop(dfc[(dfc[variate] > c2v_drop[0]) &
					 (dfc[variate] < c2v_drop[1])].index,
					 inplace=True)

		dfc.sort_values(by=covariate,inplace=True)

		tck_c2v = splrep(dfc[covariate],dfc[variate],k=degspl,s=stol_c2v)
		#-------------------------------------------------------------

		#------- Covariate grid ---------------------------------
		min_covariate = dfc[covariate].min()
		max_covariate = dfc[covariate].max()
		grid = np.linspace(min_covariate,max_covariate,num=steps_out)
		#--------------------------------------------------------

		# ----- Plot Inverse spline and its derivative -----------
		# ax1.plot(splev(grid,tck_c2v),grid,c="red",
		# 				lw=0.5,zorder=1,label="Inverse")

		# ax2.plot(splev(grid,tck_c2v),grid,c="red",
		# 				lw=0.5,zorder=1,label="Inverse derivative")
		# -------------------------------------------------------

		#----- Save -----
		ax2.set_xlabel(variate)
		ax2.set_ylabel("Derivative of " + covariate)
		ax1.set_ylabel(covariate)
		ax1.set_xscale("log")
		ax1.legend(loc="lower right")
		ax1.invert_yaxis()
		ax1.grid(which="both",linewidth=0.5)
		ax2.grid(which="both",linewidth=0.5)
		pdf_v2c.savefig()
		plt.close()
		#----------------

		#------- Append tcks -----
		tcks_v2c.append(tck_v2c)
		tcks_c2v.append(tck_c2v)
		#-------------------------


	#=====================================================================

	#---------- Write output file --------------------------
	base_out = "{0}{1}_{2}_Myr_{3}".format(dir_out,base_,age["value"],variate)
	dfo.to_csv(base_out + ".csv")
	#-------------------------------------------------------

	for cmd in cmds:
		plt.figure()
		ax = plt.gca()
		ax.plot(
			dfo[cmd["color"][0]]-dfo[cmd["color"][1]],
			dfo[cmd["magnitude"]],
			linewidth=1,color="black",
			label="Interpolation")

		#----- Plot original models -----------------
		for g in dfgm.groups:
			dfg = dfgm.get_group(g)
			ax.scatter(
				x=dfg[cmd["color"][0]]-dfg[cmd["color"][1]],
				y=dfg[cmd["magnitude"]],
				s=3,zorder=0,label=g)
		#-------------------------------------------

		ax.set_ylabel(cmd["magnitude"])
		ax.set_xlabel("{0} - {1}".format(
					cmd["color"][0],
					cmd["color"][1]))
		ax.invert_yaxis()
		ax.grid(which="both",linewidth=0.5)
		ax.legend(loc="upper right")
		pdf_cmd.savefig()
		plt.close()

	#- Close pdf ---
	pdf_v2c.close()
	pdf_cmd.close()
	#----------

	#------ Save splines -----------
	tck = {}
	tck["mass2mag"] = tcks_v2c
	tck["mag2mass"] = tcks_c2v

	with open(base_out + ".pkl", 'wb') as out_strm: 
		dill.dump(tck, out_strm,protocol=2) 
	#==================================================================================

	