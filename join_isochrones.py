import sys
import os
import dill
import numpy as np
import pandas as pn
from scipy.interpolate import splrep,splev

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pn.set_option('display.max_rows', None)

#---------- Directories -------------------------------
dir_mdl = "/home/jolivares/OCs/Models/"
dir_out = "/home/jolivares/OCs/GroupX/Sakam/iter_6/models/"
#------------------------------------------------------

variate = "Mass"

if variate == "Mass":
	# ---------------------- Mass ----------------------------------
	out_models = {
	"PARSEC":  {"join":["PARSEC"],     "lower_variate":[0.0],"upper_variate":[10.0]},
	"BT-Settl":{"join":["BT-Settl"],   "lower_variate":[0.0],"upper_variate":[10.0]},
	"PB":      {"join":["PARSEC","BT-Settl"],"lower_variate":[1.0,0.0],"upper_variate":[10.0,1.0]}
	}
	ages = [
		{"value":200,"max_variate":3.8,
		"stol_v2c":{"PARSEC":5e-2,"BT-Settl":1e-3,"PB":1e-2},
		"stol_c2v":{"PARSEC":1e-4,"BT-Settl":1e-4,"PB":1e-4}
		},
		{"value":300,"max_variate":3.2,
		"stol_v2c":{"PARSEC":5e-2,"BT-Settl":1e-3,"PB":1e-3},
		"stol_c2v":{"PARSEC":1e-4,"BT-Settl":1e-4,"PB":1e-4}
		},
		{"value":400,"max_variate":2.9,
		"stol_v2c":{"PARSEC":5e-2,"BT-Settl":1e-3,"PB":1e-2},
		"stol_c2v":{"PARSEC":1e-4,"BT-Settl":1e-4,"PB":1e-4}
		}
		]
	# --------------------------------------------------------------
elif variate == "Luminosity":
	#----------------- Luminosity ----------------------------------
	out_models = {
	"PARSEC":  {"join":["PARSEC"],     "lower_variate":[-6],"upper_variate":[10.0]},
	"BT-Settl":{"join":["BT-Settl"],   "lower_variate":[-10],"upper_variate":[10.0]},
	"PB":      {"join":["PARSEC","BT-Settl"],"lower_variate":[0,-10.0],"upper_variate":[10.0,0.0]}
	}
	ages = [
		{"value":200,"max_variate":10**(0.38),
		"stol_v2c":{"PARSEC":1e-2,"BT-Settl":1e-2,"PB":1e-2},
		"stol_c2v":{"PARSEC":1e-3,"BT-Settl":1e-3,"PB":1e-3}
		},
		{"value":300,"max_variate":10**(0.34),
		"stol_v2c":{"PARSEC":1e-2,"BT-Settl":5e-2,"PB":5e-2},
		"stol_c2v":{"PARSEC":1e-3,"BT-Settl":1e-3,"PB":5e-3}
		},
		{"value":400,"max_variate":10**(0.30),
		"stol_v2c":{"PARSEC":1e-2,"BT-Settl":5e-2,"PB":5e-2},
		"stol_c2v":{"PARSEC":1e-3,"BT-Settl":1e-3,"PB":5e-3}
		}
		]
	#------------------------------------------------------------
else:
	sys.exit("Error in variate name")

models  = {
		"PARSEC":{
		"file":dir_mdl + "PARSEC/Parsec_200-400_Myr_Gaia+PanSTARRS+2MASS.csv",
		"columns":["age_Myr","Mini","logL",
					"Gmag","G_BPmag","G_RPmag", #Gaia
					"gP1mag","rP1mag","iP1mag","zP1mag","yP1mag", #PanSTARRS
					"Jmag","Hmag","Ksmag" #2MASS
					],
		"mapper":{"age_Myr":"Age/Myr","Mini":"Mass", "logL":"Luminosity",
				"Gmag":"Gmag","G_BPftmag":"G_BPmag","G_RPmag":"G_RPmag", # Gaia
				"gP1mag":"gP1mag","rP1mag":"rP1mag","iP1mag":"iP1mag","zP1mag":"zP1mag","yP1mag":"yP1mag", #PanSTARRS
				"Jmag":"Jmag","Hmag":"Hmag","Ksmag":"Kmag" # 2MASS
				},
		"steps_interpolation":2000,
		},
		# {
		# "name":"MIST",
		# "file":dir_mdl + "MIST/MIST_1-10_Myr_Gaia+2MASS+PanSTARRS.csv",
		# "columns":["age_Myr","initial_mass","log_L",
		# 			"Gaia_G_EDR3","Gaia_BP_EDR3","Gaia_RP_EDR3", #Gaia
		# 			"PS_g","PS_r","PS_i","PS_z","PS_y", #PanSTARRS
		# 			"2MASS_J","2MASS_H","2MASS_Ks" #2MASS
		# 			],
		# "mapper":{"age_Myr":"Age/Myr","initial_mass":"Mass", "log_L":"Luminosity",
		# 		"Gaia_G_EDR3":"Gmag","Gaia_BP_EDR3":"G_BPmag","Gaia_RP_EDR3":"G_RPmag", # Gaia
		# 		"PS_g":"gP1mag","PS_r":"rP1mag","PS_i":"iP1mag","PS_z":"zP1mag","PS_y":"yP1mag", #PanSTARRS
		# 		"2MASS_J":"Jmag","2MASS_H":"Hmag","2MASS_Ks":"Kmag" # 2MASS
		# 		},
		# "steps_interpolation":10,
		# # "lower_variate":0.0,"upper_variate":4.0,
		# "lower_variate":1.0,"upper_variate":4.0
		# },
		"BT-Settl":{
		"file":dir_mdl + "BT-Settl/BT-Settl_all_Myr_Gaia+2MASS+PanSTARRS.csv",
		"columns":["age_Myr","M/Ms","L/Ls",
					"G","G_BP","G_RP", #Gaia
					"g_p1","r_p1","i_p1","z_p1","y_p1", #PanSTARRS
					"J","H","K" #2MASS
					],
		"mapper":{"age_Myr":"Age/Myr","M/Ms":"Mass","L/Ls":"Luminosity",
				"G":"Gmag","G_BP":"G_BPmag","G_RP":"G_RPmag", # Gaia
				"g_p1":"gP1mag","r_p1":"rP1mag","i_p1":"iP1mag","z_p1":"zP1mag","y_p1":"yP1mag", #PanSTARRS
				"J":"Jmag","H":"Hmag","K":"Kmag" # 2MASS
				},
		"steps_interpolation":100,
		},

		# {
		# "name":"BT-Cond",
		# "file":dir_mdl + "BT-Cond/2MASS+Gaia+Subaru+UKIDSS.csv",
		# "columns":["Age/Myr","M/Ms","L/Ls",
		# 			"G","G_BP","G_RP", #Gaia
		# 			"g","r","i","z","Y", #Subaru
		# 			"J","H","K", #2MASS
		# 			"y" #UKIDSS
		# 			],
		# "mapper":{"M/Ms":"Mass","L/Ls":"Luminosity",
		# 		"G":"Gmag","G_BP":"G_BPmag","G_RP":"G_RPmag", # Gaia
		# 		"g":"gmag","r":"rmag","i":"imag","z":"zmag","Y":"Ymag", #Subaru
		# 		"y":"Ymag_ukidss", #UKIDSS
		# 		"J":"Jmag","H":"Hmag","K":"Kmag" # 2MASS
		# 		},
		# "lower_variate":0.0,"upper_variate":10.0}
		# # "lower_variate":0.0,"upper_variate":1.4}
		}

cmds = [
	{"color":["G_BPmag","Gmag"],"magnitude":"Gmag"},
	# {"color":["G_BPmag","G_RPmag"],"magnitude":"Gmag"},
	{"color":["Gmag","Kmag"],"magnitude":"Gmag"}
	]
#------------------------------------------------------------------------

#---------- Output file ------------------------------
steps_out = 1000
base_rel  = dir_out + "{0}_{1}_Myr_{2}_relations.pdf"
base_cmd  = dir_out + "{0}_{1}_Myr_{2}_cmds.pdf"
base_csv  = dir_out + "{0}_{1}_Myr_{2}.csv"
base_pkl  = dir_out + "{0}_{1}_Myr_{2}.pkl"
#----------------------------------------------------------

#--------------------- Isochrone -----------------------------------------
degspl = int(3)
c2v_drop = [2.3,2.4]
covariates = [
			"Gmag","G_BPmag","G_RPmag", # Gaia
			"gP1mag","rP1mag","iP1mag","zP1mag","yP1mag", #PanSTARRS
			# "gmag","rmag","imag","zmag","Ymag", #PanSTARRS
			# "Ymag_ukidss", #UKIDSS
			"Jmag","Hmag","Kmag" # 2MASS
			 ]
#----------------------------------------------------------------------

#------------ Load isochrones --------------------
dfs_mod = []
for name,model in models.items():
	#------- Read input file ---------------------------------
	tmp = pn.read_csv(model["file"], usecols=model["columns"])
	#---------------------------------------------------------

	#-------------- Rename -----------------------------------
	tmp.rename(columns=model["mapper"],inplace=True)
	#---------------------------------------------------------

	#------- Ages to integers ------------------
	tmp["Age/Myr"] = tmp["Age/Myr"].astype(int)
	#-------------------------------------------

	#-------- Drop ages ----------------------------
	tmp.drop(index=tmp.loc[~np.isin(tmp["Age/Myr"],
				[a["value"] for a in ages])].index,
				inplace=True)
	#-----------------------------------------------

	#------- Add model and weight ----------------
	tmp["Model"] = name
	#------------------------------------

	dfs_mod.append(tmp)

#----- Join data frames and sort them by mass --
df_mod = pn.concat(dfs_mod,sort=False)
df_mod.sort_values(by=variate,inplace=True)
#-----------------------------------------------


os.makedirs(dir_out,exist_ok=True)
for key,values in out_models.items():
	print(10*"=" + " " + key + " " + 10*"=")

	#--------- Select models ----------------------------------------
	dfm = df_mod.loc[np.isin(df_mod["Model"],values["join"])].copy()
	dfm.reset_index(inplace=True)
	#----------------------------------------------------------------
	
	#-------- Drop variate outside limit ---------------------------------
	index = []
	for c,mod in enumerate(values["join"]):
		idx = dfm.loc[(
				(dfm["Model"] == mod ) & 
				(dfm[variate] < values["upper_variate"][c]) &
				(dfm[variate] > values["lower_variate"][c])
				)].index.values
		index.append(idx)
	index = np.concatenate(index)
	dfm = dfm.loc[index]
	#--------------------------------------------------------------------

	#================ Loop over ages ================================================
	#------ Group by ages ------------------
	dfga = dfm.groupby(by="Age/Myr")
	#-----------------------------------------
	for age in ages:
		#------ Open plot -------
		pdf_v2c = PdfPages(base_rel.format(key,age["value"],variate))
		pdf_cmd = PdfPages(base_cmd.format(key,age["value"],variate))
		#------------------------

		#--- Get age data -----------------
		dfa = dfga.get_group(age["value"])
		#----------------------------------

		#-------- Drop masses outside limit -------------
		dfa = dfa.loc[dfa[variate] < age["max_variate"]]
		#------------------------------------------------

		#---------- Min Max --------------
		max_variate  = dfa[variate].max()
		min_variate  = dfa[variate].min()
		#---------------------------------

		#------ Group by models ------------------
		dfgm = dfa.groupby(by="Model")
		#-----------------------------------------

		#--------- Create output grid ------------------------
		if variate == "Mass":
			dfo = pn.DataFrame(data={variate:np.logspace(
								start=np.log10(min_variate),
								stop=np.log10(max_variate),
								num=steps_out)})
		else:
			dfo = pn.DataFrame(data={variate:np.linspace(
								start=min_variate,
								stop=max_variate,
								num=steps_out)})
		#-----------------------------------------------------

		#========== Loop over magnitude =======================================
		tcks_v2c = {}
		tcks_c2v = {}
		for covariate in covariates:
			#----- Initialize plot ----
			fig, (ax1, ax2,ax3) = plt.subplots(nrows=3, sharex=False,figsize=(10,15))
			# plt.subplots_adjust(hspace=0.0)
			fig.suptitle(str(age["value"]) + " Myr",y=0.95)
			#--------------------------

			#========= Loop over models ==================
			tmp = []
			for name in values["join"]:
				dfg = dfgm.get_group(name)

				#---------- Min Max --------------
				gmax = dfg[variate].max()
				gmin = dfg[variate].min()
				#---------------------------------

				#--------- Grid ------------------
				if variate == "Mass":
					xg = np.logspace(start=np.log10(gmin),
								  stop=np.log10(gmax),
						num=models[name]["steps_interpolation"])
					xg_plt = np.log10(xg)
					xv_plt = np.log10(dfg[variate])
				elif variate == "Luminosity":
					xg = np.linspace(start=gmin,
								  stop=gmax,
						num=models[name]["steps_interpolation"])
					xg_plt = xg
					xv_plt = dfg[variate]
				else:
					sys.exit()
				#----------------------------------

				#------- Interpolate ---------------
				yg = np.interp(xg,xp=dfg[variate],
								  fp=dfg[covariate],
								  left=np.nan,
								  right=np.nan)
				#-----------------------------------

				#------------ Plot ---------------------
				ax1.scatter(xv_plt,dfg[covariate],
							s=3,zorder=0,marker="s",
							label="{0} original".format(name))
				ax1.scatter(xg_plt,yg,
							s=0.1,zorder=0,marker=".",
							label="{0} interpolated".format(name))
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
						xb=min_variate,xe=max_variate,k=degspl,
						s=age["stol_v2c"][key])
			#-------------------------------------------------------------

			#------ Join interpolated results  --------
			dfo[covariate] = splev(dfo[variate],tck_v2c)
			#------------------------------------------

			#========== Inverted relation ===================================
			tmp.drop(tmp[(tmp["x"] > c2v_drop[0]) &
						 (tmp["x"] < c2v_drop[1])].index,
						 inplace=True)

			tmp.sort_values(by="y",inplace=True)
			# print(tmp.describe())

			tck_c2v = splrep(tmp["y"],tmp["x"],k=degspl,
						s=age["stol_c2v"][key])
			#-------------------------------------------------------------

			#------- Append tcks --------
			tcks_v2c[covariate] = tck_v2c
			tcks_c2v[covariate] = tck_c2v
			#----------------------------

			#------- Covariate grid ---------------------------------
			min_covariate = dfo[covariate].min()
			max_covariate = dfo[covariate].max()
			grid = np.linspace(min_covariate,max_covariate,num=steps_out)
			#--------------------------------------------------------

			if variate == "Mass":
				xo_plt = np.log10(dfo[variate])
				yg_plt = np.log10(splev(grid,tck_c2v))

			elif variate == "Luminosity":
				xo_plt = dfo[variate]
				yg_plt = splev(grid,tck_c2v)
			else:
				sys.exit()

			# ----- Plot Spline and its derivative ----------------
			ax1.plot(xo_plt,dfo[covariate],c="black",
							lw=0.5,zorder=1,label="Spline fit")
			ax2.plot(xo_plt,splev(dfo[variate],tck_v2c,der=1),
							c="black",lw=0.5,zorder=1)
			#--------------------------------------------------------

			# ----- Plot Inverse spline and its derivative -----------
			ax1.plot(yg_plt,grid,c="red",
							lw=0.5,zorder=1,label="Inverse")

			ax3.plot(grid,splev(grid,tck_c2v,der=1),c="red",
							lw=0.5,zorder=1)
			# -------------------------------------------------------

			#----- Save -----
			ax1.set_xlabel(variate)
			ax1.set_ylabel(covariate)
			ax1.legend(loc="upper left")
			ax1.invert_yaxis()
			ax1.grid(which="both",linewidth=0.5)

			ax2.set_xlabel(variate)
			ax2.set_ylabel("Derivative of " + covariate)
			ax2.grid(which="both",linewidth=0.5)

			ax3.set_xlabel(covariate)
			ax3.set_ylabel("Derivative of " + variate)
			ax3.grid(which="both",linewidth=0.5)


			pdf_v2c.savefig(bbox_inches='tight')
			plt.close()
			#----------------

		#=====================================================================

		#---------- Write output file --------------------------
		dfo.to_csv(base_csv.format(key,age["value"],variate))
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

		with open(base_csv.format(key,age["value"],
			variate).replace("csv","pkl"), 'wb') as out_strm: 
			dill.dump(tck, out_strm,protocol=2)
	print(30*"=")
#==================================================================================

	