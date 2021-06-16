import sys
import dill
import numpy as np
import pandas as pn
from scipy.interpolate import splrep,splev

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#---------- Input files --------------------------------------------------
dir_   = "/home/javier/Cumulos/Perseus/Models/"
models = ["PARSEC","MIST","BT-Settl"]
ages   = ["3","5","7","10"]
#------------------------------------------------------------------------

#---------- Output file ------------------------------
base_out = "PMB"
steps = 5000
file_plot = dir_ + "join.pdf"
base_tck  = "PMB_splines_"
#----------------------------------------------------------

#--------------------- Isochrone -----------------------------------------
stol_v2c = 0.3
stol_c2v = 1e-3
degspl = int(3)
variate      = "Mass"
max_variate  = 10.0
min_variate  = 0.01
cut_weights  = 1.0
covariates   = [
            'G_BPmag','Gmag','G_RPmag', # Gaia
            'Jmag','Hmag','Kmag', # 2MASS
            'gP1mag','rP1mag','iP1mag','zP1mag','yP1mag' #PanSTARRS
            ]
#----------------------------------------------------------------------

#----- Weight function -----------------------------
def weights(model, mass):
    if model=="BT-Settl":
        w = np.where(mass < cut_weights, 1.0,0.0)
    else:
        w = np.where(mass > cut_weights, 1.0,0.0)
    return w
#--------------------------------------------------



#------ Open plot -------
pdf = PdfPages(file_plot)
#------------------------

#================ Loop over ages ================================================
for age in ages:

    #------------- Loop over models --------------------------
    dfm = []
    for model in models:
        #------- Read input file -----------------------------------------
        file = "{0}{1}_{2}Myr.csv".format(dir_,model,age)
        tmp = pn.read_csv(file, usecols=sum([[variate],covariates],[]))
        #------- Add model and weight ----------------------------------------
        tmp["Model"] = model
        tmp["Weight"] = tmp[["Model",variate]].apply(lambda x: weights(*x),axis=1)

        #-------- Drop masses outside limit ----------------------
        tmp.drop(tmp[tmp[variate]>max_variate].index,inplace=True)
        #---------------------------------------------------------

        dfm.append(tmp)

    #----- Join data frames and sort them by mass --
    dfm = pn.concat(dfm)
    dfm.sort_values(by=variate,inplace=True)
    #-----------------------------------------------

    #------ Group by models ------
    dfgs = dfm.groupby("Model")
    #-----------------------------

    #--------- Create output grid ------------------------
    df = pn.DataFrame(data={variate:np.logspace(start=np.log10(min_variate),
                        stop=np.log10(max_variate),num=steps)})
    #-----------------------------------------------------

    #========== Loop over magnitude =======================================
    tcks_v2c = []
    tcks_c2v = []
    for covariate in covariates:
        #----- Initialize plot ----
        plt.figure()
        plt.suptitle(age + "Myr")
        #--------------------------

         #----- Plot original models ----------------------------
        for g in dfgs.groups:
            dfg = dfgs.get_group(g)
            plt.scatter(dfg[variate],dfg[covariate],
                        s=1,zorder=0,label=g)
        #-----------------------------------------------------------------

        #------ Fit spline ------------------------------------------
        tck_v2c = splrep(dfm[variate],dfm[covariate],w=dfm["Weight"],
                    xb=min_variate,xe=max_variate,k=degspl,s=stol_v2c)
        #-----------------------------------------------------------

        #------ Join interpolated results  --------
        df[covariate] = splev(df[variate],tck_v2c)
        #------------------------------------------

        dfc = df[[variate,covariate]].copy()
        dfc.sort_values(by=covariate,inplace=True)

        tck_c2v = splrep(dfc[covariate],dfc[variate],k=degspl,s=stol_c2v)
        #-------------------------------------------------------------

        #------- Covariate grid ---------------------------------
        min_covariate = dfc[covariate].min()
        max_covariate = dfc[covariate].max()
        grid = np.linspace(min_covariate,max_covariate,num=steps)
        #--------------------------------------------------------

        # ----- Plot interpolated points ---------------------------
        plt.scatter(df[variate],df[covariate],c="black",
                        s=0.1,zorder=1,label="Joint forward")

        plt.scatter(splev(grid,tck_c2v),grid,c="red",
                        s=0.1,zorder=1,label="Joint backward")
        # -------------------------------------------------------

        #----- Save -----
        plt.xlabel("Mass")
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
    file_out = "{0}{1}_{2}Myr.csv".format(dir_,base_out,age)
    df.to_csv(file_out)
    #-------------------------------------------------------

    #------ Save splines -----------
    tck = {}
    tck["mass2mag"] = tcks_v2c
    tck["mag2mass"] = tcks_c2v

    file_tck = dir_ + base_tck + age + "Myr.pkl"
    with open(file_tck, 'wb') as out_strm: 
        dill.dump(tck, out_strm,protocol=2) 
#==================================================================================

#- Close pdf ---
pdf.close()
#----------