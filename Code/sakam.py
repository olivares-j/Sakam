'''
Copyright 2018 Javier Olivares Romero

This file is part of Sakam.

    Sakam is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PyAspidistra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sakam.  If not, see <http://www.gnu.org/licenses/>.
'''
#------------ LOAD LIBRARIES -------------------
from __future__ import absolute_import,division,print_function
import sys
import os

import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec

import numpy as np
import scipy.stats as st
import random
import scipy
import corner
import h5py
import pandas as pn
import progressbar

from scipy.interpolate import interp1d
from matplotlib.ticker import NullFormatter
from abs2mass import posterior_mass

############### Files and directories #########################################
dir_main        = "/home/javier/Cumulos/Stock_2/"
#---------- Input files --------------------------
file_isochrone  = dir_main + "COLIBRI.csv"
file_data       = dir_main + "Absolute_magnitudes.csv"
#---------- Output files -------------------------------
dir_out         = dir_main + "Masses/" 
file_out        = dir_out  + "Masses.h5"
file_mass_csv   = dir_out  + "Masses.csv"
file_avs_csv    = dir_out  + "Avs.csv"

#------- Creates directories -------
if not os.path.exists(dir_out):
    os.makedirs(dir_out)
#---------------------------------
######################################################################

########################### VARIABLES COLIBRI ##########################
identifier   = ['source_id']
variate      = ['Mini']
covariates   = ['G_BPmag','Gmag','G_RPmag',
                'gP1mag','rP1mag','iP1mag','zP1mag','yP1mag',
                'Jmag','Hmag','Ksmag']
observable   = ['abs_BP','abs_G','abs_RP',
                'abs_gmag','abs_rmag','abs_imag','abs_zmag','abs_ymag',
                'abs_Jmag','abs_Hmag','abs_Kmag']
uncertainty  = ['abs_BP_error','abs_G_error','abs_RP_error',
                'abs_e_gmag','abs_e_rmag','abs_e_imag','abs_e_zmag','abs_e_ymag',
                'abs_e_Jmag','abs_e_Hmag','abs_e_Kmag']
#----------- The following value transforms the visual extinction to each of the bands --------
av2al        = [1.06794,0.85926,0.65199,1.16529,0.86813,0.67659,0.51743,0.43092,0.29434,0.18128,0.11838]
########################################################################################################################

# ############################ VARIABLES MIST ##########################
# identifier   = ['ID_member']
# variate      = ['initial_mass']
# covariates   = ['BP','G','RP','g_sdss','r_sdss','i_sdss','z_sdss','J','H','K']
# observable   = ['BP','G','RP','g_sdss','r_sdss','i_sdss','z_sdss','J','H','Ks']
# uncertainty  = ['e_BP','e_G','e_RP','e_g_sdss','e_r_sdss','e_i_sdss','e_z_sdss','e_J','e_H','e_Ks']
# #----------- The following value transforms the visual extinction to each of the bands --------
# av2al        = [1.067,0.859,0.6519,1.20585,0.87122,0.68319,0.49246,0.28887,0.18353,0.11509]
# # ########################################################################################################################



######################## Parameters ################
N_iter   = 2000
nwalkers = 30
npar     = 6
prior    = "Half-Cauchy" # Prior for the mass
upper_limit_variate = 3.5 # Max mass in the isochrone where it is still produces an injective function
name_parameters = [r"Mass $[\mathrm{M_{\odot}}]$",r"Av", r"$Pb$", r"$Yb$",r"$Vb$",r"$V$"]
quantiles = [0.16,0.84]
nan_threshold = 3 # Only objects with at least three bands and their uncertainties will be used
##########################################################


############################### ############################
cols_data    = sum([identifier,observable,uncertainty],[])
cols_isoc    = sum([variate,covariates],[])
#################################################################



############################ LOAD ISOCHRONE ########################################################
isochrone = pn.read_csv(file_isochrone,usecols=cols_isoc,dtype=np.float64)
isochrone = isochrone.drop_duplicates(subset=variate)
idx_valid = np.where(isochrone[variate]< upper_limit_variate)[0]
isochrone = isochrone.iloc[idx_valid,:]
####################################################################################################

############################ LOAD DATA ####################################################
data  = pn.read_csv(file_data,usecols=cols_data)
data  = data.reindex(columns=cols_data)

#---------------------- Select only some objects -------------
# idx = np.in1d(data[variate], [456753077401059072])
# data = data.iloc[[0,1],:]
#------------------------------------------------------------

#------- index as string ------
data[cols_data[0]] = data[cols_data[0]].astype('str')


#----- put ID as row name-----
data.set_index(cols_data[0],inplace=True)

data  = data.dropna(thresh=nan_threshold*2) 

N         = data.shape[0]
N_bands   = len(observable)
############################################################################################

########################### INTERPOLATING FUNCTION OF THE MASSS ###########################
'''
This is intended to give a smooth representation of the mass.
For a given value of the mass it returns a set of true photometric values.
'''

mass2phot = [interp1d(isochrone[variate[0]],isochrone[cov],kind="cubic") for cov in covariates]

min_variate  = np.min(isochrone[variate[0]])
max_variate  = np.max(isochrone[variate[0]])

print("The range of the variate is [{0},{1}].".format(min_variate,max_variate))
############################################################################################

###########################################################################################
###### LOOP OVER STARS TO INFER MASSS
############################################################################################

#--------- Use existing sources ------------------
if os.path.exists(file_out):
    fh5       = h5py.File(file_out,'r')
    ids       = fh5.keys()
    data.drop(ids,inplace=True)
    fh5.close()
    fh5       = h5py.File(file_out,'a')

else:
    fh5       = h5py.File(file_out,'w')

av2al     = np.array(av2al)
print("Sampling the posterior ...")

# ------- Prepare plots --------------------
nullfmt = plt.NullFormatter()
left, width = 0.1, 0.4
bottom, height = 0.1, 0.4
bottom_h = left_h = left + width + 0.0
rect_scatter = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.4]
rect_histy = [left_h, bottom, 0.1, height]

#----- start the progress bar ------
bar = progressbar.ProgressBar(maxval=N).start()

i = 0
for ID,datum in data.iterrows():
    print("ID: ",ID)
    grp = fh5.create_group(ID)
    #----------------
    observed  = datum[observable].values
    uncert    = datum[uncertainty].values

    #------ Initialize the module --------
    Module = posterior_mass(observed,uncert,
            N_bands=N_bands,mass2phot=mass2phot,
            av2al=av2al,
            nwalkers=nwalkers,
            prior_mass=prior,
            min_variate=min_variate,
            max_variate=max_variate,
            quantiles=quantiles,
            burnin_frac=0.25)

    #------- run the module -------------------------------------------------------------
    MAP,Median,SD,CI,int_time,sample,mean_acceptance_fraction = Module.run(N_iter=N_iter)

    #--------- Flatten sample ----------------------------------------------------------
    sample_flatten = sample.reshape((sample.shape[0]*sample.shape[1],sample.shape[2])).T

    #------- Save sample ------------------------------------
    dset = grp.create_dataset("MAP",    data=MAP)
    dset = grp.create_dataset("Median", data=Median)
    dset = grp.create_dataset("SD",     data=SD)
    dset = grp.create_dataset("CI",     data=CI)
    dset = grp.create_dataset("sample", data=sample_flatten)
    fh5.flush()

    #---- populate arrays--------------------------------------------
    print("Acceptance fraction: ",mean_acceptance_fraction)

    pdf = PdfPages(filename=dir_out+"Object_{0}.pdf".format(str(ID)))
    y_min,y_max= 0.95*np.min(sample[:,:,0]),1.05*np.max(sample[:,:,0])

    fig = plt.figure(221, figsize=(10, 10))
    ax0 = fig.add_subplot(223, position=rect_scatter)
    ax0.set_xlabel("Iteration")
    ax0.set_ylabel(r"Mass $[\mathrm{M_{\odot}}]$")
    ax0.set_ylim(y_min,y_max)
    ax0.plot(sample[:,:,0].T, '-', color='black', alpha=0.3,linewidth=0.3)
    ax0.axhline(MAP[0],  color='blue',ls="-",linewidth=0.5,label="MAP")
    ax0.axhline(CI[0,0], color='blue',ls=":",linewidth=0.5,label="CI")
    ax0.axhline(CI[1,0], color='blue',ls=":",linewidth=0.5)
    ax0.legend(loc="upper left",ncol=4,fontsize=4)


    ax1 = fig.add_subplot(224, position=rect_histy)
    ax1.set_ylim(y_min,y_max)
    ax1.axhline(MAP[0],     color='blue',ls="-",linewidth=0.5,label="MAP")
    ax1.axhline(CI[0,0],    color='blue',ls=":",linewidth=0.5,label="CI")
    ax1.axhline(CI[1,0],    color='blue',ls=":",linewidth=0.5)

    ax1.set_xlabel("Density")
    ax1.yaxis.set_ticks_position('none') 

    xticks = ax1.xaxis.get_major_ticks() 
    xticks[0].label1.set_visible(False)

    ax1.yaxis.set_major_formatter(nullfmt)
    ax1.yaxis.set_minor_formatter(nullfmt)

    ax1.hist(sample[:,:,0].flatten(),bins=100,density=True, 
        color="k",orientation='horizontal', fc='none', histtype='step',lw=0.5)
    pdf.savefig(bbox_inches='tight')  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    idx     = np.ix_(np.where(np.isfinite(datum))[0])
    true_phot = np.array([mass2phot[j](MAP[0]) for j in range(N_bands)]) + MAP[1]*av2al
    x  = np.arange(N_bands)

    plt.scatter(x,true_phot,color="grey",label="Model")
    plt.errorbar(x,observed,yerr=uncert,fmt=".",label="Observed")
    plt.xticks(x,observable,rotation='vertical')
    plt.margins(0.2)
    # Tweak spacing to prevent clipping of tick-labels
    plt.subplots_adjust(bottom=0.15)
    plt.legend(loc="upper right",fontsize=4)
    plt.ylabel("Magnitude")
    plt.xlabel("Filter")
    pdf.savefig(bbox_inches='tight')
    plt.close()

    # Corner plot
    fig = plt.figure()
    figure = corner.corner(sample_flatten.T, labels=name_parameters,truths=MAP,truth_color="red",
               quantiles=quantiles,
               show_titles=True, title_kwargs={"fontsize": 12})
    pdf.savefig(bbox_inches='tight')  # saves the current figure into a pdf page
    plt.close() 
    pdf.close()

    #----- update bar ----
    bar.update(i+1)
    i += 1
fh5.close()

#----------- Compute statistics -----------------------
fh5  = h5py.File(file_out,'r')
ids  = fh5.keys()

#------------ Intitialize arrays and directory ----------------
N = len(ids)
maps      = np.zeros((N,npar))
medians   = np.zeros((N,npar))
cis       = np.zeros((N,2,npar))

for i,key in enumerate(ids):
    grp = fh5.get(key)
    maps[i]    = np.array(grp.get("MAP"))
    medians[i] = np.array(grp.get("Median"))
    cis[i]     = np.array(grp.get("CI"))
fh5.close()


#---------- Masses ------------------------------------------------
masses = pn.DataFrame(np.column_stack((ids,maps[:,0],medians[:,0],
    cis[:,0,0],cis[:,1,0])),
        columns=[identifier[0],'map',
        'median',
        'lower','upper'])

avs = pn.DataFrame(np.column_stack((ids,maps[:,1],medians[:,1],
    cis[:,0,1],cis[:,1,1])),
        columns=[identifier[0],'map',
        'median',
        'lower','upper'])

M_tot = np.sum(maps[:,0])
M_low = M_tot - np.sum(cis[:,0,0])
M_up  = np.sum(cis[:,1,0]) - M_tot

print("Total mass: {0:3.1f}_{1:3.1f}^{2:3.1f}".format(M_tot,M_low,M_up))

masses.to_csv(path_or_buf=file_mass_csv,index=False)
avs.to_csv(path_or_buf=file_avs_csv,index=False)


