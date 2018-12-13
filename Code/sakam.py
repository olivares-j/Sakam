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
import numpy as np
import scipy.stats as st
import random
import scipy
import corner
import h5py
import pandas as pn
import progressbar
import matplotlib.pyplot as plt

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
from scipy.interpolate import interp1d
from astroML.decorators import pickle_results
from matplotlib.ticker import NullFormatter
from abs2mass import posterior_mass

############### DIRECTORIES #########################################
dir_data        = "/home/jromero/Desktop/Rup147/"
file_isochrone  = dir_data + "Models/COLIBRI_2.5Gyr.csv"
file_data       = dir_data + "Absolute_Magnitudes_unified.csv"
dir_out         = dir_data + "Masses/COLIBRI/" 
file_out        = dir_out  + "Masses.h5"
file_out_csv    = dir_out  + "Masses.csv"

final_masses    = False

if not os.path.isdir(dir_out):
    os.mkdir(dir_out)
######################################################################

############################ VARIABLES ##########################
identifier   = ['ID']
variate      = ['Mass']
covariates   = ['BP','G','RP','u_sdss','g_sdss','r_sdss','i_sdss','z_sdss','Y','J','H','K']
observable   = ['BP','G','RP','u_sdss','g_sdss','r_sdss','i_sdss','z_sdss','Y','J','H','K']
uncertainty  = ['e_BP','e_G','e_RP','e_u_sdss','e_g_sdss','e_r_sdss','e_i_sdss','e_z_sdss','e_Y','e_J','e_H','e_K']
########################################################################################################################

######################## Parameters ################
N_iter   = 3000
nwalkers = 50
npar     = 5
prior    = "Half-Cauchy" # Prior for the mass
quantiles = [16,84]
##########################################################


################################ RED PREPROCESS ############################

cols_data    = sum([identifier,observable,uncertainty],[])
cols_isoc    = sum([variate,covariates],[])

#----------- types ----------------------------------------
def types(col_name):
    if col_name == identifier[0]:
        return(np.dtype('S19'))
    else:
        return(np.float64)

types_cols   = dict([(i, types(i)) for i in cols_data])

#################################################################



############################ LOAD ISOCHRONE ########################################################
isochrone = pn.read_csv(file_isochrone,usecols=cols_isoc,dtype=np.float64)
isochrone = isochrone.drop_duplicates(subset=variate)
####################################################################################################

############################ LOAD DATA ####################################################
data  = pn.read_csv(file_data,usecols=cols_data,dtype=types_cols,na_values="99.0")#,nrows=10)
data  = data.reindex(columns=cols_data)
#------- index as string ------
data[cols_data[0]] = data[cols_data[0]].astype('str')

#----- put ID as row name-----
data.set_index(cols_data[0],inplace=True)

data  = data.dropna(thresh=3*2) # Only objects with at least three bands and their uncertainties

N     = data.shape[0]
N_bands   = len(observable)
############################################################################################

########################### INTERPOLATING FUNCTION OF THE MASSS ###########################
'''
This is intended to give a smooth representation of the mass.
For a given value of the mass it returns a set of true photometric values.
'''
mass2phot = [interp1d(isochrone[covariates[0]],isochrone[obs],kind="cubic") for obs in observable]

min_mass  = np.min(isochrone[covariates[0]])
max_mass  = np.max(isochrone[covariates[0]])

print("The range of masses is [{0},{1}].".format(min_mass,max_mass))
############################################################################################

###########################################################################################
###### LOOP OVER STARS TO INFER MASSS
############################################################################################
#------------ Intitialize arrays and directory ----------------
fh5       = h5py.File(file_out,'w') 
maps      = np.zeros((N,npar))
medians   = np.zeros((N,npar))
times     = np.zeros(N)
cis       = np.zeros((N,2,npar))
acc_fcs   = np.zeros(N)
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
    #----------------
    observed  = datum[observable].values
    uncert    = datum[uncertainty].values

    #------ Initialize the module --------
    Module = posterior_mass(observed,uncert,
            N_bands=N_bands,mass2phot=mass2phot,nwalkers=nwalkers,
            prior_mass=prior,min_mass=min_mass,max_mass=max_mass,
            quantiles=quantiles,
            burnin_frac=0.25)
    #------- run the module ----------------------------
    MAP,Median,SD,CI,int_time,sample,mean_acceptance_fraction = Module.run(N_iter=N_iter)

    #------- Save sample -----
    dset = fh5.create_dataset(str(ID), data=sample)
    fh5.flush()

    #---- populate arrays----
    maps[i,:]      = MAP
    medians[i,:]   = Median
    times[i]       = int_time
    cis[i,:,:]     = CI
    acc_fcs[i]     = mean_acceptance_fraction

    pdf = PdfPages(filename=dir_out+"Object_{0}.pdf".format(str(ID)))
    y_min,y_max= 0.95*np.min(sample[0]),1.05*np.max(sample[0])

    fig = plt.figure(221, figsize=(10, 10))
    ax0 = fig.add_subplot(223, position=rect_scatter)
    ax0.set_xlabel("Iteration")
    ax0.set_ylabel(r"Mass $[\mathrm{M_{\odot}}]$")
    ax0.set_ylim(y_min,y_max)
    ax0.plot(sample[0], '-', color='k', alpha=0.3,linewidth=0.3)
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

    ax1.hist(sample[0],bins=100,density=True, 
        color="k",orientation='horizontal', fc='none', histtype='step',lw=0.5)
    pdf.savefig(bbox_inches='tight')  # saves the current figure into a pdf page
    plt.close()

    plt.figure()
    idx     = np.ix_(np.where(np.isfinite(datum))[0])
    true_phot = np.array([mass2phot[j](MAP[0]) for j in range(N_bands)]) + MAP[1]
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
    figure = corner.corner(sample.T, labels=[r"Mass $[\mathrm{M_{\odot}}]$", r"$Pb$", r"$Yb$",r"$Vb$",r"$V$"],
               quantiles=quantiles,
               show_titles=True, title_kwargs={"fontsize": 12})
    pdf.savefig(bbox_inches='tight')  # saves the current figure into a pdf page
    plt.close() 
    pdf.close()

    
    if final_masses :
        #transform initial masses to final masses
        sample[0]     = init2final(sample[0])
        MAP[0]        = init2final(MAP[0])
        Mean[0]       = init2final(Mean[0])

    #----- update bar ----
    bar.update(i+1)
    i += 1
fh5.close()

print("Acceptance fraction statistics:")
print("Min.: {0:.3f}, Mean: {1:.3f}, Max.: {2:.3f}".format(np.min(acc_fcs),np.mean(acc_fcs),np.max(acc_fcs)))

print("Autocorrelation times statistics:")
print("Min.: {0:.3f}, Mean: {1:.3f}, Max.: {2:.3f}".format(np.min(times),np.mean(times),np.max(times)))

#---------- return data frame----
data_out = pn.DataFrame(np.column_stack((data.index,maps[:,0],medians[:,0],
    cis[:,0,0],cis[:,1,0],times)),
        columns=[identifier[0],'map_mass',
        'median_mass',
        'ci_low_distance','ci_up_distance',
        'integrated_autocorr_time'])

data_out.to_csv(path_or_buf=file_out_csv,index=False)


