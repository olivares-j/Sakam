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

from sakam import Sakam

############### Files and directories #########################################
dir_main        = "/home/javier/Cumulos/Perseus/Data/Groups/K8/PMB/"
dir_plots       = dir_main + "plots"
#---------- Input files --------------------------
file_isochrone  = "/home/javier/Cumulos/Perseus/Models/PMB_3Myr.csv"
file_data       = dir_main + "Absolute_magnitudes.csv"
file_samples    = dir_main + "samples.h5"
file_statistics = dir_main + "statistics.csv"

os.makedirs(dir_plots,exist_ok=True)
######################################################################

########################### ISOCHRONE ##########################
variate      = "Mass"
max_variate  = 10.0
covariates   = ['G_BPmag','Gmag','G_RPmag',
                'gP1mag','rP1mag','iP1mag','zP1mag','yP1mag',
                'Jmag','Hmag','Kmag']
#-- Wave-lengths in angstroms--------
waves = [5335.42, 6422.01, 7739.17, 4907.71, 6208.38, 7531.06, 8669.70, 9619.44,
                12329.79, 16395.59, 21522.05]
############################################################################

###################### DATA ##################################################
identifier   = 'source_id'
init_variate = None # In case you have an initial guess of the value. Optional column.
observables   = ['abs_bp','abs_g','abs_rp',
                'abs_gmag','abs_rmag','abs_imag','abs_zmag','abs_ymag',
                'abs_Jmag','abs_Hmag','abs_Kmag']
uncertainties  = ['abs_bp_error','abs_g_error','abs_rp_error',
                'abs_e_gmag','abs_e_rmag','abs_e_imag','abs_e_zmag','abs_e_ymag',
                'abs_e_Jmag','abs_e_Hmag','abs_e_Kmag']

########################################################################################################################



#-- Prior and hyper-parameter -----
prior = {"variate":"Chabrier",
         "Av":"Uniform",
         "Rv":"Gaussian"
        }

hyper = {"loc_Av":0.0,
         "scl_Av":10.0,
         "loc_Rv":3.1,
         "scl_Rv":0.5,
         "alpha_Pb":[1,19],
         "beta_sd_b":1.0,
         "beta_sd_m":0.1}
#------------------------------------

#------ Preprocessing ----------------------------
n_obs_min = 3 # Minimum number of observed bands
#-------------------------------------------------

#------- Running -------
iterations = 2000
walkers_ratio = 4
burnin_fraction = 0.5
# Hyper-parameters to generate initial solution
initial_hyper = {
        "loc_variate":0.1,
        "scl_variate":0.1,
        "loc_Av":0.0,
        "scl_Av":0.1,
        "loc_Rv":3.1,
        "scl_Rv":0.01,
        }
#-----------------------

#-------- Plots and statistics -------------
plots_scale  = "log"
quantiles    = [0.16,0.84]
name_variate = r"Mass $[\mathrm{M_{\odot}}]$"
#-------------------------------------------

####################### RUN ################################
sakam = Sakam(file_samples=file_samples,
                prior=prior,
                hyperparameters=hyper,
                initial_hyper=initial_hyper,
                quantiles=quantiles,
                name_variate=name_variate)

sakam.load_isochrone(file_isochrone=file_isochrone,
                        variate=variate,
                        covariates=covariates,
                        waves=waves,
                        upper_limit_variate=max_variate)

sakam.load_data(file_data=file_data,
                    identifier=identifier,
                    bands=observables,
                    errors=uncertainties,
                    nan_threshold=n_obs_min,
                    init_variate=init_variate,
                    nrows=1)

sakam.run(iterations=iterations,
            walkers_ratio=walkers_ratio,
            burnin_fraction=burnin_fraction)

sakam.plots(dir_plots=dir_plots,scale=plots_scale)

sakam.statistics(file_statistics=file_statistics)
