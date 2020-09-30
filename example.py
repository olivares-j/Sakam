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
dir_main        = "/home/javier/Cumulos/Taurus/Sakam"
# dir_main        = "/raid/jromero/OCs/Taurus/Sakam"
dir_plots       = dir_main + "plots"
#---------- Input files --------------------------
file_isochrone  = dir_main + "/data/colibri_1Myr.csv"
file_data       = dir_main + "/data/Absolute_magnitudes.csv"
file_samples    = dir_main + "/samples.h5"
file_statistics = dir_main + "/statistics.csv"

os.makedirs(dir_plots,exist_ok=True)
######################################################################

########################### VARIABLES COLIBRI ##########################
variate      = "Mini"
max_variate  = 3.5
covariates   = ['G_BPmag','Gmag','G_RPmag',
                'gP1mag','rP1mag','iP1mag','zP1mag','yP1mag',
                'Jmag','Hmag','Ksmag']
#-- The following value transforms the visual extinction to each of the bands --------
av2al        = [1.06794,0.85926,0.65199,1.16529,0.86813,
                0.67659,0.51743,0.43092,0.29434,0.18128,0.11838]
############################################################################

###################### DATA ##################################################
identifier   = 'source_id'
init_variate = 'init_mass' # In case you have an initial guess of the value. Optional column.
observables   = ['abs_bp','abs_g','abs_rp',
                'abs_gmag','abs_rmag','abs_imag','abs_zmag','abs_ymag',
                'abs_Jmag','abs_Hmag','abs_Kmag']
uncertainties  = ['abs_bp_error','abs_g_error','abs_rp_error',
                'abs_e_gmag','abs_e_rmag','abs_e_imag','abs_e_zmag','abs_e_ymag',
                'abs_e_Jmag','abs_e_Hmag','abs_e_Kmag']

########################################################################################################################

#-- Prior and hyper-parameter -----
prior_variate="Chabrier"

hyper = {  "alpha_Av":1.0,
           "beta_Av":2.0,
           "beta_sd_b":1.0,
           "beta_sd_m":0.1}
#------------------------------------


#------ Preprocessing ------------------------------
nan_values = 99.0

# These are the filtering values of BP
label_BP = "BP" # Use None to avoid the filtering
limit_BP = 15.0

n_obs_min = 3 # Minimum number of observed bands

add_unc = 0.05
#-------------------------------------------------

#------- Running -------
iterations = 4000
walkers_ratio = 4
burnin_fraction = 0.5
# Parameters of solution
initial_hyper = {
        "loc_variate":1.0,
        "scl_variate":0.1,
        "loc_Av":0.05,
        "scl_Av":0.01,
        "loc_Pb":0.05,
        "scl_Pb":0.01
        }
#-----------------------


#-------- Plots and statistics -------------
plots_scale = "log"
quantiles   = [0.16,0.84]
name_variate= r"Mass $[\mathrm{M_{\odot}}]$"
#-------------------------------------------


####################### RUN ################################
sakam = Sakam(file_samples=file_samples,
                hyperparameters=hyper,
                quantiles=quantiles,
                name_variate=name_variate)

sakam.load_isochrone(file_isochrone=file_isochrone,
                        variate=variate,
                        covariates=covariates,
                        av2al=av2al,
                        upper_limit_variate=max_variate)

sakam.load_data(file_data=file_data,
                    identifier=identifier,
                    bands=observables,
                    errors=uncertainties,
                    nan_threshold=n_obs_min,
                    init_variate=init_variate)

sakam.run(iterations=iterations,
            walkers_ratio=walkers_ratio,
            burnin_fraction=burnin_fraction,
            prior_variate=prior_variate,
            initial_hyper=initial_hyper)

sakam.plots(dir_plots=dir_plots,scale=plots_scale)

sakam.statistics(file_statistics=file_statistics)

