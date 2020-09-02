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

from sakam import Sakam

############### Files and directories #########################################
dir_main        = "/home/javier/Cumulos/Taurus/Sakam/"
# dir_main        = "/raid/jromero/OCs/Taurus/Sakam/"
#---------- Input files --------------------------
file_isochrone  = dir_main + "colibri_1Myr.csv"
file_data       = dir_main + "Absolute_magnitudes.csv"
######################################################################

########################### VARIABLES COLIBRI ##########################
identifier   = 'source_id'
init_mass    = 'init_mass' # In case you have an initial guess of the value. Optional column.
variate      = "Mini"
covariates   = ['G_BPmag','Gmag','G_RPmag',
                'gP1mag','rP1mag','iP1mag','zP1mag','yP1mag',
                'Jmag','Hmag','Ksmag']
observables   = ['abs_bp','abs_g','abs_rp',
                'abs_gmag','abs_rmag','abs_imag','abs_zmag','abs_ymag',
                'abs_Jmag','abs_Hmag','abs_Kmag']
uncertainties  = ['abs_bp_error','abs_g_error','abs_rp_error',
                'abs_e_gmag','abs_e_rmag','abs_e_imag','abs_e_zmag','abs_e_ymag',
                'abs_e_Jmag','abs_e_Hmag','abs_e_Kmag']
#----------- The following value transforms the visual extinction to each of the bands --------
av2al        = [1.06794,0.85926,0.65199,1.16529,0.86813,0.67659,0.51743,0.43092,0.29434,0.18128,0.11838]
########################################################################################################################
#----- Hyperparameters ----------------------
hyper = {  "alpha_Av":2.0,
            "beta_Av":1.5,
            "beta_sd_b":1.0,
            "beta_sd_m":0.1}
#-----------------------------

sakam = Sakam(dir_output=dir_main)

sakam.load_isochrone(file_isochrone=file_isochrone,
                        variate=variate,
                        covariates=covariates,
                        av2al=av2al)
sakam.load_data(file_data=file_data,
                    identifier=identifier,
                    bands=observables,
                    errors=uncertainties,
                    init_mass=init_mass)

sakam.run()
sakam.plots(scale="log") # Log ("log)") or linear ("lin") scale
sakam.statistics()

