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
dir_main        = "/home/javier/Cumulos/Stock_2/"
#---------- Input files --------------------------
file_isochrone  = dir_main + "COLIBRI.csv"
file_data       = dir_main + "Absolute_magnitudes.csv"
#---------- Output directory -------------------------------
dir_out         = dir_main + "Masses" 
######################################################################

########################### VARIABLES COLIBRI ##########################
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


sak = Sakam(file_data=file_data,
            file_isochrone=file_isochrone,
            dir_output=dir_out,
            identifier='source_id',
            variate='Mini',
            covariates=covariates,
            bands=observables,
            errors=uncertainty,
            av2al=av2al)

sak.run()
