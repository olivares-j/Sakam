import dill

#----------- Split ------------------------------
# Number of processes/parts to partition the data
size = 1
#------------------------------------------------

#-------------- Directories ----------------------------------
dir_sakam       = "/home/javier/Repositories/Sakam"
dir_base        = "/home/javier/Cumulos/Perseus/Data/Groups/K7/"
dir_main        = "/home/javier/Cumulos/Perseus/Data/Groups/K7/Sakam/"
#-------------------------------------------------------------

#------- Global variables ------------------
name_globals   = "globals.pkl"
file_globals   = dir_main + name_globals
#--------------------------------------------

#---------- Input files --------------------------------------------------
file_isochrone  = "/home/javier/Cumulos/Perseus/Models/BT-Settl_5Myr.csv"
file_photometry = dir_base + "members+2MASS+PanSTARRS.csv"
file_distances  = dir_base + "Kalkayotl/Sources_statistics.csv"
#------------------------------------------------------------------------

#---------- Output files ------------------------------
# This will be created after running the unify.py routine
file_samples    = dir_main + "samples.h5"
file_statistics = dir_main + "statistics.csv"
#----------------------------------------------------------

#-------------------- Data -------------------------------------------------------------------
identifier   = 'source_id'
init_variate = None # Use the initial guess of the value. Otherwise None
kalkayotl_dimension = 6
#------------------------------------------------------------------------------------

#---------------- Observed photometry  ----------------------------------
gaia_bands       = ["bp","g","rp"]
gaia_errors      = ["bp_error","g_error","rp_error"]
twomass_bands    = ["Jmag","Hmag","Kmag"]
twomass_errors   = ["e_Jmag","e_Hmag","e_Kmag"]
panstarrs_bands  = ["gmag","rmag","imag","zmag","ymag"]
panstarrs_errors = ["e_gmag","e_rmag","e_imag","e_zmag","e_ymag"]
#-------------------------------------------------------------------------

#--------------------- Isochrone -----------------------------------------
variate      = "Mass"
max_variate  = 3.5
# THe following are the covariates in the isochrone
# !!!!!! KEEP THE SAME ORDER AS IN OBSERVED PHOTOMETRY !!!!!!!!!!!!!
gaia_covs      = ['G_BPmag','Gmag','G_RPmag',]
twomass_covs   = ['Jmag','Hmag','Kmag']
panstarrs_covs = ['gP1mag','rP1mag','iP1mag','zP1mag','yP1mag']
# The following values transform the visual extinction 
# to the extinction in each bands 
gaia_av2al     = [1.06794,0.85926,0.65199]
twomass_av2al  = [0.29434,0.18128,0.11838]
panstarrs_av2al= [1.16529,0.86813,0.67659,0.51743,0.43092]
#----------------------------------------------------------------------

#------ Preprocessing ------------------------------
nan_values = 99.0

# These are the filtering values of BP
label_BP = None # Use None to avoid the filtering
limit_BP = 15.0

n_obs_min = 3 # Minimum number of observed bands

add_unc = 0.05
#-------------------------------------------------


#-- Prior and hyper-parameter -----
prior_variate = "Chabrier"

hyper = {  "alpha_Av":1.0,
           "beta_Av":2.0,
           "beta_sd_b":1.0,
           "beta_sd_m":0.1}
#------------------------------------

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
plots_scale  = "log"
quantiles    = [0.16,0.84]
name_variate = r"Mass $[\mathrm{M_{\odot}}]$"
#-------------------------------------------

#------- Save globals -----------
dill.dump_session(file_globals)
