import dill
############ NOTE ##################################
#### DO NOT ADD VARIABLES NEITHER COMPOSE THEM #####
#### THE BASH FILE WILL NOT RECOGNIZE THEM     #####
####################################################

#----------- Split ------------------------------
# Number of processes/parts to partition the data
size = 4
#------------------------------------------------

#------ Ages --------------------------
ages = [1,3,5,7,10]
age  = "YYY"
#-------------------------------------

#------ Models --------------------------------------
# NOTE The following variable is read by the bash script
models = ["PMB"] 
#,"PARSEC","MIST","PMB"]
model = "XXX"
#----------------------------------------------------

#-------------- Directories ----------------------------------
dir_sakam  = "/home/javier/Repositories/Sakam"
dir_models = "/home/javier/Cumulos/Perseus/Models/"
dir_base   = "/home/javier/Cumulos/Perseus/Data/Groups/K8/"
dir_main   = dir_base + "/" + model + "/" + age + "/"
#-------------------------------------------------------------

#------- Global variables ------------------
name_globals   = "globals.pkl"
file_globals   = dir_main + name_globals
#--------------------------------------------

#---------- Input files --------------------------------------------------
file_isochrone  = dir_models + model + "_" + age + "Myr.csv"
file_photometry = dir_base   + "members+2MASS+PanSTARRS.csv"
file_distances  = dir_base   + "Kalkayotl/Sources_statistics.csv"
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
max_variate  = 10.0
# THe following are the covariates in the isochrone
# !!!!!! KEEP THE SAME ORDER AS IN OBSERVED PHOTOMETRY !!!!!!!!!!!!!
gaia_covs      = ['G_BPmag','Gmag','G_RPmag',]
twomass_covs   = ['Jmag','Hmag','Kmag']
panstarrs_covs = ['gP1mag','rP1mag','iP1mag','zP1mag','yP1mag']
# Effective wave lengths  
gaia_waves     = [5335.42, 6422.01, 7739.17]
twomass_waves  = [12329.79, 16395.59, 21522.05]
panstarrs_waves= [4907.71, 6208.38, 7531.06, 8669.70, 9619.44] #6395.40]
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
prior = {"variate":"LogNorm",
         "Av":"Uniform",
         "Rv":"Gaussian"
        }

hyper = {"loc_Av":0.0,
         "scl_Av":10.0,
         "loc_Rv":3.1,
         "scl_Rv":0.5,
         "beta_sd_b":1.0,
         "beta_sd_m":0.1}
#------------------------------------

#------- Running -------
iterations = 4000
walkers_ratio = 4
burnin_fraction = 0.5
# Hyper-parameters to generate initial solution
initial_hyper = {
		"loc_variate":0.2,
        "scl_variate":0.1,
        "loc_Av":0.05,
        "scl_Av":0.01,
        "loc_Rv":3.1,
        "scl_Rv":0.1,
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
