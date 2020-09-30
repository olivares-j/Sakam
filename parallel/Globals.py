#----------- Split ------------------------------
# Number of processes/parts to partition the data
size = 4
#------------------------------------------------
############### Files and directories #########################################
dir_sakam       = "/home/javier/Repositories/Sakam"
dir_main        = "/home/javier/Cumulos/Taurus/Sakam"
# dir_main        = "/raid/jromero/OCs/Taurus/Sakam"
dir_data        = dir_main + "/data"
dir_plots       = dir_main + "/plots"
dir_chunks      = dir_main + "/chunks"
#---------- Input files --------------------------
file_isochrone  = dir_data + "/colibri_1Myr.csv"
file_data       = dir_data + "/Absolute_magnitudes.csv"
#---------- Output files ------------------------------
# This will be created after running the unify.py routine
file_samples    = dir_main + "/samples_all.h5"
file_statistics = dir_main + "/statistics_all.csv"
######################################################################

########################### VARIABLES COLIBRI ##########################
variate      = "Mini"
covariates   = ['G_BPmag','Gmag','G_RPmag',
                'gP1mag','rP1mag','iP1mag','zP1mag','yP1mag',
                'Jmag','Hmag','Ksmag']

# max_variate indicates the maximum value allowed to the variate.
# Useful if you need to cut the sequence.              
max_variate  = 3.5

#----------- The following value transforms the visual extinction to each of the bands --------
av2al        = [1.06794,0.85926,0.65199,1.16529,0.86813,
				0.67659,0.51743,0.43092,0.29434,0.18128,0.11838]
###########################################################################################

################### Data ########################################
identifier   = 'source_id'

observables   = ['abs_bp','abs_g','abs_rp',
                'abs_gmag','abs_rmag','abs_imag','abs_zmag','abs_ymag',
                'abs_Jmag','abs_Hmag','abs_Kmag']
uncertainties  = ['abs_bp_error','abs_g_error','abs_rp_error',
                'abs_e_gmag','abs_e_rmag','abs_e_imag','abs_e_zmag','abs_e_ymag',
                'abs_e_Jmag','abs_e_Hmag','abs_e_Kmag']

# In case you have an initial guess of the value. Optional column.
init_variate = 'init_mass' 

# Remove init_variate if you do not have it in your data
columns_data = sum([[identifier],observables,uncertainties,[init_variate]],[])
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
label_BP = None # Use None to avoid the filtering
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
plots_scale  = "log"
quantiles    = [0.16,0.84]
name_variate = r"Mass $[\mathrm{M_{\odot}}]$"
#-------------------------------------------


