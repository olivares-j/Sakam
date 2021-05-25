#------------ LOAD LIBRARIES -------------------
from __future__ import absolute_import,division,print_function
import sys
import dill

#----- Global variables ----------
path_globals  = str(sys.argv[1]) 
dill.load_session(path_globals)
#---------------------------------

#------ Sakam ------------
sys.path.append(dir_sakam)
from sakam import Sakam
#--------------------------

#------------ Files -------------------------------------------------
file_chunk = dir_chunks + "/data_{0}_of_{1}.csv".format(XXX,size)
file_samp  = dir_chunks + "/samples_{0}_of_{1}.h5".format(XXX,size)
file_stat  = dir_chunks + "/statistics_{0}_of_{1}.csv".format(XXX,size)
#-------------------------------------------------------------------

#======================================================================
sakam = Sakam(file_samples=file_samp,
				hyperparameters=hyper,
				quantiles=quantiles,
				name_variate=name_variate)

sakam.load_isochrone(file_isochrone=file_isochrone,
                        variate=variate,
                        covariates=covariates,
                        av2al=av2al,
                        upper_limit_variate=max_variate)

sakam.load_data(file_data=file_chunk,
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

sakam.statistics(file_statistics=file_stat)
