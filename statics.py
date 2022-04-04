import sys
import os
import dill

#----- Global variables ----------
path_globals  = str(sys.argv[1]) 
dill.load_session(path_globals)
#--------------------------------

prefix     = "abs_"
dir_plots  = dir_main + "plots/"
dir_chunks = dir_main + "chunks/"
file_abs   = dir_main + "Absolute_magnitudes.csv"

bands = []
errors = []
covariates = []
waves = []
for instrument in instruments:
  bands.append(instrument["observables"])
  errors.append(instrument["uncertainties"])
  covariates.append(instrument["covariates"])
  waves.append(instrument["waves"])


bands  = sum(bands,[])
errors = sum(errors,[])
covariates = sum(covariates,[])
waves = sum(waves,[])

observables   = [ prefix + band for band in bands]
uncertainties = [ prefix + error for error in errors]

if init_variate is not None:
  columns_abs = sum([[identifier],observables,uncertainties,[init_variate]],[])
else:
  columns_abs = sum([[identifier],observables,uncertainties],[])


columns_photometry = sum([[identifier],bands,errors],[])

if kalkayotl_dimension == 1:
  # Kalkayotl first version
  columns_distances = sum([[identifier],["mean","lower","upper"]],[])
else:
  # Kalkayotl second version
  lowers = ["hdi_2.5%_"  + var for var in ["X","Y","Z"]]
  uppers = ["hdi_97.5%_" + var for var in ["X","Y","Z"]]
  columns_distances = sum([[identifier],["mode_distance"],lowers,uppers],[])

#-------- Create directories ------------
os.makedirs(dir_plots,exist_ok=True)
os.makedirs(dir_chunks,exist_ok=True)
#---------------------------------------

#------- Save globals -----------
dill.dump_session(file_globals)
