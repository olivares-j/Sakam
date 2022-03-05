# Sakam
Infer masses from magnitudes and isochrones

## Installation:
----------------
create a conda environement:
```conda create -n sakam python=3.8```

and activate it with ``conda activate sakam``
and install the following packages:
```
conda install -c conda-forge numpy pandas scipy h5py corner progressbar matplotlib emcee dill extinction numba
```

Note: The [emcee](https://github.com/dfm/emcee) version from the repository shows slightly better convergence.

## How to run:
--------------
1. Run [kalkayotl](https://github.com/olivares-j/Kalkayotl) to obtain individual distances to the stars in your sample.
2. Download the isochrone (e.g. Parsec) with the photometric bands you intent to use (save it as a csv file)
4. Modify the ``globals_example.py`` file to provide:

   1. Input files. 
   2. Variable names. 
   3. MCMC parameters. 
   4. Additional parameters. 

5. Activate the environment (``conda activate sakam``) and run the script with: ``launch_sakam globals_example.py``

## Simple version
-------------------
A simple example on how to use Sakam is in ``example.py``. Use it if you already have the absolute magnitudes.
