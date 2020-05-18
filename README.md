# Sakam
Infer masses from luminosities and isochrones

## Installation:
----------------
create a conda environement:
```conda create -n sakam python=3.8```

and activate it with ``conda activate sakam``
and install the following packages:
```
conda install -c conda-forge numpy
conda install -c conda-forge pandas
conda install -c conda-forge scipy
conda install -c conda-forge h5py
conda install -c conda-forge corner
conda install -c conda-forge progressbar
conda install -c conda-forge matplotlib
```
Install also [emcee](https://github.com/dfm/emcee) from the repository.
Note: The new version has an improved mixing of walkers.

## How to run:
--------------
1. Run [kalkayotl](https://github.com/olivares-j/Kalkayotl) to obtain individual distances to the stars in your sample.
2. Modify the code apparent2absolute.py to transform apparent magnitudes into absolute ones using the computed distances.
3. Download the isochrone (e.g. Parsec) with the photometric bands you intent to use (save it as a csv file)
4. Modify the example.py file to provide:

   1. Input files (absolute magnitudes and isochrone). 
   2. Variable names. 
   3. MCMC parameters. 
   4. Additional parameters. 
5. Activate the environment (``conda activate sakam``) and run the script with: ``python example.py``
