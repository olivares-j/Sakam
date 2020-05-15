# Sakam
Infer masses from luminosities and isochrones

Steps:
1. Run kalkayotl to obtain individual distances to the stars in your sample
2. Modify the code apparent2absolute.py to transform apparent magnitudes into absolute ones using the computed distances
3. Download the isochrone with the photometric bands you intent to use (save it as a csv file)
4. Modify the sakam.py file to provide:
	a. Input files (absolute magnitudes and isochrone)
	b. Variable names
	c. MCMC parameters
	d. Additional parameters
5. Run the script
