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
from __future__ import absolute_import,division,print_function
import sys
import myemcee as emcee
# import emcee
import numpy as np
import scipy.stats as st

class posterior_mass:
	"""
	This class provides flexibility to infer the posterior distribution of the mass
	"""
	def __init__(self,datum,uncert,loc_mu,scale_mu,N_bands,mass2phot,nwalkers=12,prior_mass="Uniform",min_mass=0.1,max_mass=10,
					burnin_frac=0.2):

		self.nwalkers    = nwalkers
		self.burnin_frac = burnin_frac
		self.ndim        = 6
		self.max_mass    = max_mass
		self.min_mass    = min_mass
		self.N_bands     = N_bands
		self.mass2phot   = mass2phot
		self.a           = 3.5


		#------------------------------------------------------
		idx     = np.ix_(np.where(np.isfinite(datum))[0])
		o_phot  = datum[idx]
		u_phot  = uncert[idx]

		self.idx = idx
		self.obs = o_phot
		self.unc = u_phot


		################################# PRIORS ######################################################

		def log_prior_mu_Pb_Yb_Vb(mu,Pb,Yb,Vb,V):
			'''
			Normal  prior for distance modulus
			Uniform prior for Pb between 0 and 1
			Normal  prior for Yb 
			Half-Cauchy  prior for Vb scale of 10 
			Half-Cauchy  prior for V scale of 10 
		    '''
			lp_mu = st.norm.logpdf(mu,loc=loc_mu,scale=scale_mu)
			lp_Pb = st.truncnorm.logpdf(Pb,a=0.0,b=1.0,loc=0,scale=1.0)
			lp_Yb = st.norm.logpdf(Yb,loc=np.mean(o_phot),scale=5.0*np.std(o_phot))
			lp_Vb = st.halfcauchy.logpdf(Vb,loc=1e-3,scale=10)
			lp_V  = st.halfcauchy.logpdf(V,loc=1e-6,scale=10)
			return lp_Pb + lp_Yb + lp_Vb + lp_V


		if prior_mass=="Uniform" :
			def lnprior(theta):
				uniform_mass = st.uniform.logpdf(theta[0],loc=min_mass,scale=max_mass-min_mass)
				uniform_aux  = log_prior_mu_Pb_Yb_Vb(theta[1],theta[2],theta[3],theta[4],theta[5])
				return(uniform_mass+uniform_aux)

		if prior_mass=="Half-Cauchy" :
			def lnprior(theta):
				pri_mass = st.halfcauchy.logpdf(theta[0],loc=0.0,scale=100.0)
				pri_aux  = log_prior_mu_Pb_Yb_Vb(theta[1],theta[2],theta[3],theta[4],theta[5])
				return(pri_mass+pri_aux)

		self.pos0 = [np.array([st.norm.rvs(loc=min_mass + 0.5*(max_mass-min_mass),scale=0.1,size=1)[0],
				st.norm.rvs(loc=loc_mu,scale=scale_mu,size=1)[0],
                st.uniform.rvs(loc=0,scale=0.1,size=1)[0],
                st.uniform.rvs(loc=0,scale=30,size=(1))[0],
                st.uniform.rvs(loc=1e-3,scale=10,size=(1))[0],
                st.uniform.rvs(loc=1e-6,scale=0.1,size=(1))[0]]) for i in range(self.nwalkers)]

		self.lnprior = lnprior

		############################################################################################

	############################# LIKELIHOOD ###############################################

	def log_likelihood(self,parameters):
	    '''
	    This log likelihood depends on the mass, photometry and uncertainty of the later
	    '''
	    mass = parameters[0] # The mass of the star
	    mu   = parameters[1] # Distance modulus
	    Pb   = parameters[2] # The probability of being an outlier
	    Yb   = parameters[3] # The mean position of the outlier distribution
	    Vb   = parameters[4] # The variance of the outlier distribution
	    V    = parameters[5] # The variance added to the photometry

	    true_phot = np.array([self.mass2phot[i](mass) for i in range(self.N_bands)]) + mu

	    t_phot    = true_phot[self.idx]

	    good = (1-Pb)*(1.0/np.sqrt(2.0*np.pi*(self.unc**2+V)))*np.exp(-0.5*((self.obs-t_phot)**2)/(2.0*(self.unc**2 + V)))
	    bad  = (Pb)*(1.0/np.sqrt(2.0*np.pi*(self.unc**2 + Vb)))*np.exp(-0.5*((self.obs-Yb)**2)/(2.0*(self.unc**2 + Vb))) +1e-200
	    
	    return(np.sum(np.log(good+bad)))


	################ POSTERIOR#######################
	def lnprob(self,theta):
		not_good_values = (theta[0] > self.max_mass or theta[0] < self.min_mass or 
						   theta[1] < 0.0 or theta[1] > 30 or
						   theta[2] < 0.0 or theta[2] > 1.0 or 
						   theta[4] < 1e-3  or theta[4]> 50 or
						   theta[5] < 1e-6  or theta[5]> 10) 
		if not_good_values:
			return(-np.inf)

		return(self.lnprior(theta) + self.log_likelihood(theta))

	#################### RUN THE SAMPLER ####################################

	def run(self,N_iter):
		sampler = emcee.EnsembleSampler(self.nwalkers,self.ndim, self.lnprob)
		sampler.a = self.a
		sampler.run_mcmc(self.pos0,N_iter)
		# print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))

		sample = sampler.chain[:,int(self.burnin_frac*N_iter):,:]

		#-----MAP ----
		ind_map = np.unravel_index(np.argmax(sampler.lnprobability),np.shape(sampler.lnprobability))
		MAP  = sampler.chain[ind_map[0],ind_map[1],:]
		#-----------
		Mean = np.mean(sample,axis=(0,1))
		#----- SD ------
		SD  = np.std(sample,axis=(0,1))
		#---- CI 95%
		CI  = np.percentile(sample,axis=(0,1),q=[2.5,97.5])
		#------ autocorrelation time
		int_time = emcee.autocorr.integrated_time(sample[:,:,0].flatten(),axis=0)#,c=1)

		return MAP,Mean,SD,CI,int_time,sample,np.mean(sampler.acceptance_fraction)
