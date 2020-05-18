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
import os
import numpy as np
import pandas as pd
import scipy.stats as st
from scipy.interpolate import interp1d
import h5py
import corner
import emcee
import progressbar


#---------------- Matplotlib ----------
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import gridspec
from matplotlib.ticker import NullFormatter


class posterior_mass:
    """
    This class provides flexibility to infer the posterior distribution of the mass given the absolute magnitudes
    It also infers the extinction nu
    """
    def __init__(self,observed,uncert,N_bands,mass2phot,av2al,hyper,walkers_ratio,prior_mass="Uniform",min_variate=0.1,max_variate=10,
                    burnin_frac=0.2,quantiles=[2.5,97.5]):

        
        self.burnin_frac = burnin_frac
        self.ndim        = 6
        self.nwalkers    = walkers_ratio * self.ndim
        self.max_variate = max_variate
        self.min_variate = min_variate
        self.N_bands     = N_bands
        self.mass2phot   = mass2phot
        self.av2al       = av2al
        self.quantiles   = quantiles


        #------------------------------------------------------
        idx     = np.ix_(np.where(np.isfinite(observed))[0])[0]
        o_phot  = observed[idx]
        u_phot  = uncert[idx]

        self.idx = idx
        self.obs = o_phot
        self.unc = u_phot
        ################################# PRIORS ######################################################

        prior_Av   = st.gamma(a=hyper["alpha_Av"],scale=hyper["beta_Av"])
        prior_Pb   = st.uniform(loc=0,scale=1.0)
        prior_Yb   = st.norm(loc=np.mean(o_phot),scale=5.0*np.std(o_phot))
        prior_sd_b = st.gamma(a=2.0,scale=hyper["beta_sd_b"])
        prior_sd_m = st.gamma(a=2.0,scale=hyper["beta_sd_m"])

        if prior_mass=="Uniform":
            prior_Ms = st.uniform(loc=min_variate,scale=max_variate-min_variate)
        elif prior_mass == "Half-Cauchy":
            prior_Ms = st.halfcauchy(loc=0.0,scale=100.0)
        elif prior_mass == "Chabrier":
            prior_Ms = st.lognorm(s=0.55,loc=np.log(0.2))
        else:
            sys.exit("Incorrect prior type")


        def lnprior(theta):
            lp_Ms   = prior_Ms.logpdf(theta[0])
            lp_Av   = prior_Av.logpdf(theta[1])
            lp_Pb   = prior_Pb.logpdf(theta[2])
            lp_Yb   = prior_Yb.logpdf(theta[3])
            lp_sd_b = prior_sd_b.logpdf(theta[4])
            lp_sd_m = prior_sd_m.logpdf(theta[5])
            return lp_Ms+lp_Av+lp_Pb+lp_Yb+lp_sd_b+lp_sd_m


        self.pos0 = [
                np.array([
                st.norm.rvs(loc=1.0,scale=0.1,size=1)[0],
                st.halfnorm.rvs(loc=0.05,scale=0.01,size=1)[0],
                st.halfnorm.rvs(loc=0.05,scale=0.01,size=1)[0],
                prior_Yb.rvs(size=1)[0],
                prior_sd_b.rvs(size=1)[0],
                prior_sd_m.rvs(size=1)[0]
                ]) 
                for i in range(self.nwalkers)]

        self.lnprior = lnprior

        ############################################################################################

    ############################# LIKELIHOOD ###############################################

    def log_likelihood(self,parameters):
        '''
        This log likelihood depends on the mass, photometry and uncertainty of the later
        '''
        mass = parameters[0] # The mass of the star
        Av   = parameters[1] # Extinction
        Pb   = parameters[2] # The probability of being an outlier
        Yb   = parameters[3] # The mean position of the outlier distribution
        sd_b = parameters[4] # The variance of the outlier distribution
        sd_m = parameters[5] # The variance added to the photometry

        true_phot = np.array([self.mass2phot[i](mass) for i in range(self.N_bands)])

        redden_phot = true_phot + Av*self.av2al

        t_phot    = redden_phot[self.idx]

        good = (1-Pb)*(1.0/np.sqrt(2.0*np.pi*(self.unc**2+sd_m**2)))*np.exp(-((self.obs-t_phot)**2)/(2.0*(self.unc**2 + sd_m**2)))
        bad  = (Pb)*(1.0/np.sqrt(2.0*np.pi*(self.unc**2 + sd_b**2)))*np.exp(-((self.obs-Yb)**2)/(2.0*(self.unc**2 + sd_b**2))) +1e-200
        
        return np.sum(np.log(good+bad))


    ################ POSTERIOR#######################
    def lnprob(self,theta):
        not_good_values = (theta[0] > self.max_variate or theta[0] < self.min_variate or 
                           theta[1] < 0.0 or
                           theta[2] < 0.0 or theta[2] >= 1.0 or 
                           theta[4] < 0.0 or
                           theta[5] < 0.0 ) 
        if not_good_values:
            return -np.inf

        return self.lnprior(theta) + self.log_likelihood(theta)

    #################### RUN THE SAMPLER ####################################

    def run(self,N_iter):
        sampler = emcee.EnsembleSampler(self.nwalkers,self.ndim, self.lnprob)
        sampler.run_mcmc(self.pos0,N_iter)
        sample = sampler.chain[:,int(self.burnin_frac*N_iter):,:]

        #-----MAP ----
        ind_map = np.unravel_index(np.argmax(sampler.lnprobability),np.shape(sampler.lnprobability))
        MAP  = sampler.chain[ind_map[0],ind_map[1],:]
        #-----------
        Median = np.median(sample,axis=(0,1))
        #----- SD ------
        SD  = np.std(sample,axis=(0,1))
        #---- CI 95% and median ---------
        CI  = np.quantile(sample,axis=(0,1),q=self.quantiles)

        return MAP,Median,SD,CI,sample,np.mean(sampler.acceptance_fraction)

class Sakam:
    def __init__(self,dir_output,hyperparameters=None,
                quantiles=[0.16,0.84]):

        self.name_parameters = [r"Mass $[\mathrm{M_{\odot}}]$",r"Av", r"$Pb$", r"$Yb$",r"$Sdb$",r"$Sdm$"]
        self.quantiles       = quantiles
        self.hyper           = hyperparameters
        self.n_parameters    = len(self.name_parameters)

        #-------- Hyper-parameters --------------
        if hyperparameters is None:
            self.hyper = {  "alpha_Av":1.0,
                            "beta_Av":2.0,
                            "beta_sd_b":1.0,
                            "beta_sd_m":0.1}

        #-------------- Files --------------------------
        os.makedirs(dir_output,exist_ok=True)
        self.file_h5    = dir_output + "/samples.h5"
        self.file_mass  = dir_output + "/"
        self.dir_plots  = dir_output + "/Plots"
        self.file_stats = dir_output + "/statistics.csv"
        os.makedirs(self.dir_plots,exist_ok=True)
        
        
        #------- Plot parameters ---------------
        self.nullfmt = plt.NullFormatter()
        left, width = 0.1, 0.4
        bottom, height = 0.1, 0.4
        bottom_h = left_h = left + width + 0.0
        self.rect_scatter = [left, bottom, width, height]
        self.rect_histx   = [left, bottom_h, width, 0.4]
        self.rect_histy   = [left_h, bottom, 0.1, height]


    def load_isochrone(self,file_isochrone,variate,covariates,av2al,upper_limit_variate=3.5):
        cols_isoc    = sum([[variate],covariates],[])
        
        print("Loading isochrone")
        isochrone = pd.read_csv(file_isochrone,usecols=cols_isoc,dtype=np.float64)
        print("Dropping duplicated values")
        isochrone = isochrone.drop_duplicates(subset=[variate])
        print("Dropping variate values larger than upper limit")
        idx_valid = np.where(isochrone[variate]< upper_limit_variate)[0]
        isochrone = isochrone.iloc[idx_valid,:]
        ####################################################################################################

        ########################### INTERPOLATING FUNCTION OF THE MASSS ###########################
        '''
        This is intended to give a smooth representation of the mass.
        For a given value of the mass it returns a set of true photometric values.
        '''

        self.mass2phot = [interp1d(isochrone[variate],isochrone[cov],kind="cubic") for cov in covariates]

        self.min_variate  = np.min(isochrone[variate])
        self.max_variate  = np.max(isochrone[variate])

        print("The range of the variate is [{0:2.2f},{1:2.2f}].".format(self.min_variate,self.max_variate))
        ############################################################################################

        self.av2al = np.array(av2al)


    def load_data(self,file_data,identifier,bands,errors,nan_threshold=3):
        columns_data = sum([[identifier],bands,errors],[])
        data         = pd.read_csv(file_data,usecols=columns_data)
        data         = data.reindex(columns=columns_data)

        #------- index as string ------
        data[identifier] = data[identifier].astype('str')

        #----- put ID as row name-----
        data.set_index(identifier,inplace=True)

        self.data  = data.dropna(thresh=nan_threshold*2) 

        self.N       = self.data.shape[0]
        self.n_bands = len(bands)
        self.bands   = bands
        self.errors  = errors
        self.identifier = identifier
        ############################################################################################


    def run(self,iterations=3000,walkers_ratio=4,burnin_fraction=0.3,prior_mass="Chabrier"):
        #--------- Use existing sources ------------------
        if os.path.exists(self.file_h5):
            print("Using existing samples")
            fh5       = h5py.File(self.file_h5,'r')
            ids       = fh5.keys()
            self.data.drop(ids,inplace=True)
            fh5.close()
            self.fh5  = h5py.File(self.file_h5,'a')

        else:
            self.fh5  = h5py.File(self.file_h5,'w')

        print("Sampling the posterior")

        bar = progressbar.ProgressBar(maxval=self.N).start()
        i   = 0

        for ID,datum in self.data.iterrows():
            grp = self.fh5.create_group(ID)

            #------ Initialize the module --------
            Module = posterior_mass(
                    datum[self.bands].values,
                    datum[self.errors].values,
                    N_bands=self.n_bands,
                    mass2phot=self.mass2phot,
                    av2al=self.av2al,
                    hyper=self.hyper,
                    walkers_ratio=walkers_ratio,
                    prior_mass=prior_mass,
                    min_variate=self.min_variate,
                    max_variate=self.max_variate,
                    quantiles=self.quantiles,
                    burnin_frac=burnin_fraction)

            #------- run the module -------------------------------------------------------------
            MAP,Median,SD,CI,sample,acceptance_fraction = Module.run(N_iter=iterations)

            #------- Save ------------------------------------
            dset = grp.create_dataset("MAP",    data=MAP)
            dset = grp.create_dataset("Median", data=Median)
            dset = grp.create_dataset("SD",     data=SD)
            dset = grp.create_dataset("CI",     data=CI)
            dset = grp.create_dataset("sample", data=sample)
            dset = grp.create_dataset("acc", data=acceptance_fraction)
            self.fh5.flush()

            #----- update bar ----
            bar.update(i+1)
            i += 1

        self.fh5.close()


    def plots(self,):
        print("Plotting posterior samples")
        fh5  = h5py.File(self.file_h5,'r')
        for ID,datum in self.data.iterrows():
            grp = fh5.get(ID)
            MAP     = np.array(grp.get("MAP"))
            CI      = np.array(grp.get("CI"))
            sample  = np.array(grp.get("sample"))
            acc     = np.array(grp.get("acc"))

            if acc < 0.2:
                print("Warning: source {0} has low acceptance fraction!".format(ID))

            self.plot_source(ID,
                    datum[self.bands].values,
                    datum[self.errors].values,
                    sample,MAP,CI)

        fh5.close()


    def statistics(self,):
        #----------- Compute statistics -----------------------
        fh5  = h5py.File(self.file_h5,'r')
        ids  = list(fh5.keys())

        #------------ Intitialize arrays and directory ----------------
        N = len(ids)
        maps      = np.zeros((N,self.n_parameters))
        medians   = np.zeros((N,self.n_parameters))
        cis       = np.zeros((N,2,self.n_parameters))

        for i,key in enumerate(ids):
            grp = fh5.get(key)
            maps[i]    = np.array(grp.get("MAP"))
            medians[i] = np.array(grp.get("Median"))
            cis[i]     = np.array(grp.get("CI"))
        fh5.close()

        M_tot = np.sum(maps[:,0])
        M_low = M_tot - np.sum(cis[:,0,0])
        M_up  = np.sum(cis[:,1,0]) - M_tot

        print("Total mass: {0:3.1f}_{1:3.1f}^{2:3.1f}".format(M_tot,M_low,M_up))

        #---------- output -----------
        data = {
                self.identifier:ids,
                "lower_mass": cis[:,0,0],
                "map_mass":    maps[:,0],
                "upper_mass": cis[:,1,0],
                "lower_av":   cis[:,0,1],
                "map_av":      maps[:,1],
                "upper_av":   cis[:,1,1]
                }
        #-------- Save data -------------------------------------
        out = pd.DataFrame(data)
        out.to_csv(path_or_buf=self.file_stats,index=False)

    def plot_source(self,ID,observed,uncert,sample,MAP,CI):
            file_plot = self.dir_plots+"/source_{0}.pdf".format(str(ID))

            pdf = PdfPages(filename=file_plot)
            y_min,y_max= 0.95*np.min(sample[:,:,0]),1.05*np.max(sample[:,:,0])

            fig = plt.figure(221, figsize=(10, 10))
            ax0 = fig.add_subplot(223, position=self.rect_scatter)
            ax0.set_xlabel("Iteration")
            ax0.set_ylabel(r"Mass $[\mathrm{M_{\odot}}]$")
            ax0.set_ylim(y_min,y_max)
            ax0.plot(sample[:,:,0].T, '-', color='black', alpha=0.3,linewidth=0.3)
            ax0.axhline(MAP[0],  color='blue',ls="-",linewidth=0.5,label="MAP")
            ax0.axhline(CI[0,0], color='blue',ls=":",linewidth=0.5,label="CI")
            ax0.axhline(CI[1,0], color='blue',ls=":",linewidth=0.5)
            ax0.legend(loc="upper left",ncol=4,fontsize=4)


            ax1 = fig.add_subplot(224, position=self.rect_histy)
            ax1.set_ylim(y_min,y_max)
            ax1.axhline(MAP[0],     color='blue',ls="-",linewidth=0.5,label="MAP")
            ax1.axhline(CI[0,0],    color='blue',ls=":",linewidth=0.5,label="CI")
            ax1.axhline(CI[1,0],    color='blue',ls=":",linewidth=0.5)

            ax1.set_xlabel("Density")
            ax1.yaxis.set_ticks_position('none') 

            xticks = ax1.xaxis.get_major_ticks() 
            xticks[0].label1.set_visible(False)

            ax1.yaxis.set_major_formatter(self.nullfmt)
            ax1.yaxis.set_minor_formatter(self.nullfmt)

            ax1.hist(sample[:,:,0].flatten(),bins=100,density=True, 
                color="k",orientation='horizontal', fc='none', histtype='step',lw=0.5)
            pdf.savefig(bbox_inches='tight')  # saves the current figure into a pdf page
            plt.close()

            plt.figure()
            true_phot = np.array([self.mass2phot[j](MAP[0]) for j in range(self.n_bands)]) + MAP[1]*self.av2al
            x  = np.arange(self.n_bands)

            # plt.scatter(x,true_phot-true_phot,yerr=,color="grey",label="Model")
            plt.errorbar(x,observed-true_phot,yerr=uncert,fmt=".",label="Observed")
            plt.xticks(x,self.bands,rotation='vertical')
            plt.margins(0.2)
            # Tweak spacing to prevent clipping of tick-labels
            plt.subplots_adjust(bottom=0.15)
            plt.legend(loc="upper right",fontsize=4)
            plt.ylabel("$\\Delta$ Magnitude (observed-modelled)")
            plt.xlabel("Filter")
            pdf.savefig(bbox_inches='tight')
            plt.close()


            # Corner plot
            sample_flatten = sample.reshape((sample.shape[0]*sample.shape[1],sample.shape[2]))
            fig = plt.figure()
            figure = corner.corner(sample_flatten, labels=self.name_parameters,truths=MAP,
                        truth_color="red",
                        quantiles=self.quantiles,
                        show_titles=True, 
                        title_kwargs={"fontsize": 12})
            pdf.savefig(bbox_inches='tight')
            plt.close() 
            pdf.close()


