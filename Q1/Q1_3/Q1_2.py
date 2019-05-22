# CTA200H Final Project 
# last edited May 17, 2019 1:20am
# Q1_2 varying ombh2, omch2, mnu

import sys, platform, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower

# Indicating which parameter to interate over
H0     = 67.5    # standard value

#ombh2 = np.array([0.022, 0.001, 0.01, 0.033]) # max 0.033 
omch2 = 0.122  # standard value

ombh2  = 0.022   # default
#omch2  = np.array([0.122, 0.4, 0.05]) # use when varying omch2: also try 0.03, 0.01
#mnu    = 1.6 # 0.1 ~ 0.3

#mnu = np.array([0.1, 0.01, 0.8, 1.8])
#tau   = 0.06

par_itr = mnu #set as the parameter that we want to iterate over

# Setting up plot figure
plt.clf()
plt.figure(0, figsize=(10,6))

# Making a list of colours to iterate over
colors = np.array(['k', 'b', 'r', 'g','orange','c', 'm', 'y'])

# Iterating over this indicated list of parameters
for i in range(len(par_itr)):
    # Set initial paramters
    print(i)
    pars  = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2, mnu=par_itr[i])    # CHANGE HERE
    pars.InitPower.set_params(ns=0.965)

    # Note non-linear corrections
    pars.set_matter_power(kmax=2.0) # default redshift z = [0.0]
    
    # Linear spectra
    pars.NonLinear = model.NonLinear_none
    results        = camb.get_results(pars)
    kh, z, pk      = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
    s8             = np.array(results.get_sigma8())

    # Non-Linear spectra (Halofit)
    pars.NonLinear = model.NonLinear_both
    results.calc_power_spectra(pars)
    kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
    
    # Plot matter power spectra at z = 0
    
    # For varying baryonic matter density
    #plt.loglog(kh, pk[0], color=colors[i], ls='-', label="$\Omega_m h ^2$ = {0}, linear".format(par_itr[i]))
    #plt.loglog(kh_nonlin, pk_nonlin[0], color=colors[i], ls=':', label="$\Omega_m h ^2$ = {0}, nonlinear".format(par_itr[i]))

    # For varying CDM matter density with fixed mnu
    #plt.loglog(kh, pk[0], color=colors[i], ls='-', label="$\Omega_c h ^2$ = {0}, linear".format(par_itr[i]))
    #plt.loglog(kh_nonlin, pk_nonlin[0], color=colors[i], ls=':', label="$\Omega_c h ^2$ = {0}, nonlinear".format(par_itr[i]))
    
    # For varying mnu with fixed standard CDM density 
    plt.loglog(kh, pk[0], color=colors[i], ls='-', label=r"$M_\nu$={0}, linear".format(par_itr[i]))
    plt.loglog(kh_nonlin, pk_nonlin[0], color=colors[i], ls=':', label=r"$M_\nu$={0}, nonlinear".format(par_itr[i]))
    
    
plt.xlabel('k/h Mpc');
plt.legend(loc='lower left', title="$\Omega_c h ^2$ = {0}".format(omch2))
plt.grid(True)

# For varying baryonic matter density
#plt.title('Matter Power $P_m (k)$ at z={0}\n$H_0$={1}, $\Omega_c h ^2$={2}'.format(z[0], H0, omch2))    
#plt.savefig("powerspectrum_var_baryonic_density.pdf".format(H0))                                        

# For varying CDM matter density
plt.title('Matter power at z={0}\n$H_0$={1}, $\Omega_b h ^2$={2}'.format(z[0], H0, ombh2)) 
plt.savefig("powerspectrum_var_cdm_mnu_{0}.pdf".format(mnu))                                 
    

