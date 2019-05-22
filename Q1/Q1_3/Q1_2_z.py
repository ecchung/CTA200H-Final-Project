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

# Define parameters
H0    = 67.5    
omch2 = 0.122  
ombh2 = 0.022  

# Setting up plot figure
plt.clf()
plt.figure(0, figsize=(10,6))

# Making a list of colours to iterate over
colors = np.array(['k', 'b', 'r', 'g','orange','c', 'm', 'y'])

# Iterating over this indicated list of parameters

# Set initial paramters
pars  = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2) 
pars.InitPower.set_params(ns=0.965)

# Note non-linear corrections
pars.set_matter_power(redshifts=[0., 0.8, 1.6, 2.4, 3.0, 5.0, 9.0], kmax=2.0) # default redshift z = [0.0]

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
for i in range(len(z)): 
    plt.loglog(kh, pk[i],color=colors[i], ls='-', label="z={0}, linear".format(z[i]))
    plt.loglog(kh_nonlin, pk_nonlin[i],color=colors[i], ls='--', label="z={0}, nonlinear".format(z[i]))
    
plt.xlabel('k/h Mpc');
plt.legend(loc='lower center')
plt.grid(True)
plt.title('Matter Power $P_m (k)$ at $H_0$={0}\n$\Omega_b h ^2$={1}, $\Omega_c h ^2$={2}'.format(H0, ombh2, omch2))    
plt.savefig("powerspectrum_var_z.pdf")                                 
    

