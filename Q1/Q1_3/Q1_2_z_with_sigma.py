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
from scipy.integrate import simps, trapz

# ================================================================
# Defining functions that calculate sigma R

def W(k, R):
    return 3 * (np.sin(k*R) / (k*R) - np.cos(k*R)) / (k*R)**3 
 
def D2(k, pmk):
    return k**3 * pmk / (2 * np.pi**2)
    
def integrand(R, k, pmk):
    return W(k,R)**2 * D2(k,pmk) / k

def sR2(R, k, pmk):
    return simps(integrand(R, k, pmk), k)

def sR(R, k, pmk):
    return np.sqrt(sR2(R, k, pmk))
    
# ================================================================

# Define parameters
H0    = 67.5    
omch2 = 0.122  
ombh2 = 0.022  

# Making a list of colours to iterate over
colors = np.array(['k', 'b', 'r', 'g','orange','c', 'm', 'y'])

# Iterating over this indicated list of parameters

# Set initial paramters
pars  = camb.CAMBparams()
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2) 
pars.InitPower.set_params(ns=0.965)

# Note non-linear corrections
pars.set_matter_power(redshifts=[0., 0.8, 1.6, 5.0, 15.0, 30.], kmax=2.0) # default redshift z = [0.0]

# Linear spectra
pars.NonLinear = model.NonLinear_none
results        = camb.get_results(pars)
kh, z, pk      = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
s8             = sR(8, kh, pk)
s200           = sR(200, kh, pk)
s8_camb        = np.array(results.get_sigma8())
print("s8 me: {0}\ns8 camb: {1}\ndifference: {2}".format(s8, s8_camb, abs(s8-s8_camb)))

# Non-Linear spectra (Halofit)
pars.NonLinear = model.NonLinear_both
results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
s8_n                           = sR(8, kh_nonlin, pk_nonlin)
s200_n                         = sR(200, kh_nonlin, pk_nonlin)
    
# Plot matter power spectra at z = 0
plt.clf()
plt.figure(figsize=(10,6))
for i in range(len(z)): 
    plt.loglog(kh, pk[i],color=colors[i], ls='-', label="$z$ = {0}   \tlinear  \t   $\sigma_8$={1} \t  $\sigma200$={2}".format(z[i], round(s8[i], 4), round(s200[i],4)))
    plt.loglog(kh_nonlin, pk_nonlin[i],color=colors[i], ls='--', label="$z$ = {0}    \tnonlinear \t$\sigma_8$={1} \t$\sigma200$={2}".format(z[i], round(s8_n[i], 4), round(s200_n[i],4)))
    
plt.xlabel('k/h Mpc');
plt.legend(loc='lower left')
plt.grid(True)
plt.title('Matter Power $P_m (k)$ at $H_0$={0}\n$\Omega_b h ^2$={1}, $\Omega_c h ^2$={2}'.format(H0, ombh2, omch2))    
plt.savefig("powerspectrum_var_z_with_sigma.pdf")                                 
    

