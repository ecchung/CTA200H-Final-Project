# CTA200H Final Project 
# last edited May 18, 2019 11:50pm
# Q1_2 varying ombh2, omch2, mnu with sigma

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

# Indicating which parameter to interate over
H0     = 67.5    

ombh2 = np.array([0.01, 0.022, 0.033]) # max 0.033 
omch2 = 0.122  # standard value
mnu    = 0.1 # 0.1 ~ 0.3

#ombh2  = 0.022  
#omch2  = np.array([0.122, 0.4, 0.05]) # use when varying omch2: also try 0.03, 0.01
#mnu    = 0.1 

#ombh2  = 0.022   
#omch2 = 0.122  
#mnu = np.array([0.1, 0.8, 1.8])

par_itr = ombh2 # set as the parameter that we want to iterate over         CHANGE HERE 

# Setting up figures
plt.clf()
plt.figure(figsize=(10,6))
# Making a list of colours to iterate over
colors = np.array(['k', 'b', 'r', 'g','orange','c', 'm', 'y'])

# Iterating over this indicated list of parameters
for i in range(len(par_itr)):
    # Set initial paramters
    print(i)
    pars  = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=par_itr[i], omch2=omch2, mnu=mnu) # add mnu=par_itr[i] for mnu varying     # CHANGE HERE
    pars.InitPower.set_params(ns=0.965)

    # Note non-linear corrections
    pars.set_matter_power(kmax=2.0) # default redshift z = [0.0]
    
    # Linear spectra
    pars.NonLinear = model.NonLinear_none
    results        = camb.get_results(pars)
    kh, z, pk      = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
    s8             = sR(8, kh, pk)[0]
    s200           = sR(200, kh, pk)[0]
    s8_camb        = np.array(results.get_sigma8())[0]
    print("s8 me: {0}, s8 camb: {1}, difference: {2}".format(s8, s8_camb, abs(s8-s8_camb)))

    # Non-Linear spectra (Halofit)
    pars.NonLinear                 = model.NonLinear_both
    results.calc_power_spectra(pars)
    kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
    s8_n                           = sR(8, kh_nonlin, pk_nonlin)[0]
    s200_n                         = sR(200, kh_nonlin, pk_nonlin)[0]
    
    # Plot matter power spectra at z = 0
    
    # For varying baryonic matter density
    plt.loglog(kh, pk[0], color=colors[i], ls='-', label="$\Omega_m h ^2$ = {0}   \tlinear  \t   $\sigma_8$={1} \t  $\sigma200$={2}".format(par_itr[i], round(s8, 4), round(s200,4)))
    plt.loglog(kh_nonlin, pk_nonlin[0], color=colors[i], ls=':', label="$\Omega_m h ^2$ = {0}    \tnonlinear \t$\sigma_8$={1} \t$\sigma200$={2}".format(par_itr[i], round(s8_n, 4), round(s200_n,4)))

    # For varying CDM matter density with fixed mnu
    #plt.loglog(kh, pk[0], color=colors[i], ls='-', label="$\Omega_c h ^2$ = {0}   \tlinear  \t   $\sigma_8$={1} \t  $\sigma200$={2}".format(par_itr[i], round(s8, 4), round(s200,4)))
    #plt.loglog(kh_nonlin, pk_nonlin[0], color=colors[i], ls=':', label="$\Omega_c h ^2$ = {0}    \tnonlinear \t$\sigma_8$={1} \t$\sigma200$={2}".format(par_itr[i], round(s8_n, 4), round(s200_n,4)))
    
    # For varying mnu with fixed standard CDM density 
    #plt.loglog(kh, pk[0], color=colors[i], ls='-', label=r"$M_\nu$={0}   linear        $\sigma_8$={1}    $\sigma200$={2}".format(par_itr[i], round(s8, 4), round(s200,4)))
    #plt.loglog(kh_nonlin, pk_nonlin[0], color=colors[i], ls=':', label=r"$M_\nu$={0}   nonlinear  $\sigma_8$={1}    $\sigma200$={2}".format(par_itr[i], round(s8_n, 4), round(s200_n,4)))
    
plt.xlabel('k/h Mpc')
plt.ylabel('$P_m(k)$')
plt.ylim(bottom=1)
plt.legend(loc='lower left')
plt.grid(True)

# For varying baryonic matter density
plt.title('Matter Power $P_m (k)$ at z={0}\n$H_0$={1}, $\Omega_c h ^2$={2}'.format(z[0], H0, omch2))    
plt.savefig("powerspectrum_var_baryonic_density_with_sigma.pdf")                                        

# For varying CDM matter density
#plt.title('Matter power at z={0}\n$H_0$={1}, $\Omega_b h ^2$={2}'.format(z[0], H0, ombh2)) 
#plt.savefig("powerspectrum_var_cdm_with_sigma.pdf")

# For varying neutrino mass
#plt.title('Matter power at z={0}\n$H_0$={1}, $\Omega_b h ^2$={2}, $\Omega_c h ^2$={3}'.format(z[0], H0, ombh2, omch2)) 
#plt.savefig("powerspectrum_var_mnu_with_sigma.pdf")
        
    

