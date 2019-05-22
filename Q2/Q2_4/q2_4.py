# CTA200H Final Project
# Question 2 part 4
# last edited May 21, 2019 5:10pm

import numpy as np
import matplotlib.pyplot as plt
from nbodykit.lab import *

# =========================================================================
#                         4 a)
# =========================================================================

# Define constants and stuff
redshift = 0
cosmo    = cosmology.Planck15
Plin     = cosmology.LinearPower(cosmo, redshift, transfer='EisensteinHu')
nbar     = 1e-4 # units of h^3 Mpc^-3
nbar_str = '1e-4'

seed    = np.array([1,2,3,4,5,6,7,8,9,10])  # array of seed values for random generation
#seed     = np.array([1])

power_spectra        = []   # define list of linear power_spectra values
power_spectra_nowind = []   # define list of linear power_spectra values

for i in range(len(seed)):
    print(i)
    catalog      = LogNormalCatalog(Plin=Plin, nbar=nbar, BoxSize=1000., Nmesh=256, seed=seed[i])
    r            = FFTPower(catalog.to_mesh(resampler='cic', Nmesh=512, compensated=True), mode='1d')
    Pk           = r.power
    k            = Pk['k']
    power        = Pk['power'].real
    power_spectra.append(power)
    
    r_nowind     = FFTPower(catalog.to_mesh(), mode='1d')
    Pk_nowind    = r_nowind.power
    k_nowind     = Pk_nowind['k']
    power_nowind = Pk_nowind['power'].real
    power_spectra_nowind.append(power_nowind)    
    

color = ['k', 'r', 'b', 'orange', 'g', 'c', 'm', 'y', 'grey', 'limegreen']    
plt.clf()
for i in range(len(seed)):
    print(i)
    #plt.loglog(k, power_spectra[i], color=color[i], label="Linear {0}".format(seed[i]), ls='-')
    #plt.loglog(k_nowind, power_spectra_nowind[i], color=color[i], label="Linear No Window Deconvolution {0}".format(seed[i]), ls='--')
    plt.loglog(k, power_spectra[i], color=color[i], ls='-')
    plt.loglog(k_nowind, power_spectra_nowind[i], color=color[i], ls='--')
    
    
plt.title("Power Spectra for 10 Lognormal Density Field \n nbar = {0}".format(nbar_str))
plt.xlabel(r"$k$ $[h Mpc^{-1}]$")
plt.ylabel(r"$P_m(k)$ $[h^{-3} Mpc{^3}]$")
plt.legend(['Linear','Linear No Window Deconvolution'])
#plt.savefig("powerspectra_grf.pdf")
plt.savefig(r"powerspectra_lognormal_with_nowind_CIC_nbar{0}.pdf".format(nbar_str))

power_spectraT    = np.array(power_spectra).T
mean_power        = np.mean(power_spectraT, axis=1)
power_std         = np.std(power_spectraT, axis=1)

power_spectra_nowindT = np.array(power_spectra_nowind).T
mean_power_nowind     = np.mean(power_spectra_nowindT, axis=1)
power_std_nowind      = np.std(power_spectra_nowindT, axis=1)

plt.clf()
plt.loglog(k, mean_power, ls='none', marker='o', color='k', label='Linear Mean')
plt.errorbar(k, mean_power, yerr=power_std, ls='none', color='k', capsize=1)
plt.loglog(k_nowind, mean_power_nowind, ls='none', marker='.', color='r', label='Linear No Window Deconvolution Mean')
plt.errorbar(k_nowind, mean_power_nowind, yerr=power_std_nowind, ls='none', color='r')
plt.title("Mean Power Spectra of 10 Lognormal Density Field \n nbar= {0}".format(nbar_str))
plt.grid(True)
plt.xlabel(r"$k$ $[h Mpc^{-1}]$")
plt.ylabel(r"$P_m(k)$ $[h^{-3} Mpc{^3}]$")
plt.legend()
plt.savefig(r"mean_power_lognormal_with_nowind_CIC_nbar{0}.pdf".format(nbar_str))
#plt.savefig("mean_power_lognormal_with_nowind_tsc.pdf")