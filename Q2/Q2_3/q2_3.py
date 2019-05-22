# CTA200H Final Project
# Question 2 part 3
# last edited May 19, 2019 4:41pm

import numpy as np
import matplotlib.pyplot as plt
from nbodykit.lab import *

# =========================================================================
#                         3 a)
# =========================================================================
# Generate 10 realizations of a Gaussian random field in a box of length L = 1000 Mpc/h on a (256)^3 mesh

cosmo = cosmology.Planck15
Plin  = cosmology.LinearPower(cosmo, redshift=0., transfer='CLASS')
Pnl   = cosmology.HalofitPower(cosmo, redshift=0.)

# initialize the 10 meshes for the gaussian realizations
seed = np.array([1,2,3,4,5,6,7,8,9,10])  # array of seed values for random generation
#seed  = np.array([1,2])

power_spectra    = []   # define list of linear power_spectra values
power_spectra_nl = []   # define list of nonlinear power_spectra values

for i in range(len(seed)):
    print(i)
    # get the linear power spectrum
    mesh = LinearMesh(Plin, Nmesh=256, BoxSize=1000, seed=seed[i])
    # preview the density field
    #plt.clf()
    #plt.figure(0)
    #plt.imshow(mesh.preview(axes=[0,1]))
    #plt.savefig("grf_seed{0}.pdf".format(seed[i]))
    r     = FFTPower(mesh, mode='1d', dk=0.005, kmin=1e-4)
    Pk    = r.power
    k     = Pk['k']
    print(np.shape(k))
    power = Pk['power'].real
    power_spectra.append(power)
    
    # get nonlinear power spectrum
    mesh = LinearMesh(Pnl, Nmesh=256, BoxSize=1000, seed=seed[i])
    rnl     = FFTPower(mesh, mode='1d', dk=0.005, kmin=1e-4)
    Pknl    = rnl.power
    knl     = Pknl['k']
    print(np.shape(knl))
    power_nl = Pknl['power'].real
    power_spectra_nl.append(power_nl)
    
color = ['k', 'r', 'b', 'orange', 'g', 'c', 'm', 'y', 'grey', 'limegreen']    
plt.clf()
#plt.figure(1)
for i in range(len(seed)):
    plt.loglog(k, power_spectra[i], color=color[i], label="Linear {0}".format(seed[i]), ls='-')
    plt.loglog(knl, power_spectra_nl[i], color=color[i], label="Nonlinear {0}".format(seed[i]), ls='--')
    
plt.title("Power Spectra")
#plt.title("Power Spectra {0}".format(seed[i]))
plt.xlabel(r"$k$ $[h Mpc^{-1}]$")
plt.ylabel(r"$P_m(k)$ $[h^{-3} Mpc{^3}]$")
plt.legend()
#plt.savefig("powerspectra_grf_seed{0}.pdf".format(seed[i]))
plt.savefig("powerspectra_grf.pdf")

power_spectraT    = np.array(power_spectra).T
mean_power        = np.mean(power_spectraT, axis=1)
power_std         = np.std(power_spectraT, axis=1)

power_spectraT_nl = np.array(power_spectra_nl).T
mean_power_nl     = np.mean(power_spectraT_nl, axis=1)
power_std_nl      = np.std(power_spectraT_nl, axis=1)

plt.clf()
plt.loglog(k, mean_power, ls='none', marker='.', color='k', label='Linear Mean')
plt.errorbar(k, mean_power, yerr=power_std, ls='none', color='k')
plt.loglog(knl, mean_power_nl, ls='none', marker='.', color='r', label='Nonlinear Mean')
plt.errorbar(knl, mean_power_nl, yerr=power_std_nl, ls='none', color='r')
plt.title("Power Spectra of 10 Gaussian Random Fields")
plt.xlabel(r"$k$ $[h Mpc^{-1}]$")
plt.ylabel(r"$P_m(k)$ $[h^{-3} Mpc{^3}]$")
plt.legend()
plt.savefig("mean_power.pdf")


