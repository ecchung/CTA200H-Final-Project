# CTA200H Final Project 
# last edited May 16, 2019
# Question 1.2 final version

import sys, platform, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower

# Get matter power spectra at z = 0
pars  = camb.CAMBparams()
H0    = 67.5
ombh2 = 0.022
omch2 = 0.122
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2) # default : H0=67.5, ombh2=0.022, omch2=0.122
pars.InitPower.set_params(ns=0.965)
 
 # ombh2  = (float64) Omega_baryon h^2, (phyiscal matter density) ^2
 # omch2  = (float64) Omega_cdm h^2,    (cold dark matter density) ^2
 # omk    = (float64) Omega_K,          just set as 0
 # omnuh2 = (float64) Omega_massive_neutrino h^2
 # H0     = (float64) Hubble parameter is km/s/Mpc units

# Note non-linear corrections couples to smaller scales than you want
pars.set_matter_power(kmax=2.0) # default redshits = [0.0]

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
# Linear scale plot
plt.clf()
plt.figure(0)
plt.loglog(kh, pk[0], color='k', label="linear")
plt.loglog(kh_nonlin, pk_nonlin[0], color='r', label='non-linear')
plt.xlabel('k/h Mpc');
plt.legend(loc='lower left')
plt.grid(True)
plt.title('Matter power at z={0}\nH0 = {1}, ombh2 = {2}, omch2 = {3}\nlogscale'.format(z[0], H0, ombh2, omch2))
plt.savefig("Q1_plots/powerspectrum_H0={0}_logscale.pdf".format(H0))

# Log scale plot
plt.clf()
plt.figure(1)
plt.plot(kh, pk[0], color='k', label="linear")
plt.plot(kh_nonlin, pk_nonlin[0], color='r', label='non-linear')
plt.xlabel('k/h Mpc');
plt.legend(loc='upper right')
plt.grid(True)
plt.title('Matter power at z={0}\nH0 = {1}, ombh2 = {2}, omch2 = {3}\nlinearscale'.format(z[0], H0, ombh2, omch2))
plt.savefig("Q1_plots/powerspectrum_H0={0}_linearscale.pdf".format(H0))
