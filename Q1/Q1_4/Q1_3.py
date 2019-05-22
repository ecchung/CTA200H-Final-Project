# CTA200H Final Project 
# last edited May 16, 2019

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
pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2)
pars.InitPower.set_params(ns=0.965)
 
# Define Delta^2 function
def D2(k, pk):
    return k**3 * pk / (2 * np.pi**2)
 
# Note non-linear corrections couples to smaller scales than you want
pars.set_matter_power(kmax=2.0) # default redshift = [0.0]

# Linear spectra
pars.NonLinear = model.NonLinear_none
results        = camb.get_results(pars)
kh, z, pk      = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
s8             = np.array(results.get_sigma8())
Delta2         = D2(kh, pk[0])

# Non-Linear spectra (Halofit)
pars.NonLinear                 = model.NonLinear_both
results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
Delta2_nonlin                  = D2(kh_nonlin, pk_nonlin[0])

# Find k where Delta2 is 1
k1 = (np.abs(Delta2 - 1.)).argmin()
k1_nonlin = (np.abs(Delta2_nonlin - 1.)).argmin()
print("Delta^2(k) = 1 at k_linear = {0} and k_nonlinear = {1}".format(kh[k1], kh_nonlin[k1_nonlin]))

# Plot matter power spectra at z = 0
plt.clf()
plt.figure(0)
plt.loglog(kh, Delta2, color='k', label="linear")
plt.loglog(kh_nonlin, Delta2_nonlin, color='r', label='non-linear')
plt.loglog([kh[k1], kh_nonlin[k1_nonlin]], [Delta2[k1], Delta2_nonlin[k1_nonlin]], color='yellow', ls='none', marker='.')
plt.xlabel('k/h Mpc')
plt.grid(True)
plt.legend(['linear','non-linear', '$\Delta ^2(k)$=1'], loc='lower right')
plt.title('Dimensionless matter power $\Delta ^2(k)$ at z={0}'.format(z[0]))
plt.savefig('Dimensionless Matter Power Spectra.pdf')



