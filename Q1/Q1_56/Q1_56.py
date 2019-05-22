# CTA200H Final Project
# Last edited May 17, 2019 10:48pm
# Q1 part 5 and 6

import sys, platform, os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import camb
from camb import model, initialpower
from scipy.integrate import simps, trapz

#===========================================================================================
#                                       Q1 part 5                                          #
#===========================================================================================


# Defining the W factor in the integrad of the variance of the smoothed cosmic density field
def W(k, R):
    # k = wavevector
    # R = radius of the sphere (in units Mpc/h)
    return 3 * (np.sin(k*R) / (k*R) - np.cos(k*R)) / (k*R)**3
 
# Define parameter values and array of variable values
R = np.array([8., 200.]) 
k = np.linspace(1, 5 , 300)   
   
# Plot this W(k,R) for R = 8, 20
plt.clf()
plt.figure(figsize=(10,6))
plt.plot(k, W(k, R[0]), label=r"R = {0} $h^-1$ Mpc".format(R[0]))
plt.plot(k, W(k, R[1]), label=r"R = {0} $h^-1$ Mpc".format(R[1]))
plt.xlabel('k')
plt.ylabel('$W (k,R) = 3 j_1  (kR) / (kR)$')
plt.grid(True)
plt.legend(loc='lower right')
plt.title(r'$W (k,R)$')
plt.savefig('W_kR.pdf')

#===========================================================================================
#                                       Q1 part 6                                          #
#===========================================================================================

# Define the values of R
R = np.arange(1.,201.,1.)

# Get matter power spectra at z = 0
pars  = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
pars.InitPower.set_params(ns=0.965)
 
# Define Delta^2 function
def D2(k, pmk):
    # Dimensionless power spectra given Pm(k)
    return k**3 * pk / (2 * np.pi**2)
    
# Define the integrand of the variance sigma_R^2
def integrand(R, k, pmk):
    return W(k,R) * D2(k,pmk) / k
    
# Computing sigma_R^2
def sR2(R, k, pmk):
    # Return sigma_R^2 at a given R - integral from 0 to infinity of the integrand
    # R      = float, radius of sphere
    # k      = array, given by kh
    # pmk    = array corresponding to kh, matter power spectra
    return simps(integrand(R, k, pmk), k)

# Computing sigma_R
def sR(R, k, pmk):
    # square root of sR2
    return np.sqrt(sR2(R, k, pmk))
 
# Note non-linear corrections couples to smaller scales than you want
pars.set_matter_power(kmax=2.0) # default redshift = [0.0]

# Linear spectra
pars.NonLinear = model.NonLinear_none
results        = camb.get_results(pars)
kh, z, pk      = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
pk             = pk[0]
s8             = np.array(results.get_sigma8())

# Non-Linear spectra (Halofit)
pars.NonLinear                 = model.NonLinear_both
results.calc_power_spectra(pars)
kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=1, npoints = 200)
pk_nonlin                      = pk_nonlin[0]

# Plot matter power spectra at z = 0
plt.clf()
plt.figure(0)
for r in R:
    plt.plot(r, sR(r, kh, pk), color='k', marker='.')
    plt.plot(r, sR(r, kh_nonlin, pk_nonlin), color='r', marker='.')
plt.xlabel('Radius of Sphere, R (Mpc/h)')
plt.ylabel('Variance, $\sigma_R$')
plt.grid(True)
plt.legend(['linear','nonlinear'],loc='upper right')
plt.title('$\sigma_R$ at z={0}'.format(z[0]))
plt.savefig('sigma_R.pdf')

print('s8 = {0}, s200 = {1}'.format(sR(8, kh, pk), sR(200, kh, pk)))
print('s8 calculated from camb = {0}'.format(s8))
