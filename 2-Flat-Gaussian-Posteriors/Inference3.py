# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 15:25:13 2017

@author: axw404
"""

import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt

#define constants and inputs for equation
rmin = 1e24
rmax = 1e26
n = int(1e6)
E = 1e44
N = 8
k = 1
r = np.linspace(rmin, rmax, n, endpoint = True)

#define funtional and integrate in sensible region to yield normalising factor
denominator = lambda x: (np.exp(-((E/1e-9)*(1/(4*np.pi*(x**2)))))*(((E/1e-9)*
                                  (1/(4*np.pi*(x**2))))**N)*(4*k*np.pi*(x**2)))/(np.math.factorial(N))
normalise = integ.quad(denominator, 1e24, 1e26)

#Define N(average) and Posterior 
Nav = (E/1e-9)*(1/(4*np.pi*(r**2)))
P = (np.exp(-Nav)*(Nav**N)*(4*k*np.pi*(r**2)))/((np.math.factorial(N))*normalise[0])

#Calculate modal value
#np.argmax returns index value of maximum value in array
minusmode = r[(np.argmax(P))- 1]
mode = r[(np.argmax(P))]
plusmode = r[(np.argmax(P))+ 1]

#calculate mean
#sum all the probabilities multiplied by the valuees, then multiply by the density of each slice of distribution
total = 0
for i in range(0,n):
    total += P[i] * r[i]
mean = total*((rmax-rmin)/n)

#Sigma estimate from calculation
a = ((E/1e-9)*(1/(4*np.pi)))
sigma = ((6*a/(mode**4)) - (2/(mode**2)*(N-1)))**(-0.5)

print_string = "Mode  = {0:.4g} \nMean  = {1:.4g} \nSigma = {2:.4g}".format(mode, mean, sigma)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(r, P, linewidth=2)
ax.text(0.6,0.75, print_string, transform=ax.transAxes, fontsize=12)
plt.title('Posterior Probability Density Function for N=8', y=1.08)
plt.ylabel('Posterior Density P(N|r,E,GRB)')
plt.xlabel('Distance to GRB (m)')
ax.set_ylim(bottom=0)
L1 = plt.axvline(x=mode, linestyle='--', linewidth=1, label="Mode")
L2 = plt.axvline(x=mode-sigma, linestyle=':', linewidth=1.5, label="Sigma")
L3 = plt.axvline(x=mode+sigma, linestyle=':', linewidth=1.5, label="L3")
plt.legend(handles=[L1,L2], loc=4)
fig.savefig("Q3.pdf", format = 'pdf')
plt.show