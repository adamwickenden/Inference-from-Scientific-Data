# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 17:11:24 2017

@author: axw404
"""
#
#
import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt

#Set upper limit to flat posterior, define number of intervals over which to plot function
# and define variance of gaussian prior
Smax = 30
n = int(1e6)
sigma = 0.3

#Define the function in the normalising numerator, integrate using Scipy.Integrate.Quad
denominator = lambda x: np.exp(-10*x)*(10*x**50)*np.exp((-(x-7.7)**2)/(2*(sigma**2)))
normalise = integ.quad(denominator,0,Smax)

#define array to hold range
S = np.linspace(0, Smax, n, endpoint = True)
#define posterior probability density function
P = ((10*S**50)*np.exp(-10*S)*(np.exp((-(S-7.7)**2)/(2*(sigma**2)))))/normalise[0]

#Calculate modal value
#np.argmax returns index value of maximum value in array
minusmode = S[(np.argmax(P))- 1]
mode = S[(np.argmax(P))]
plusmode = S[(np.argmax(P))+ 1]
  
#calculate mean, sum all the probabilities multiplied by the valuees, then multiple by the density of each point
total = 0
tot = 0
for i in range(0,n):
    total += S[i] * P[i]
    tot += P[i]*(Smax/n)
mean = total*(Smax/n)

print_string = "Mode = {0:.4g} \nMean = {1:.4g}".format(mode, mean)
#plot and save figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(S, P, linewidth=2)
ax.text(0.7,0.8, print_string, transform=ax.transAxes, fontsize=12)
plt.title('Posterior Probability Density Function with Updated Data for n=50')
plt.ylabel('Posterior Density P(S|n)')
plt.xlabel('Density of Stars per Square Degree (S)')
ax.set_ylim(bottom=0)
fig.savefig("Q2.3.pdf", format = 'pdf')
plt.show