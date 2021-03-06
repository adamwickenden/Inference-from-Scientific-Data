# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 19:39:37 2017

@author: axw404
"""


import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt
import urllib.request as ur

#Import data, strip whitespace characters and convert to float numbers instead of strings
#Define number of data points taken
N = 24
url = "http://www.sr.bham.ac.uk/~imandel/Inference/data_24hr.txt"
im = ur.urlopen(url)
data = np.zeros(N)
for i ,line in enumerate(im):
    line = float(line.decode('utf8').strip('b\n'))
    data[i] = line    

#Initialise given variables
B = 5
t = np.linspace(0,N, N, endpoint = True)

#Calculate the constant in all 3 probabilities
C = (1/(np.sqrt(2*np.pi)))**N
#Calculate PM0
PM0 = (C)*(np.exp((-0.5)*np.sum((data-B)**2)))

#Calculate PM1
PMa1 = []
L0 = np.linspace(0, 10, 200)
for l in L0:
    expon = 1
    for d in data:
        expon *= np.exp((-0.5)*(d - B - l)**2)
    PMa1.append(expon*C)

PM1 = integ.simps(PMa1, L0)/10

#Calculate PM2
PMa2t = []
PMa2a = []
A0 = np.linspace(0, 10, 200, endpoint = True)
TAU = np.linspace(1, 100, 200, endpoint = True)
for a in A0:
    PMa2t = []
    for tau in TAU:
        expon1 = 1
        for j, d in enumerate(data):
           expon1 *= np.exp((-0.5)*(d - B - a*np.exp(-j/tau))**2)
           
        PMa2t.append(expon1*C)
        
    PMa2a.append(integ.simps(PMa2t, TAU)/99)
    
PM2 = integ.simps(PMa2a, A0)/10

#Calculate Odds Ratios
O10 = PM1/PM0
O20 = PM2/PM0
O21 = PM2/PM1
print_odds = "O10 = {0:.3g} \nO20 = {1:.3g} \nO21 = {2:.3g}\n".format(O10, O20, O21)
print(print_odds)

#Plot Marginalised Posterior Distribution for M1 and find most probable L0:
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(L0, (PMa1/(PM1*10)), linewidth=2)
plt.title('Marginalised Pobability Density Function for M1', y=1.08)
plt.ylabel('Posterior Density P(L0|data,M1)')
plt.xlabel('Excess Luminosity, L0 (arb)')
ax.set_ylim(bottom=0)
fig.savefig("Q3.2-M1.pdf", format = 'pdf')
plt.show

modeL0 = np.argmax(PMa1)/20

#Plot Marginalised Posterior Distributions for M2 for A0 and Tau and fine most probable A0 and Tau

#Marginalise out Tau and plot against A0, just use the calculation from above
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(A0, (PMa2a/(PM2*10)), linewidth=2)
plt.title('Marginalised Pobability Density Function for A0 in M2', y=1.08)
plt.ylabel('Posterior Density P(A0|data,M2)')
plt.xlabel('Initial Amplitude, A0 (arb)')
ax.set_ylim(bottom=0)
fig.savefig("Q3.2-M2-A0.pdf", format = 'pdf')
plt.show

modeA0 = np.argmax(PMa2a)/20

#Marginalise out A0 and plot against Tau, also create 2d array for contour plot
PMa2t = []
PM2ta = np.zeros((len(TAU), len(A0)))
X, Y = np.meshgrid(A0, TAU)
for i, tau in enumerate(TAU):
    PMa2a = []
    for a in A0:
        expon1 = 1
        for j, d in enumerate(data):
           expon1 *= np.exp((-0.5)*(d - B - a*np.exp(-j/tau))**2)         
        PMa2a.append(expon1*C)
    PMa2t.append(integ.simps(PMa2a, A0)/10)
    PM2ta[i] = PMa2a
    
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(TAU, (PMa2t/(PM2*10)), linewidth=2)
plt.title('Marginalised Pobability Density Function for Tau in M2', y=1.08)
plt.ylabel('Posterior Density P(Tau|data, M2)')
plt.xlabel('Decay Constant, Tau (Hours)')
ax.set_ylim(bottom=0)
fig.savefig("Q3.2-M2-Tau.pdf", format = 'pdf')
plt.show

fig = plt.figure()
ax = fig.add_subplot(111)
plot = ax.pcolor(X, Y, PM2ta/(PM2*990))
fig.colorbar(plot, ax=ax)
plt.title('Joint Pobability Density Function for M2', y=1.08)
plt.ylabel('Decay Constant, Tau (Hours)')
plt.xlabel('Initial Amplitude, A0 (arb)')
fig.savefig("Q3.2-M2-Tau&A0.pdf", format = 'PDF')
plt.show

modeTau = np.argmax(PMa2t)/2

print_modes = "L0  = {0:.3g} \nA0  = {1:.3g} \nTau = {2:.3g}".format(modeL0, modeA0, modeTau)
print(print_modes)