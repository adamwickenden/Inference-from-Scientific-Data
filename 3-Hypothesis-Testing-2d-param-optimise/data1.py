# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import scipy.integrate as integ
import matplotlib.pyplot as plt
import urllib.request as ur

#Define dataset URL's
N1 = 10
url1 = "http://www.sr.bham.ac.uk/~imandel/Inference/data_10hr.txt"
N2 = 24
url2 = "http://www.sr.bham.ac.uk/~imandel/Inference/data_24hr.txt"
N3 = 100
url3 = "http://www.sr.bham.ac.uk/~imandel/Inference/data_100hr.txt"

#Choose Dataset
key = input("Choose dataset 1,2 or 3: " )
if(key=='1'):
    N=N1
    url=url1
if(key=='2'):
    N=N2
    url=url2
if(key=='3'):
    N=N3
    url=url3

#Import data, strip whitespace characters and convert to float numbers instead of strings
im = ur.urlopen(url)
data = np.zeros(N)
for i, line in enumerate(im):
    line = float(line.decode('utf8').strip('b\n'))
    data[i] = line  

#Initialise given variables
B = 5
t = np.linspace(0,N, N, endpoint = True)

#Calculate the constant in all 3 probabilities
C = (1/(np.sqrt(2*np.pi)))**N

#Calculate PM0
PM0 = (C)*(np.exp((-0.5)*np.sum((data-B)**2)))

#Calculate PM1, create array of values for different L0's, 
#then use Simpsons integration on the data array
PMa1 = []
L0 = np.linspace(0, 10, 200, endpoint = True)
for l in L0:
    expon = 1
    for d in data:
        expon *= np.exp((-0.5)*(d - B - l)**2)
    PMa1.append(expon*C)
PM1 = integ.simps(PMa1, L0)/10

#Calculate PM2, cover all possible values, for each A0, do simpsons 
#integration across Tau, and then integrate collapsed data array 
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
#Plot array of probailities for each L0, as calculated above for PM1, against the values of L0,
#normalised using the above calculated value, PM1, multiplied by 10 to cancel out tophat normalisation
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(L0, (PMa1/(PM1*10)), linewidth=2)
plt.title('Marginalised Pobability Density Function for M1', y=1.08)
plt.ylabel('Posterior Density P(L0|data,M1)')
plt.xlabel('Excess Luminosity, L0 (arb)')
ax.set_ylim(bottom=0)
fig.savefig("Q3.1-M1.pdf", format = 'pdf')
plt.show

#Find most likely L0
modeL0 = np.argmax(PMa1)/20

#Plot Marginalised Posterior Distributions for M2 for A0 and Tau and fine most probable A0 and Tau

#Marginalise out Tau and plot against A0, just use the calculation from above
#Using the array of collapsed data for A0 from the calculation of PM2, plot it agianst all A0,
#Normalise using PM2 multiplied by the relevant tophat normalisation constant
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(A0, (PMa2a/(PM2*10)), linewidth=2)
plt.title('Marginalised Pobability Density Function for A0 in M2', y=1.08)
plt.ylabel('Posterior Density P(A0|data,M2)')
plt.xlabel('Initial Amplitude, A0 (arb)')
ax.set_ylim(bottom=0)
fig.savefig("Q3.1-M2-A0.pdf", format = 'pdf')
plt.show

#Find most likely A0
modeA0 = np.argmax(PMa2a)/20

#Marginalise out A0 and plot against Tau, also create 2d array for contour plot
#Do simpsons integration across A0 for each possible value of Tau, then plot this against Tau
#Normalise using PM2 multiplied by the relevant tophat normalisation constant
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

#Plot Tau posterior probability density
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(TAU, (PMa2t/(PM2*10)), linewidth=2)
plt.title('Marginalised Pobability Density Function for Tau in M2', y=1.08)
plt.ylabel('Posterior Density P(Tau|data, M2)')
plt.xlabel('Decay Constant, Tau (Hours)')
ax.set_ylim(bottom=0)
fig.savefig("Q3.1-M2-Tau.pdf", format = 'pdf')
plt.show

#plot joint probability density
fig = plt.figure()
ax = fig.add_subplot(111)
plot = ax.pcolor(X, Y, PM2ta/(PM2*990))
fig.colorbar(plot, ax=ax)
plt.title('Joint Pobability Density Function for M2', y=1.08)
plt.ylabel('Decay Constant, Tau (Hours)')
plt.xlabel('Initial Amplitude, A0 (arb)')
fig.savefig("Q3.1-M2-Tau&A0.pdf", format = 'pdf')
plt.show

#Find most likely Tau
modeTau = np.argmax(PMa2t)/2

print_modes = "L0  = {0:.3g} \nA0  = {1:.3g} \nTau = {2:.3g}".format(modeL0, modeA0, modeTau)
print(print_modes)



