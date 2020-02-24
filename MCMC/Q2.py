# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 20:17:10 2017

@author: axw404
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ

'*******************************************************************************************************************'
'(2a)'
#Calculate expectation value
a = -0.7
Ls = (1.4)*1e10
#Normalise
Norm = integ.quad(lambda x: (1/Ls)*((x/Ls)**(a))*np.exp(-x/Ls), 1e9, 1e14)
Lav = integ.quad(lambda L: ((L/Ls)**(a+1))*np.exp(-L/Ls), 1e9, 1e14 )[0]/Norm[0]


'*******************************************************************************************************************'
'(2b)'

#Input and sort data
data = 1e9*np.array([1.39, 1.40, 1.29, 5.95, 2.97, 1.63, 1.17, 2.06, 4.69, 2.48])
data = np.sort(data)
#Define array of Luminosity values
LBase = np.linspace(1e9, 5e10, 10000, endpoint = True)

#Define empty arrays to hold values of Probability density and L values
#for integration over for cumulative distributions 
cumAr = []
cumArm = []
PL = []
PLm = []
Lm = []
#Loop over the data and calulate the cumulative distribution
for i, item in enumerate(data):
    cumAr.append((i+1)/10)

#Loop over a large number L values for the model creating a predicted cumulative distribution
for l in LBase:
    PLm.append((1/Ls)*((l/Ls)**(a))*np.exp(-l/Ls))
    Lm.append(l)
    cumArm.append(integ.simps(PLm, Lm)/Norm[0])
    
'*******************************************************************************************************************'
'(2c)'
#Using arrays generated previously in 2b, for each data point find the index of the nearest value in the model's
#luminosity array. Then using that index calculate the difference between the cumulative distribution values of
#the model and the data. Then find the max value of the differences.
DAr = []
j = 0
for d in data:
    k = min(range(len(Lm)), key=lambda i: abs(Lm[i]-d))
    q = cumAr[j]
    w = cumArm[k]
    DAr.append(abs(q-w))
    j += 1
D = max(DAr)

print_D = "D = {0:.4g}".format(D)
    
#Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(Lm, cumArm, color = 'red', label = 'Model')
ax.step(data, cumAr,  label = 'Data')
ax.plot([Lav, Lav],[0, 1], color = 'black', linestyle = '--', label = '<L>')
ax.text(0.76,0.35, print_D, transform=ax.transAxes, fontsize=12)
ax.set_ylim(bottom=0)
plt.title('Graph Showing Cumulative Distribution of Data \n Against Model Prediction')
plt.ylabel('Cumulative Distribution (Arb)')
plt.xlabel('Luminisoty (' r'$L_{\odot}$' ')')
plt.legend(loc=4, fontsize=12)
fig.savefig("2.b.pdf", format = 'pdf', bbox_inches='tight')
plt.show
