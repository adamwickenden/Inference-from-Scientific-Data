# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

#Import the data from text files
with open('Q1T.txt', 'r') as f:
    tempBase = [float(line.strip()) for line in f if line]   
with open('Q1D.txt', 'r') as f:
    dispBase = [float(line.strip()) for line in f if line]
'*******************************************************************************************************************'
'(1a)'

#Define function to calculate average
def av(data):
    return sum(data)/float(len(data))

#Define function to calculate sigma
def sigma(data):
    N = len(data)
    val = 0
    for i in range(0, N):
        val += (data[i] - av(data))**2
    return np.sqrt((1/N)*val)

#Calculate pearson co-efficient, r
var = 0
temp = tempBase[:]
disp = dispBase[:]
for i, data in enumerate(temp):
    var += (temp[i]-av(temp))*(disp[i]-av(disp))
r = (1/float(len(temp)))*var*(1/(sigma(disp)*sigma(temp)))

#Calculate students t-distribution with the pearson co-efficient value for N = 10
tp = (r*np.sqrt(len(temp)-2))/(np.sqrt(1-r**2))

#Calculate the spearman rank co-efficient, ps
#assign rank values to the temperature data
temp2 = temp[:]
temp2.sort(reverse=True)
for i, data in enumerate(temp):
    for j, value in enumerate(temp2):
        if value == data:
            temp[i] = [temp[i], j+1]
#assign rank values to the displacement data
disp2 = disp[:]
disp2.sort(reverse=True)
for i, data in enumerate(disp):
    for j, value in enumerate(disp2):
        if value == data:
            disp[i] = [disp[i], j+1]
#calculate rank difference sum
diff = 0
for i in range(0, len(temp)):
    diff += ((temp[i][1]-disp[i][1])**2)
#Calculate ps
ps = 1 - (6*diff)/(len(temp)*(len(temp)**2 - 1))

#Calculate students t-distribution with the spearman co-efficient value for N = 10
ts = (ps*np.sqrt(len(temp)-2))/(np.sqrt(1-ps**2))
del temp2
del disp2

'*******************************************************************************************************************'
'(1b)'
print_t = "r   = {0:.4g}  &  tp  = {1:.4g} \nps = {2:.4g}  &  ts  = {3:.4g}".format(r, tp, ps, ts)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.text(0.05,0.8, print_t, transform=ax.transAxes, fontsize=12)
ax.scatter(tempBase, dispBase)
plt.title('Graph Showing the Relationship Between Temperature and Displacement', y=1.08)
plt.ylabel('Displacement (Microns)')
plt.xlabel('Temperature (Celsius)')
fig.savefig("1.b.pdf", format = 'pdf', bbox_inches='tight')
plt.show

'*******************************************************************************************************************'
'(1c)'
temp = tempBase[:]
tempsq = tempBase[:]
disp = dispBase[:]
dispsq = dispBase[:]
#Perform least squares fit with temperature as the independent variable (1), m = a1x + b1
#Calulate numerator sum in the equation for average of a
a1 = 0
for i, data in enumerate(temp):
    a1 += (disp[i]-av(disp))*(temp[i]-av(temp))
    tempsq[i] = tempsq[i]**2
#Calculate average squared of temp
avsqT = av(temp)**2
#Calculated average of squares of temp
sqavT = av(tempsq)
#Calculate final a1
a1 = a1/(10*(sqavT - avsqT))

#Calculate value for b1
b1 = av(disp) - a1*av(temp)

#Create array of fit 1 and then plot fit 1 against temperature and data
M1 = []
for T in temp:
    M1.append(a1*T + b1)
print_M1 = "a1  = {0:.4g}\nb1  = {1:.4g}".format(a1, b1)
#Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.text(0.05,0.8, print_M1, transform=ax.transAxes, fontsize=12)
ax.scatter(tempBase, dispBase, label = 'Data')
ax.plot(tempBase, M1, color = 'red', label = 'Best Fit')
plt.title('Graph Showing the Best Fit Between Temperature and Displacement', y=1.08)
plt.ylabel('Displacement (Microns)')
plt.xlabel('Temperature (Celsius)')
plt.legend(loc=4, fontsize=12)
fig.savefig("1.c1.pdf", format = 'pdf', bbox_inches='tight')
plt.show

#Perform least squares fit with displacement as the independent variable (2), m = a2x + b2
#Calulate numerator sum in the equation for average of a
a2 = 0
for i, data in enumerate(disp):
    a2 += (temp[i]-av(temp))*(disp[i]-av(disp))
    dispsq[i] = dispsq[i]**2
#Calculate average squared of temp
avsqD = av(disp)**2
#Calculated average of squares of temp
sqavD = av(dispsq)
#Calculate final a1
a2 = a2/(10*(sqavD - avsqD))

#Calculate value for b1
b2 = av(temp) - a2*av(disp)

#Create array of fit 2 and then plot fit 2 against displacement and data
M2 = []
disp.sort()
for d in disp:
    M2.append(a2*d + b2)
print_M2 = "a2  = {0:.4g}\nb2  = {1:.4g}".format(a2, b2)
#Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.text(0.05,0.8, print_M2, transform=ax.transAxes, fontsize=12)
ax.scatter(dispBase, tempBase, label = 'Data')
ax.plot(disp, M2, color = 'green', label = 'Best Fit')
plt.title('Graph Showing the Best Fit Between Displacement and Temperature', y=1.08)
plt.xlabel('Displacement (Microns)')
plt.ylabel('Temperature (Celsius)')
plt.legend(loc=4, fontsize=12)
fig.savefig("1.c2.pdf", format = 'pdf', bbox_inches='tight')
plt.show

#Create array to show fit 2 rearranged to be shown in terms of temperature as independent variable
M2T = []
for T in temp:
    M2T.append(T/a2 - b2/a2)
#Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.text(0.05,0.8, print_M1, transform=ax.transAxes, fontsize=12)
ax.text(0.05,0.65, print_M2, transform=ax.transAxes, fontsize=12)
ax.scatter(tempBase, dispBase, label = 'Data')
ax.plot(tempBase, M1, color = 'red', label = 'Fit 1')
ax.plot(tempBase, M2T, color = 'green', label = 'Fit 2')
plt.title('Graph Showing Both Best Fits Between Temperature and Displacement', y=1.08)
plt.ylabel('Displacement (Microns)')
plt.xlabel('Temperature (Celsius)')
plt.legend(loc=4, fontsize=12)
fig.savefig("1.c3.pdf", format = 'pdf', bbox_inches='tight')
plt.show