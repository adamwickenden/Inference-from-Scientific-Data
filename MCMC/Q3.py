# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 14:38:06 2017

@author: axw404
"""

import numpy as np
import matplotlib.pyplot as plt
import random

#Import the data from text file
with open('Q3Data.txt', 'r') as f:
    vBase = [float(line.strip()) for line in f if line]  
  
'*******************************************************************************************************************'
'(3a)'
#Data is gaussian distributed so minimise S over a range of means for the gaussian
means = np.linspace(9695, 9705, 1000)
S = []
for m in means:
    s = 0
    for v in vBase:
        s += ((v-m)**2)
    S.append(s)
#find minimum value of chi and the value of v that it corresponds to
chi = np.min(S)
minchiv = means[np.argmin(S)]

#From Avni 1976 the 90% confidence interval of a single model parameter is given by S(a) = S(min) + 2.71
chi90 = chi + 2.71
#array of chi's from index 0-476
chidown = S[:np.argmin(S)-1]
#array of chi's from index 478-1000
chiup = S[np.argmin(S)+1: ]

#search arrays for nearest values
i = min(range(len(chidown)), key=lambda i: abs(chidown[i]-chi90))
j = min(range(len(chiup)), key=lambda i: abs(chiup[i]-chi90))

#find corresponding velocity value
vdown = means[i]
vup = means[477+j]

print_v = "\u03C7\u00b2 = {0:.5g}\nv  = {1:.5g} ".format(chi, minchiv)
print_bound = "90% Interval: {0:.5g} < v < {1:.5g}".format(vdown, vup)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(means[400:554], S[400:554], label = 'Fitting Data')
ax.text(0.39,-0.32, print_v, transform=ax.transAxes, fontsize=12)
ax.text(0.2,-0.39, print_bound, transform=ax.transAxes, fontsize=12)
ax.plot([means[400],means[554] ],[chi90, chi90], color = 'yellow', linestyle = '-', label = '\u03C7\u00b2 Intercept')
ax.plot([vup, vup],[60, 100], color = 'red', linestyle = '--', label = '90% Interval')
ax.plot([vdown, vdown],[60, 100], color = 'red', linestyle = '--')
ax.set_ylim(bottom=60)
plt.title('Graph Showing ' r'${\chi}^{2}$' ' Value for Different Star Velocities ')
plt.ylabel(''r'${\chi}^{2}$''(Arb)')
plt.xlabel('Star Velocity''(' r'$ms^{-1}$' ')''')
plt.legend(loc=1, fontsize=12)
fig.savefig("3.a.pdf", format = 'pdf', bbox_inches='tight')
plt.show()        

'*******************************************************************************************************************'
'(3b)'
#define probability function
def P(value):
    N = len(vBase)
    if (value > 9705):
        return 0
    elif (value < 9675):
        return 0
    else:
        sums = 1
        for item in vBase:
            sums *= np.exp((-0.5)*(item-value)**2)
        return (1/10)*((1/np.sqrt(2*np.pi))**N)*sums

Probs = []
Vels = []
#Random starting value
vi = ((random.random()*10)+9695)

#burn in
sigmaj = 1
acc = 0
j = 0
for i in range(0,10000):
    Pi = P(vi)    
    vj = np.random.normal(vi, sigmaj)
    Pj = P(vj)
    a = (Pj/Pi)
    u = random.random()
    if u < a:
        vi = vj
        acc += 1
    else:
        vi = vi
    if(j==100):
        j = 0
        ar = acc/i
        sigmaj = sigmaj*(ar/0.25)
    j += 1
    

#Chain
for j in range(0,1000000):
        Pi = P(vi)    
        vj = np.random.normal(vi, sigmaj)
        Pj = P(vj)
        a = (Pj/Pi)
        u = random.random()
        if u < a:
            Probs.append(Pj)
            Vels.append(vj)
            vi = vj
        else:
            vi = vi
#Mean velocity
vmean = sum(Vels)/len(Vels)

#Calculate fraction of PDF contained within 90% interval
i = 0
for v in Vels:
    if(v > vdown and v < vup):
        i += 1
frac = i/len(Vels)


fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(Vels, Probs, label = 'Probabilities', s = 1)
ax.set_ylim([0.8*np.min(Probs),1.2*np.max(Probs)])
plt.title('Scatter Showing the Non-Normalised PDF of the Stars Velocity ', y=1.05)
plt.ylabel('P(v) (Arb)')
plt.xlabel('Velocity''(' r'$ms^{-1}$' ')''')
plt.legend(loc=2, fontsize=12)
fig.savefig("3.b.1.pdf", format = 'pdf', bbox_inches='tight')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
plot = ax.hist(Vels, bins = 100, normed = True, label = 'Frequency')
plt.title('Histogram Showing the Normalised PDF of the Stars Velocity  ')
ax.plot([vmean, vmean],[0, max(plot[0])], color = 'red', linestyle = '-', label = 'Mean')
ax.plot([vup, vup],[0, max(plot[0])], color = 'red', linestyle = '--', label = '90% Confidence')
ax.plot([vdown, vdown],[0, max(plot[0])], color = 'red', linestyle = '--')
plt.ylabel('P(v) (Arb)')
plt.xlabel('Velocity''(' r'$ms^{-1}$' ')''')
plt.legend(loc=2, fontsize=12)
fig.savefig("3.b.2.pdf", format = 'pdf', bbox_inches='tight')
plt.show()

'*******************************************************************************************************************'
'(3d)'
def P2d(value, sigma):
    N = len(vBase)
    if (value >= 9705):
        return 0
    elif (value <= 9675):
        return 0
    elif (sigma >= 3):
        return 0
    elif (sigma <= 0):
        return 0
    else:
        sums = 1
        for item in vBase:
            sums *= np.exp(((-0.5)*(item-value)**2)/(sigma**2))
        return (1/10)*((1/(np.sqrt(2*np.pi)*sigma)**N)*sums)
    
Probs2D = []
Points2D = []
v2di = ((random.random()*10)+9695)
sigmai = (random.random()*3)

sigmajv = 0.2
sigmajs = 0.1

vacc = 0
l = 0
#Burn-in for v
for j in range(0,10000):
    #Step in v
    P2di = P2d(v2di,sigmai)   
    v2dj = np.random.normal(v2di, sigmajv)
    P2dj = P2d(v2dj, sigmai)
    a = (P2dj/P2di)
    u = random.random()
    if u < a:
        v2di = v2dj
        vacc += 1
    else:
        v2di = v2di
    
    if (l==100):
        var = vacc/j
        sigmajv = sigmajv*(var/0.25)
        l=0
    l += 1

sacc = 0
l = 0
#Burn in for sigma
for j in range(0,10000):
    #Step in sigma
    P2di = P2d(v2di,sigmai)   
    sigmaj = np.random.normal(sigmai, sigmajs)
    P2dj = P2d(v2di, sigmaj)
    a = (P2dj/P2di)
    u = random.random()
    if u < a:
        sigmai = sigmaj
        sacc += 1
    else:
        sigmai = sigmai
        
    if (l==100):
        sar = sacc/j
        sigmajs = sigmajs*(sar/0.25)
        l=0
    l += 1

#Chain
for j in range(0,1000000):
    #Step in v
    P2di = P2d(v2di,sigmai)   
    v2dj = np.random.normal(v2di, sigmajv)
    P2dj = P2d(v2dj, sigmai)
    a = (P2dj/P2di)
    u = random.random()
    if u < a:
        Probs2D.append(P2dj)
        Points2D.append([v2dj, sigmai])
        v2di = v2dj
    else:
        v2di = v2di
        
    #Step in sigma
    P2di = P2d(v2di,sigmai)   
    sigmaj = np.random.normal(sigmai, sigmajs)
    P2dj = P2d(v2di, sigmaj)
    a = (P2dj/P2di)
    u = random.random()
    if u < a:
        Probs2D.append(P2dj)
        Points2D.append([v2di, sigmaj])
        sigmai = sigmaj
    else:
        sigmai = sigmai

Vels2D = [x[0] for x in Points2D]
Sigmas2D = [x[1] for x in Points2D]

hist2d = np.histogram2d(Vels2D, Sigmas2D, bins = 100, normed = True)
velocity = np.histogram(Vels2D, bins = 100, normed = True)
sigma = np.histogram(Sigmas2D, bins = 100, normed = True)

#Marginalise velocity from 2d histogram
velocitymarg = []
velstep = hist2d[1][2] - hist2d[1][1]
for i in range(0,100):
    value = 0
    for j in range(0,100):
        value += hist2d[0][i][j]
    velocitymarg.append(value)
velocitymarg = np.array(velocitymarg)/(velstep * sum(velocitymarg))

#Velocity mean:
vmean2 = sum(Vels2D)/len(Vels2D)
#Velocity standard deviation
value = 0
for v in Vels2D:
    value += (v - vmean2)**2
vdev2 = np.sqrt(value/(len(Vels2D)))

#Marginalise sigma from 2d histogram
sigmamarg = []
sigstep = hist2d[2][2] - hist2d[2][1]
for i in range(0,100):
    value = 0
    for j in range(0,100):
        value += hist2d[0][j][i]
    sigmamarg.append(value/100)
sigmamarg = np.array(sigmamarg)/(sigstep * sum(sigmamarg))

#Sigma mean:
smean = sum(Sigmas2D)/len(Sigmas2D)
#Velocity standard deviation
value = 0
for s in Sigmas2D:
    value += (s - smean)**2
sdev = np.sqrt(value/(len(Sigmas2D)))

#Chi squared using new sigma value
Chinew = 0
for v in vBase:
    Chinew += ((v-vmean2)**2)/(smean**2)


fig = plt.figure()
ax = fig.add_subplot(111)
plot2 = ax.hist2d(Vels2D, Sigmas2D, bins = 100, normed = True)
fig.colorbar(plot2[3], ax=ax)
plt.title('Histogram Showing the Normalised 2D PDF of Velocity vs Sigma')
plt.ylabel('Sigma''(' r'$ms^{-2}$' ')''')
plt.xlabel('Velocity''(' r'$ms^{-1}$' ')''')
fig.savefig("3.d.1.pdf", format = 'pdf', bbox_inches='tight')
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(Vels2D, bins = 100, normed = True, label = 'Raw Histogram')
ax.step(hist2d[1][1:101], velocitymarg, color = 'green', label = 'Marginalised')
ax.plot([vmean2, vmean2],[0, max(velocitymarg)], color = 'red', linestyle = '-', label = 'Mean')
ax.plot([vmean2 + vdev2, vmean2 + vdev2],[0, max(velocitymarg)], color = 'red', linestyle = '--', label = "1$\sigma$")
ax.plot([vmean2 - vdev2, vmean2 - vdev2],[0, max(velocitymarg)], color = 'red', linestyle = '--')
ax.plot([vmean2 + 2*vdev2, vmean2 + 2*vdev2],[0, max(velocitymarg)], color = 'red', linestyle = ':', label = "2$\sigma$")
ax.plot([vmean2 - 2*vdev2, vmean2 - 2*vdev2],[0, max(velocitymarg)], color = 'red', linestyle = ':')
plt.title('Histogram Showing the Normalised Marginalised PDF of Velocity')
plt.ylabel('P(v) (arb)')
plt.xlabel('Velocity''(' r'$ms^{-1}$' ')''')
ax.set_xlim([ (9.7*1e3)-1, max(Vels2D)])
plt.legend(loc=2, fontsize=12)
fig.savefig("3.d.2.pdf", format = 'pdf', bbox_inches='tight')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(Sigmas2D, bins = 100, normed = True, label = 'Raw Histogram')
ax.step(hist2d[2][1:101], sigmamarg, color = 'green', label = 'Marginalised')
ax.plot([smean, smean],[0, max(sigmamarg)], color = 'red', linestyle = '-', label = 'Mean')
ax.plot([smean + sdev, smean + sdev],[0, max(sigmamarg)], color = 'red', linestyle = '--', label = '1$\sigma$')
ax.plot([smean - sdev, smean - sdev],[0, max(sigmamarg)], color = 'red', linestyle = '--')
ax.plot([smean + 2*sdev, smean + 2*sdev],[0, max(sigmamarg)], color = 'red', linestyle = ':', label = '2$\sigma$')
ax.plot([smean - 2*sdev, smean - 2*sdev],[0, max(sigmamarg)], color = 'red', linestyle = ':')
ax.set_xlim([min(Sigmas2D)-0.1 ,1.8])
plt.title('Histogram Showing the Normalised Marginalised PDF of Sigma')
plt.ylabel('P(v) (arb)')
plt.xlabel('Sigma''(' r'$ms^{-1}$' ')''')
plt.legend(loc=1, fontsize=12)
fig.savefig("3.d.3.pdf", format = 'pdf', bbox_inches='tight')
plt.show()
