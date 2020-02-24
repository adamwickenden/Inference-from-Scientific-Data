# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 12:14:38 2017

@author: axw404
"""
import numpy as np
import math as m
import matplotlib.pyplot as plt

#Set lambda value
lmb = 6
#set number of counts to run for
N = 25
#define arrays to hold data
probs = np.zeros(N)
count = np.zeros(N)
#iterate through equation and fill arrays
for k in range(0,N):
    P = m.exp(-lmb)*((lmb**k)/(m.factorial(k)))
    probs[k] = P
    count[k] = k
    print(k, probs[k])

#calculate question 3c
P10 = 0
for k in range(10, probs.size):
    P10 += probs[k]
   
print(P10)

#plot and save figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.bar(count, probs)
plt.title('Probability Distribution of Number of Counts After 6 Hours')
plt.ylabel('Probability of Count')
plt.xlabel('Number of Counts')
fig.savefig("Q3.pdf", format = 'pdf')
plt.show