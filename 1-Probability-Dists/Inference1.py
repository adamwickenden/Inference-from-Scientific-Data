# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

#No of days to run for
N = 100

#Initial Probabililties
Ps = 1
Pw = 0
#array of sunny probs
sunny = np.zeros(N+1)
sunny[0] = Ps

#loop across array
for n in range(1, sunny.size):
    
    #calculate probabilities given different initial conditions
    Ps_N = Ps*0.9 + Pw*0.4
    Pw_N = Ps*0.1 + Pw*0.6
    #Update probabilities
    Ps = Ps_N
    Pw = Pw_N
    #fill array
    sunny[n] = Ps

#plot array
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(sunny)
plt.title('Probability of Day Being Sunny After N Days')
plt.ylabel('Probability of Day Being Sunny')
plt.xlabel('Number of Days')
fig.savefig("Q1.pdf", format = 'pdf')
plt.show