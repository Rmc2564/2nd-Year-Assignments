# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 19:15:55 2023

@author: alexa
"""

'Simulatews the random radioactive decay of r nuclei A, B, C with the following'
'half lives:    A -> B, t = 10.1 hrs. B->C, t = 15.7hrs'
'Under neutron flux, C goes back to A with an effective half life of 3.2 hours'

import numpy
from matplotlib import pyplot as plt
import random

def has_transitioned(prob):
    r = random.random()
    if r<prob:
        return True
    else:
        return False
    
def evolveOne(currentState, rules):
    for (initial,fin,prob) in rules:
        if currentState == initial:
            if has_transitioned(prob):
                return fin
            else:
                return initial
    if currentState != initial:
        return currentState
    
def evolveMany(states, rules):
    newState = []
    for i in range(len(states)):
        state = states[i]
        newstate = evolveOne(state, rules)
        newState.append(newstate)
    return newState

def evolve_system(NA, NB, NC, rules, n_step):  #Given a system of A,B,C this will simulate how it changes with time
    state = (['A'] * NA)+(['B'] * NB)+(['C'] * NC)


    A_count = numpy.empty(n_step + 1, dtype=int)
    B_count = numpy.empty(n_step + 1, dtype=int)
    C_count = numpy.empty(n_step + 1, dtype=int)
    
    A_count[0] = NA
    B_count[0] = NB
    C_count[0] = NC
    
    for i in range(1,n_step+1):
        state = evolveMany(state,rules)
        A_count[i] = state.count('A')
        B_count[i] = state.count('B')
        C_count[i] = state.count('C')
    return A_count, B_count, C_count

def simulate():
    nsteps = 200
    t_total = 100
    dt = t_total/nsteps
    t_half_A = 10.1
    t_half_B = 15.7
    t_half_C = 3.2

    countA = numpy.empty(401, dtype = 'int')
    countB = numpy.empty(401, dtype = 'int')
    countC = numpy.empty(401, dtype = 'int')

    probA = (numpy.log(2)*dt)/t_half_A
    probB = (numpy.log(2)*dt)/t_half_B
    probC = (numpy.log(2)*dt)/t_half_C

    rules = [('A','B',probA),('B','C',probB),('C','A',probC)]


    countA1, countB1, countC1 = evolve_system(0,0,250,rules,200)

    rules = [('A','B',probA),('B','C',probB)]

    countA2, countB2, countC2 = list(evolve_system(countA1[-1],countB1[-1],countC1[-1],rules,200))

    for i in range(len(countA1)):
        countA[-1*i] = countA2[-1*i]  #must be in this order so the 0th index isn't changed
        countB[-1*i] = countB2[-1*i]
        countC[-1*i] = countC2[-1*i]

        countA[i] = countA1[i]
        countB[i] = countB1[i]
        countC[i] = countC1[i]
    
    counts = (countA, countB, countC)
    return counts

time = numpy.linspace(0,200,401)
counts = simulate()
    
fig = plt.figure(figsize = (1,2))
ax = fig.add_subplot(111)
plt.plot(time,counts[0], label = 'A')
plt.plot(time,counts[1], label = 'B')
plt.plot(time,counts[2], label = 'C')
plt.legend()
plt.xlabel('Time (hours)')
plt.ylabel('Number of nuclei')
plt.title('Plot of the numbers of A, B and C nuclei initially under neutron flux (turned off at 100hrs)')
plt.show()


nsim = 20
counts = numpy.zeros((nsim,401))
print(numpy.shape(counts))
for k in range(nsim):
    count = simulate()
    count = count[0]
    counts[k] = count
    
averages = numpy.average(counts,axis = 0) #uncertainty in A at each time step
errors = numpy.std(counts, axis = 0)
fig2 = plt.figure(figsize = (1,2))
ax = fig2.add_subplot(111)
plt.errorbar(time,averages, yerr = errors)
plt.xlabel('Time (hours)')
plt.ylabel('Average number of A atoms')
plt.title('The average number of atoms of A over 20 simulations, with the error given as the standard deviation taken as the error')
plt.show()