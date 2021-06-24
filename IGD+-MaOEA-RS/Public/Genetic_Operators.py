"""
Generate offspring with simulated binary crossover and polynomial mutation.

K. Deb and R. Agrawal, "Simulated binary crossover for continuous search 
space," in Complex Systems, vol. 9, no. 2, pp. 115–148, 1995.

Y. Tian, R. Cheng, X. Zhang, and Y. Jin, "PlatEMO: A MATLAB platform for 
evolutionary multi-objective optimization," in IEEE Computational Intelligence 
Magazine, vol. 12, no. 4, pp. 73-87, 2017.
"""

import numpy as np

#Performs simulated binary crossover
def simulatedBinaryCrossover(parentA, parentB, lb, ub, pc, nc):
    n = len(parentA)
    beta = np.zeros(n)
    mu = np.random.rand(n)
    beta[mu<=0.5] = (2*mu[mu<=0.5])**(1/(nc+1))
    beta[mu>0.5] = (2-2*mu[mu>0.5])**(-1/(nc+1))
    beta *= ((-1)**np.random.randint(0,2,n))
    beta[np.random.rand(n)<=0.5] = 1
    beta[np.matlib.repmat(np.random.rand()>pc,1,n)[0]] = 1
    offspringA = (parentA+parentB)/2+beta*(parentA-parentB)/2
    offspringB = (parentA+parentB)/2-beta*(parentA-parentB)/2
    offspringA = np.minimum(np.maximum(offspringA,lb),ub)
    offspringB = np.minimum(np.maximum(offspringB,lb),ub)
    return offspringA,offspringB

#Performs polynomial mutation
def polynomialMutation(individual, lb, ub, pm, nm):
    n = len(individual)
    Site = np.random.rand(n) <= pm
    mu = np.random.rand(n)
    temp = Site & (mu <= 0.5)
    individual[temp] += (ub[temp]-lb[temp])*((2*mu[temp]+(1-2*mu[temp])*(1-(individual[temp]-lb[temp])/(ub[temp]-lb[temp]))**(nm+1))**(1/(nm+1))-1)
    temp = Site & (mu > 0.5)
    individual[temp] += (ub[temp]-lb[temp])*(1-(2*(1-mu[temp])+2*(mu[temp]-0.5)*(1-(ub[temp]-individual[temp])/(ub[temp]-lb[temp]))**(nm+1))**(1/(nm+1)))
    individual = np.minimum(np.maximum(individual,lb),ub)
    return individual
