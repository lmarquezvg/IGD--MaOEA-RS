"""
Reference set generated with the systematic approach proposed by Das and Dennis.

I. Das and J. E. Dennis, "Normal-Boundary Intersection: A New Method for 
Generating the Pareto Surface in Nonlinear Multicriteria Optimization Problems," 
in SIAM Journal on Optimization, vol. 8, no. 3, pp. 631-657, 1998.

Y. Tian, R. Cheng, X. Zhang, and Y. Jin, "PlatEMO: A MATLAB platform for 
evolutionary multi-objective optimization," in IEEE Computational Intelligence 
Magazine, vol. 12, no. 4, pp. 73-87, 2017.
"""

import numpy as np
import numpy.matlib
import itertools
import math

#Returns a matrix with all combinations of the elements of v taken k at a time
def nchoosek(v,k):
    return np.array(list(itertools.combinations(v,k)))

#Generates a set of uniformly distributed points on the unit hyperplane
def uniformReferencePoints(N,m):
    H1 = 1
    while (math.comb(H1+m,m-1) <= N):
        H1 += 1
    W = nchoosek(list(range(1,H1+m)),m-1)-np.matlib.repmat(np.arange(0,m-1),math.comb(H1+m-1,m-1),1)-1
    W = (np.concatenate((W,np.zeros([len(W),1])+H1),axis=1)-np.concatenate((np.zeros([len(W),1]),W),axis=1))/H1
    if (H1 < m):
        H2 = 0
        while (math.comb(H1+m-1,m-1)+math.comb(H2+m,m-1) <= N):
            H2 += 1
        if (H2 > 0):
            W2 = nchoosek(list(range(1,H2+m)),m-1)-np.matlib.repmat(np.arange(0,m-1),math.comb(H2+m-1,m-1),1)-1
            W2 = (np.concatenate((W2,np.zeros([len(W2),1])+H2),axis=1)-np.concatenate((np.zeros([len(W2),1]),W2),axis=1))/H2
            W = np.concatenate((W,W2/2+1/(2*m)),axis=0)
    W[W<1e-6] = 1e-6
    return W
