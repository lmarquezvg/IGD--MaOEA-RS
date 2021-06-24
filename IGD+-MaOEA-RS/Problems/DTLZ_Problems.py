"""
DTLZ test problems.

K. Deb, L. Thiele, M. Laumanns, and E. Zitzler, "Scalable Test Problems for 
Evolutionary Multiobjective Optimization," in Evolutionary Multiobjective 
Optimization. Theoretical Advances and Applications, pp. 105-145, 2005.
"""

import numpy as np

#Returns number of decision variables and their lower and upper bounds
def parameters(problem,m):
    if (problem == 'DTLZ1'):
        k = 5
    elif (problem == 'DTLZ7'):
        k = 20
    else:
        k = 10
    n = m-1+k
    lb = np.zeros(n)
    ub = np.ones(n)
    return n,m,lb,ub

#Evaluates an individual for a given problem
def evaluate(problem,individual,m):
    if (problem == 'DTLZ1'):
        return DTLZ1(individual,m)
    if (problem == 'DTLZ2'):
        return DTLZ2(individual,m)
    if (problem == 'DTLZ3'):
        return DTLZ3(individual,m)
    if (problem == 'DTLZ4'):
        return DTLZ4(individual,m)
    if (problem == 'DTLZ5'):
        return DTLZ5(individual,m)
    if (problem == 'DTLZ6'):
        return DTLZ6(individual,m)
    if (problem == 'DTLZ7'):
        return DTLZ7(individual,m)

#Evaluates an individual for the DTLZ1 problem
def DTLZ1(individual,m):
    s = 0
    n = len(individual)
    evaluation = np.zeros(m)
    for i in range(m-1,n):
        s += ((individual[i]-0.5)**2-np.cos(20*np.pi*(individual[i]-0.5)))
    g = 100*(n-m+1+s)
    for i in range(0,m):
        if (i == 0):
            mult = 1
            for j in range(0,m-1):
                mult *= individual[j]
            evaluation[i] = (1/2)*mult*(1+g)
        else:
            mult = 1
            for j in range(0,m-i-1):
                mult *= individual[j]
            evaluation[i] = (1/2)*mult*(1-individual[m-i-1])*(1+g)
    return evaluation

#Evaluates an individual for the DTLZ2 problem
def DTLZ2(individual,m):
    s = 0
    n = len(individual)
    evaluation = np.zeros(m)
    for i in range(m-1,n):
        s += (individual[i]-0.5)**2
    g = s
    for i in range(0,m):
        if (i == 0):
            mult = 1
            for j in range(0,m-1):
                mult *= np.cos(individual[j]*np.pi/2)
            evaluation[i] = (1+g)*mult
        else:
            mult = 1
            for j in range(0,m-i-1):
                mult *= np.cos(individual[j]*np.pi/2)
            evaluation[i] = (1+g)*mult*np.sin(individual[m-i-1]*np.pi/2)
    return evaluation

#Evaluates an individual for the DTLZ3 problem
def DTLZ3(individual,m):
    s = 0
    n = len(individual)
    evaluation = np.zeros(m)
    for i in range(m-1,n):
        s += ((individual[i]-0.5)**2-np.cos(20*np.pi*(individual[i]-0.5)))
    g = 100*(n-m+1+s)
    for i in range(0,m):
        if (i == 0):
            mult = 1
            for j in range(0,m-1):
                mult *= np.cos(individual[j]*np.pi/2)
            evaluation[i] = (1+g)*mult
        else:
            mult = 1
            for j in range(0,m-i-1):
                mult *= np.cos(individual[j]*np.pi/2)
            evaluation[i] = (1+g)*mult*np.sin(individual[m-i-1]*np.pi/2)
    return evaluation

#Evaluates an individual for the DTLZ4 problem
def DTLZ4(individual,m):
    s = 0
    n = len(individual)
    evaluation = np.zeros(m)
    for i in range(m-1,n):
        s += (individual[i]-0.5)**2
    g = s
    alpha = 100
    for i in range(0,m):
        if (i == 0):
            mult = 1
            for j in range(0,m-1):
                mult *= np.cos((individual[j]**alpha)*np.pi/2)
            evaluation[i] = (1+g)*mult
        else:
            mult = 1
            for j in range(0,m-i-1):
                mult *= np.cos((individual[j]**alpha)*np.pi/2)
            evaluation[i] = (1+g)*mult*np.sin((individual[m-i-1]**alpha)*np.pi/2)
    return evaluation

#Evaluates an individual for the DTLZ5 problem
def DTLZ5(individual,m):
    s = 0
    n = len(individual)
    evaluation = np.zeros(m)
    for i in range(m-1,n):
        s += (individual[i]-0.5)**2
    g = s
    theta = np.zeros(m-1)
    for i in range(0,m-1):
        if (i == 0):
            theta[i] = individual[i]
        else:
            theta[i] = 2/(4*(1+g))*(1+2*g*individual[i])
    for i in range(0,m):
        if (i == 0):
            mult = 1
            for j in range(0,m-1):
                mult *= np.cos(theta[j]*np.pi/2)
            evaluation[i] = (1+g)*mult
        else:
            mult = 1
            for j in range(0,m-i-1):
                mult *= np.cos(theta[j]*np.pi/2)
            evaluation[i] = (1+g)*mult*np.sin(theta[m-i-1]*np.pi/2)
    return evaluation

#Evaluates an individual for the DTLZ6 problem
def DTLZ6(individual,m):
    s = 0
    n = len(individual)
    evaluation = np.zeros(m)
    for i in range(m-1,n):
        s += individual[i]**0.1
    g = s
    theta = np.zeros(m-1)
    for i in range(0,m-1):
        if (i == 0):
            theta[i] = individual[i]
        else:
            theta[i] = 2/(4*(1+g))*(1+2*g*individual[i])
    for i in range(0,m):
        if (i == 0):
            mult = 1
            for j in range(0,m-1):
                mult *= np.cos(theta[j]*np.pi/2)
            evaluation[i] = (1+g)*mult
        else:
            mult = 1
            for j in range(0,m-i-1):
                mult *= np.cos(theta[j]*np.pi/2)
            evaluation[i] = (1+g)*mult*np.sin(theta[m-i-1]*np.pi/2)
    return evaluation

#Evaluates an individual for the DTLZ7 problem
def DTLZ7(individual,m):
    s = 0
    n = len(individual)
    evaluation = np.zeros(m)
    for i in range(m-1,n):
        s += individual[i]
    g = 1+9/(n-m+1)*s
    s = 0
    for i in range(0,m-1):
        evaluation[i] = individual[i]
        s += evaluation[i]/(1+g)*(1+np.sin(3*np.pi*evaluation[i]))
    h = m-s
    evaluation[m-1] = (1+g)*h
    return evaluation
