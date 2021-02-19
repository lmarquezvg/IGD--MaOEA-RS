import numpy as np

#Returns number of decision variables and their lower and upper bounds
def parameters(problem,m):
    if ((problem == 'IMOP1') or (problem == 'IMOP2') or (problem == 'IMOP3')):
        m = 2
    else:
        m = 3
    K = 5
    L = 5
    n = K+L
    lb = np.zeros(n)
    ub = np.ones(n)
    return n,m,lb,ub

#Evaluates an individual for a given problem
def evaluate(problem,individual,m):
    K = 5
    L = 5
    a1 = 0.05
    a2 = 0.05
    a3 = 10
    if (problem == 'IMOP1'):
        return IMOP1(individual,K,L,a1)
    if (problem == 'IMOP2'):
        return IMOP2(individual,K,L,a1)
    if (problem == 'IMOP3'):
        return IMOP3(individual,K,L,a1)
    if (problem == 'IMOP4'):
        return IMOP4(individual,K,L,a1)
    if (problem == 'IMOP5'):
        return IMOP5(individual,K,L,a2,a3)
    if (problem == 'IMOP6'):
        return IMOP6(individual,K,L,a2,a3)
    if (problem == 'IMOP7'):
        return IMOP7(individual,K,L,a2,a3)
    if (problem == 'IMOP8'):
        return IMOP8(individual,K,L,a2,a3)

#Defines function y1
def y1Function(individual,K,a1):
    s = 0
    for i in range(0,K):
        s += individual[i]
    return (s/K)**a1

#Defines function y2
def y2Function(individual,K,a2):
    s = 0
    n = 0
    for i in range(0,K,2):
        s += individual[i]
        n += 1
    return (s/n)**a2

#Defines function y3
def y3Function(individual,K,a3):
    s = 0
    n = 0
    for i in range(1,K,2):
        s += individual[i]
        n += 1
    return (s/n)**a3

#Defines function g
def gFunction(individual,K,L):
    s = 0
    for i in range(K,K+L):
        s += (individual[i]-0.5)**2
    return s

#Evaluates an individual for the IMOP1 problem
def IMOP1(individual,K,L,a1):
    evaluation = np.zeros(2)
    y1 = y1Function(individual,K,a1)
    g = gFunction(individual,K,L)
    evaluation[0] = g+(np.cos(np.pi/2*y1)**8)
    evaluation[1] = g+(np.sin(np.pi/2*y1)**8)
    return evaluation

#Evaluates an individual for the IMOP2 problem
def IMOP2(individual,K,L,a1):
    evaluation = np.zeros(2)
    y1 = y1Function(individual,K,a1)
    g = gFunction(individual,K,L)
    evaluation[0] = g+(np.cos(np.pi/2*y1)**0.5)
    evaluation[1] = g+(np.sin(np.pi/2*y1)**0.5)
    return evaluation

#Evaluates an individual for the IMOP3 problem
def IMOP3(individual,K,L,a1):
    evaluation = np.zeros(2)
    y1 = y1Function(individual,K,a1)
    g = gFunction(individual,K,L)
    evaluation[0] = g+1+1/5*np.cos(10*np.pi*y1)-y1
    evaluation[1] = g+y1
    return evaluation

#Evaluates an individual for the IMOP4 problem
def IMOP4(individual,K,L,a1):
    evaluation = np.zeros(3)
    y1 = y1Function(individual,K,a1)
    g = gFunction(individual,K,L)
    evaluation[0] = (1+g)*y1
    evaluation[1] = (1+g)*(y1+1/10*np.sin(10*np.pi*y1))
    evaluation[2] = (1+g)*(1-y1)
    return evaluation

#Evaluates an individual for the IMOP5 problem
def IMOP5(individual,K,L,a2,a3):
    evaluation = np.zeros(3)
    y2 = y2Function(individual,K,a2)
    y3 = y3Function(individual,K,a3)
    g = gFunction(individual,K,L)
    h1 = 0.4*np.cos(np.pi/4*np.ceil(8*y2))+0.1*y3*np.cos(16*np.pi*y2)
    h2 = 0.4*np.sin(np.pi/4*np.ceil(8*y2))+0.1*y3*np.sin(16*np.pi*y2)
    evaluation[0] = g+h1
    evaluation[1] = g+h2
    evaluation[2] = g+0.5-h1-h2
    return evaluation

#Evaluates an individual for the IMOP6 problem
def IMOP6(individual,K,L,a2,a3):
    evaluation = np.zeros(3)
    y2 = y2Function(individual,K,a2)
    y3 = y3Function(individual,K,a3)
    g = gFunction(individual,K,L)
    r = max(0,min(np.sin(3*np.pi*y2)**2,np.sin(3*np.pi*y3)**2)-0.05)
    evaluation[0] = (1+g)*y2+np.ceil(r)
    evaluation[1] = (1+g)*y3+np.ceil(r)
    evaluation[2] = (0.5+g)*(2-y2-y3)+np.ceil(r)
    return evaluation

#Evaluates an individual for the IMOP7 problem
def IMOP7(individual,K,L,a2,a3):
    evaluation = np.zeros(3)
    y2 = y2Function(individual,K,a2)
    y3 = y3Function(individual,K,a3)
    g = gFunction(individual,K,L)
    h1 = (1+g)*np.cos(np.pi/2*y2)*np.cos(np.pi/2*y3)
    h2 = (1+g)*np.cos(np.pi/2*y2)*np.sin(np.pi/2*y3)
    h3 = (1+g)*np.sin(np.pi/2*y2)
    r = min(min(abs(h1-h2),abs(h2-h3)),abs(h3-h1))
    evaluation[0] = h1+10*max(0,r-0.1)
    evaluation[1] = h2+10*max(0,r-0.1)
    evaluation[2] = h3+10*max(0,r-0.1)
    return evaluation

#Evaluates an individual for the IMOP8 problem
def IMOP8(individual,K,L,a2,a3):
    evaluation = np.zeros(3)
    y2 = y2Function(individual,K,a2)
    y3 = y3Function(individual,K,a3)
    g = gFunction(individual,K,L)
    evaluation[0] = y2
    evaluation[1] = y3
    evaluation[2] = (1+g)*(3-(y2*(1+np.sin(19*np.pi*y2))/(1+g))-(y3*(1+np.sin(19*np.pi*y3))/(1+g)))
    return evaluation
