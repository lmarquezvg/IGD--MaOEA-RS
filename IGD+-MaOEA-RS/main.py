import sys
import math
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from scipy.spatial import distance
from Public.Uniform_Points import uniformReferencePoints
from Public.Genetic_Operators import simulatedBinaryCrossover,polynomialMutation

#Uploads needed functions for a given problem
def uploadBenchmark(problem):
    global parameters, evaluate
    for i in range(1,8):
        if (problem == 'DTLZ'+str(i)):
            from Problems.DTLZ_Problems import parameters,evaluate
    for i in range(1,8):
        if (problem == 'DTLZ'+str(i)+'_MINUS'):
            from Problems.DTLZMinus_Problems import parameters,evaluate
    for i in range(1,9):
        if (problem == 'IMOP'+str(i)):
            from Problems.IMOP_Problems import parameters,evaluate
    return

#Creates a class for population
class population:
    def __init__(Pop, dec, obj):
        Pop.dec = dec
        Pop.obj = obj

#Generates a random population
def randomPopulation(N,n,m,lb,ub,problem):
    decision = lb+(np.random.rand(N,n)*(ub-lb))
    objective = np.zeros([N,m])
    for i in range(0,N):
        objective[i] = evaluate(problem,decision[i],m)
    return population(decision,objective)

#Selects two random parents from a population
def matingSelection(population):
    indexA = np.random.randint(0, len(population.dec))
    while (True):
        indexB = np.random.randint(0, len(population.dec))
        if (indexA != indexB):
            break
    return population.dec[indexA],population.dec[indexB]

#Generates offspring population
def generateOffspring(P,lb,ub,pc,nc,pm,nm,problem):
    N,n = np.shape(P.dec)
    N,m = np.shape(P.obj)
    decision = np.zeros([N,n])
    objective = np.zeros([N,m])
    s = 0
    for i in range(0,N//2):
        parentA,parentB = matingSelection(P)
        offspringA,offspringB = simulatedBinaryCrossover(parentA, parentB, lb, ub, pc, nc)
        offspringA = polynomialMutation(offspringA, lb, ub, pm, nm)
        offspringB = polynomialMutation(offspringB, lb, ub, pm, nm)
        decision[s] = offspringA
        decision[s+1] = offspringB
        s += 2
    if (N%2 == 1):
        parentA,parentB = matingSelection(P)
        offspringA,offspringB = simulatedBinaryCrossover(parentA, parentB, lb, ub, pc, nc)
        offspringA = polynomialMutation(offspringA, lb, ub, pm, nm)
        decision[s] = offspringA
    for i in range(0,N):
        objective[i] = evaluate(problem,decision[i],m)
    return population(decision,objective)

#Divides the population by fronts
def fastNonDominatedSort(evaluation):
    N = len(evaluation)
    S = []
    Fronts = []
    Front1 = []
    eta = np.zeros(N)
    ranks = np.zeros(N)
    for p in range(0,N):
        Sp = []
        for q in range(0,N):
            if (all(evaluation[p] <= evaluation[q]) and any(evaluation[p] < evaluation[q])):
                Sp.append(q)
            elif (all(evaluation[q] <= evaluation[p]) and any(evaluation[q] < evaluation[p])):
                eta[p] += 1
        S.append(Sp)
        if (eta[p] == 0):
            Front1.append(p)
            ranks[p] = 1
    Fronts.append(Front1)
    i = 0
    while (Fronts[i] != []):
        Q = []
        for p in Fronts[i]:
            for q in S[p]:
                eta[q] -= 1
                if (eta[q] == 0):
                    ranks[q] = i+2
                    Q.append(q)
        i += 1
        Fronts.append(Q)
    Fronts.pop(len(Fronts)-1)
    return Fronts

#Returns best fronts and critical front
def findCriticalFront(R,Fronts,N):
    decision = []
    objective = []
    s = 0
    for i in range(0,len(Fronts)):
        s += len(Fronts[i])
        if (s > N):
            Fl = population(R.dec[Fronts[i]],R.obj[Fronts[i]])
            break
        else:
            if (len(decision) == 0):
                decision = np.copy(R.dec[Fronts[i]])
                objective = np.copy(R.obj[Fronts[i]])
            else:
                decision = np.concatenate((decision,R.dec[Fronts[i]]),axis=0)
                objective = np.concatenate((objective,R.obj[Fronts[i]]),axis=0)
    return decision,objective,Fl

#Calculates d+ between two vectors
def dplus(a,z):
    m = len(a)
    s = 0
    for i in range(0,m):
        s += (max(a[i]-z[i],0))**2
    return np.sqrt(s)

#Calculates IGD+ indicator
def IGDplus(A,Z,Memoization):
    N,m = np.shape(A)
    M = len(Z)
    eps = 1e-6
    s = 0
    d = np.zeros([M,N])
    for i in range(0,M):
        first = np.inf
        second = np.inf
        for j in range(0,N):
            d[i,j] = dplus(A[j],Z[i])
            if (d[i,j] < first):
                second = first
                first = d[i,j]
                Memoization[i,1] = Memoization[i,0]
                Memoization[i,0] = [first,j]
            else:
                if (d[i,j] < second and (d[i,j] <= first-eps or d[i,j] >= first+eps)):
                    second = d[i,j]
                    Memoization[i,1] = [second,j]
        s += first
    return s/M

#Returns a vector of IGD+ contributions
def IGDplusDE(A,Z):
    N,m = np.shape(A)
    M = len(Z)
    Memoization = np.empty([M,2],dtype=object)
    total = IGDplus(A,Z,Memoization)
    C = np.zeros(N)
    for i in range(0,N):
        psi = 0
        for j in range(0,M):
            if (Memoization[j,1] == None):
                Memoization[j,1] = Memoization[j,0]
            if (Memoization[j,0][1] == i):
                psi += Memoization[j,1][0]
            else:
                psi += Memoization[j,0][0]
        psi /= M
        C[i] = abs(total-psi)
    return C

#Adjusts location of reference set based on process of AR-MOEA 
def adjustLocation(P,W):
    N,m = np.shape(P)
    M = len(W)
    Pcopy = np.copy(P)
    Wcopy = np.copy(W)
    Pcopy[Pcopy < 1e-6] = 1e-6
    Wcopy[Wcopy < 1e-6] = 1e-6
    Cosine = 1-distance.cdist(Pcopy,Wcopy,'cosine')
    Pnorm = np.sqrt(np.sum(Pcopy**2,axis=1))[:,np.newaxis]
    Wnorm = np.sqrt(np.sum(Wcopy**2,axis=1))[:,np.newaxis]
    d1 = np.matlib.repmat(Pnorm,1,M)*Cosine
    d2 = np.matlib.repmat(Pnorm,1,M)*np.sqrt(1-Cosine**2)
    best = np.argmin(d2,axis=0)
    return Wcopy*np.matlib.repmat(d1[best,np.arange(M)][:,np.newaxis]/Wnorm,1,m)

#Returns population with best individuals
def environmentalSelection(R,W,N,zmin,zmax,flag):
    Fronts = fastNonDominatedSort(R.obj)
    decision,objective,Fl = findCriticalFront(R,Fronts,N)
    if (len(decision) < N):
        zmin = np.min(np.concatenate(([zmin],R.obj[Fronts[0]]),axis=0),axis=0)
        Padj = Fl.obj-zmin
        if (flag != 3):
            Wsca = W*(zmax-zmin)
        else:
            Wsca = np.copy(W)
        Wadj = adjustLocation(Padj,Wsca)
        C = IGDplusDE(Padj,Wadj)
        while (len(decision)+len(Fl.dec) > N):
            Fl.dec = np.delete(Fl.dec,np.argmin(C),axis=0)
            Fl.obj = np.delete(Fl.obj,np.argmin(C),axis=0)
            C = np.delete(C,np.argmin(C))
        if (len(decision) == 0):
            P = population(Fl.dec,Fl.obj)
        else:
            decision = np.concatenate((decision,Fl.dec),axis=0)
            objective = np.concatenate((objective,Fl.obj),axis=0)
            P = population(decision,objective)  
    else:
        P = population(decision,objective)
    return P

#Normalizes population based on process of NSGA-III
def normalization(F):
    N,m = np.shape(F)
    z_min = np.min(F,axis=0)
    Fprime = F-z_min
    Extreme = np.zeros(m,dtype=int)
    w = np.zeros([m,m])+1E-6+np.eye(m)
    for i in range(0,m):
        Extreme[i] = np.argmin(np.max(Fprime/np.matlib.repmat(w[i,:],N,1),axis=1))
    if (np.linalg.cond(Fprime[Extreme]) > 10 or np.isnan(np.linalg.cond(Fprime[Extreme]))):
        a = np.max(Fprime,axis=0)
        a = a[:, np.newaxis]
    else:
        Hyperplane = np.linalg.solve(F[Extreme].T.dot(F[Extreme]), F[Extreme].T.dot(np.ones([m,1])))
        a = 1/Hyperplane
        if (np.any(np.isnan(a)) or np.any(np.isinf(a)) or np.any(a <= 1e-6)):
            a = np.max(Fprime,axis=0)
            a = a[:, np.newaxis]
    Fn = Fprime/np.matlib.repmat(np.transpose(a),N,1)
    return Fn

#Returns the number of associated solutions of each reference point
def associate(F,W):
    M = len(W)
    Fnorm = np.sqrt(np.sum(F**2,axis=1))
    Fnorm = Fnorm[:,np.newaxis]
    Cosine = 1-distance.cdist(F,W,'cosine')
    Distance = np.matlib.repmat(Fnorm,1,M)*np.sqrt(1-Cosine**2)
    pi = np.argmin(Distance,axis=1)
    rho = np.zeros(M)
    for i in range(0,M):
        rho[i] = len(pi[pi==i])
    return rho

#Performs the reference set adaptation method from A-NSGA-III
def adaptive(F,W,N):
    interval = W[0,-1]-W[1,-1]
    Fn = normalization(F)
    m = np.shape(Fn)[1]
    rho = associate(Fn,W)
    Wold = np.array([])
    while (any(rho >= 2) and not(np.array_equal(Wold,W))):
        Wold = np.copy(W)
        for i in numpy.flatnonzero(rho >= 2):
            p = np.matlib.repmat(W[i,:],m,1)-interval/m
            p[np.eye(m,dtype=bool)] += interval
            W = np.concatenate((W,p),axis=0)
        W = W[np.any(W < 0,axis=1) == False]
        u,index = np.unique(np.around(W,4),return_index=True,axis=0)
        W = W[np.sort(index)]
        rho = associate(Fn,W)
    index = np.intersect1d(np.arange(N,len(W)),np.where(rho==0)[0])
    if (len(index) > 0):
        W = np.delete(W,index,axis=0)
    return W

#Performs the reference vector regeneration strategy from RVEA*
def referenceVectorRegeneration(F,W):
    M,m = np.shape(W)
    Fronts = fastNonDominatedSort(F)
    Fnon = F[Fronts[0]]
    Range = np.array([np.min(Fnon,axis=0),np.max(Fnon,axis=0)])
    Fadj = Fnon-Range[0]
    Wsca = W*(Range[1]-Range[0])
    associate = np.argmax(1-distance.cdist(Fadj,Wsca,'cosine'),axis=1)
    inValid = np.setdiff1d(np.arange(0,M),associate)
    W[inValid] = np.random.rand(len(inValid),m)
    return W

#Performs the reference point adaptation method from AR-MOEA
def RefPointAdaption(A,W,zmin,zmax):
    Fronts = fastNonDominatedSort(A)
    Archive = A[Fronts[0]]
    Archive = Archive[np.sort(np.unique(Archive,return_index=True,axis=0)[1])]
    NA = len(Archive)
    NW = len(W)
    zmin = np.min(np.concatenate(([zmin],Archive),axis=0),axis=0)
    if (NA <= 1):
        RefPoint = np.copy(W)
    else:
        Aadj = Archive-zmin
        Wsca = W*(zmax-zmin)
        Wadj = adjustLocation(Aadj,Wsca)
        Distance = distance.cdist(Aadj,Wadj)
        nearestP = np.argmin(Distance,axis=0)
        ContributingS = nearestP[np.unique(nearestP,return_index=True,axis=0)[1]]
        nearestW = np.argmin(Distance,axis=1)
        ValidW = nearestW[ContributingS][np.unique(nearestW[ContributingS],return_index=True,axis=0)[1]]
        
        Choose = np.in1d(np.arange(0,NA),ContributingS)
        Cosine = 1-distance.cdist(Aadj,Aadj,'cosine')
        Cosine[np.eye(len(Cosine),dtype=bool)] = 0
        while (sum(Choose) < min(2*NW,len(Aadj))):
            unSelected = np.where(Choose==False)[0]
            x = np.argmin(np.max(Cosine[Choose==False][:,Choose],axis=1))
            Choose[unSelected[x]] = True
        Archive = Archive[Choose]
        Aadj = Aadj[Choose]
        
        RefPoint = np.concatenate((Wsca[ValidW],Aadj),axis=0)
        Choose = np.concatenate((np.ones(len(ValidW),dtype=bool),np.zeros(len(Aadj),dtype=bool)))
        Cosine = 1-distance.cdist(RefPoint,RefPoint,'cosine')
        Cosine[np.eye(len(Cosine),dtype=bool)] = 0
        while (sum(Choose) < min(NW,len(RefPoint))):
            Selected = np.where(Choose==False)[0]
            x = np.argmin(np.max(Cosine[Choose==False][:,Choose],axis=1))
            Choose[Selected[x]] = True
        RefPoint = RefPoint[Choose]
    return Archive,RefPoint

#Saves approximation set
def saveApproximationSet(P,flag,problem):
    N,m = np.shape(P)
    if (m == 2):
        plt.scatter(P[:,0],P[:,1],c='blue')
        plt.xlabel('$f_1$',rotation=0)
        plt.ylabel('$f_2$',rotation=0)
    elif (m == 3):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(30,45)
        ax.scatter(P[:,0],P[:,1],P[:,2],c='blue')
        ax.xaxis.set_rotate_label(False)
        ax.yaxis.set_rotate_label(False)
        ax.zaxis.set_rotate_label(False)
        ax.set_xlabel('$f_1$')
        ax.set_ylabel('$f_2$')
        ax.set_zlabel('$f_3$')
    else:
        for i in range(0,N):
            plt.plot(P[i])
        plt.xlabel('Objective function',rotation=0)
        plt.ylabel('Objective value')
        x = []
        labels = []
        for i in range(0,m):
            x.append(i)
            labels.append(str(i+1))
        plt.xticks(x,labels)
    if (flag == 0):
        algorithm = 'IGD+-MaOEA-RS'
    if (flag == 1):
        algorithm = 'A-IGD+-MaOEA-RS'
    if (flag == 2):
        algorithm = 'R-IGD+-MaOEA-RS'
    if (flag == 3):
        algorithm = 'AR-IGD+-MaOEA-RS'
    plt.title(algorithm+' on '+problem)
    plt.tight_layout()
    plt.savefig('Results/'+algorithm+'_on_'+problem+'.png')
    plt.close()
    np.savetxt('Results/'+algorithm+'_on_'+problem+'.txt',P,fmt="%.6f")
    return

#Runs main framework
def main(H,flag,problem,m,max_evaluations):
    uploadBenchmark(problem)
    n,m,lb,ub = parameters(problem,m)
    pc = 1   #Crossover probability
    nc = 20  #Crossover distribution index
    pm = 1/n #Mutation probability
    nm = 20  #Mutation distribution index
    N = math.comb(H+m-1,m-1)
    P = randomPopulation(N,n,m,lb,ub,problem)
    W = uniformReferencePoints(N,m)
    zmin = np.min(P.obj,axis=0)
    zmax = np.max(P.obj,axis=0)
    zmax[zmax-zmin <= 1e-6] = 1
    if (flag == 2):
        Wprime = np.random.rand(N,m)
    if (flag == 3):
        A,Wprime = RefPointAdaption(P.obj,W,zmin,zmax)
    evaluations = 0
    while (evaluations < max_evaluations):
        Q = generateOffspring(P,lb,ub,pc,nc,pm,nm,problem)
        R = population(np.concatenate((P.dec,Q.dec),axis=0),np.concatenate((P.obj,Q.obj),axis=0))
        if (flag == 0):
            P = environmentalSelection(R,W,N,zmin,zmax,flag)
        if (flag == 1):
            P = environmentalSelection(R,W,N,zmin,zmax,flag)
            W = adaptive(P.obj,W,N)
        if (flag == 2):
            P = environmentalSelection(R,np.concatenate((W,Wprime),axis=0),N,zmin,zmax,flag)
            Wprime = referenceVectorRegeneration(P.obj,Wprime)
        if (flag == 3):
            A,Wprime = RefPointAdaption(np.concatenate((A,Q.obj),axis=0),W,zmin,zmax)
            P = environmentalSelection(R,Wprime,N,zmin,zmax,flag)
        zmin = np.min(np.concatenate(([zmin],P.obj),axis=0),axis=0)
        zmax = np.max(P.obj,axis=0)
        zmax[zmax-zmin <= 1e-6] = 1
        evaluations += N
        if (flag == 0):
            print('Algorithm: IGD+-MaOEA-RS, 'f'Problem: {problem}, 'f'Evaluations: {evaluations}')
        if (flag == 1):
            print('Algorithm: A-IGD+-MaOEA-RS, 'f'Problem: {problem}, 'f'Evaluations: {evaluations}')
        if (flag == 2):
            print('Algorithm: R-IGD+-MaOEA-RS, 'f'Problem: {problem}, 'f'Evaluations: {evaluations}')
        if (flag == 3):
            print('Algorithm: AR-IGD+-MaOEA-RS, 'f'Problem: {problem}, 'f'Evaluations: {evaluations}')
    return P

if __name__ == '__main__':
    if (len(sys.argv) == 1):
        sys.exit("Incorrect number of arguments. For more information, please use: main.py --help")
    if (str(sys.argv[1]) == '--help'):
        f = open('../README.txt',"r")
        contents = f.read()
        f.close()
        print(contents)
    else:
        if (len(sys.argv) != 6):
            sys.exit("Incorrect number of arguments. For more information, please use: main.py --help")
        H = int(sys.argv[1])
        flag = int(sys.argv[2])
        problem = str(sys.argv[3])
        m = int(sys.argv[4])
        max_evaluations = int(sys.argv[5])
        if (H <= 0):
            sys.exit("Invalid value for the number of divisions per objective function. For more information, please use: main.py --help")
        if not(flag == 0 or flag == 1 or flag == 2 or flag == 3):
            sys.exit("Invalid value for the reference set adaptation method. For more information, please use: main.py --help")
        validProblem = 0
        for i in range(1,8):
            if (problem == 'DTLZ'+str(i)):
                validProblem = 1
        for i in range(1,8):
            if (problem == 'DTLZ'+str(i)+'_MINUS'):
                validProblem = 1
        for i in range(1,9):
            if (problem == 'IMOP'+str(i)):
                validProblem = 1
        if (validProblem == 0):
            sys.exit("Invalid MOP name. For more information, please use: main.py --help")
        if (m <= 1):
            sys.exit("Invalid value for the number of objective functions. For more information, please use: main.py --help")
        if (max_evaluations < 0):
            sys.exit("Invalid value for the maximum number of objective function evaluations. For more information, please use: main.py --help")
        P = main(H,flag,problem,m,max_evaluations)
        saveApproximationSet(P.obj,flag,problem)

