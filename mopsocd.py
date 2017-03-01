# -*- coding: utf-8 -*-
"""
Created on Thu Feb 23 20:49:46 2017

@author: lj
"""
import numpy as np
import random
from matplotlib.pylab import plt
def Dominates(X,Y):
    b = False
    if all(X<=Y) and any(X<Y):
        b = True
    return b        
def Mutate(x,pm,mu,lb,ub):
    x = list(x)
    mu = 0.05
    nVar = len(x)
    nMu = np.ceil(mu*nVar)
    j = random.sample(range(nVar),int(nMu))

    sigma =  pm * (ub - lb)
    if len(sigma)>1:
        sigma = [sigma[i] for i in j]
    y = x
    rerror= np.random.rand(len(j))*sigma
    for v in range(int(nMu)):
        y[j[v]]=x[j[v]]+rerror[v]
#        print y[v]
        if y[j[v]] < lb[j[v]]:
            y[j[v]]=lb[j[v]]
        if y[j[v]] > ub[j[v]]:
            y[j[v]] = ub[j[v]]
    return np.array(y)
    
           
def nondomSolutions(X, F):
    X = np.array(X)
    F = np.array(F)
    nF = len(F[:,0])
    nD = np.zeros((nF,1))
    index = []
    for i in xrange(0,nF):
        nD[i][0] = 0
        for j in xrange(0,nF):
            if j != i:
                if (all(F[j,:] <= F[i,:]) & any(F[j,:] < F[i,:])):
                    nD[i][0] = nD[i][0] + 1
        if nD[i][0] == 0:
            index.append(i)
    repX = X[index,:]
    repF = F[index,:]
    it = 0
    nMax = len(repX[:,0])
    while (it < nMax):
        uIdx = sum((2 == (0 == (repX - repX[it,:])).sum(1)).astype(int)) 
        if uIdx > 1:
            repX = np.delete(repX, it, axis = 0)
            repF = np.delete(repF, it, axis = 0)
            nMax = len(repX[:,0])
            it = 0
        else:
            it = it + 1
            nMax = len(repX[:,0])     
    return repX, repF
	
def crowding(archiveX,archiveF):
    nondomN=len(archiveF)
    Nobj=len(archiveF[0,:])
    crowdDist=np.zeros((nondomN,1))
    for i in range(Nobj):
        archiveF=archiveF[np.argsort(archiveF[:,i])]
        archiveX=archiveX[np.argsort(archiveF[:,i])]
        crowdDist=crowdDist[np.argsort(archiveF[:,i])]
        for j in range(1,nondomN-1):
            crowdDist[j,0]=crowdDist[j,0]+abs(archiveF[j+1,i]-archiveF[j-1,i])/abs(archiveF[0,i]-archiveF[-1,i])
#        crowdDist[0,0]=crowdDist[0,0]+np.inf
#        crowdDist[nondomN-1,0]=crowdDist[nondomN-1,0]+1
        crowdDist[0,0]=np.inf
        crowdDist[nondomN-1,0]=np.inf
    archiveX=archiveX[np.argsort(crowdDist[:,0])[::-1]]
    archiveF=archiveF[np.argsort(crowdDist[:,0])[::-1]]

    return archiveX,archiveF

def removeparticles(archiveX,archiveF,N):
	for i in range(int(N)):
		archiveX=np.delete(archiveX,-1,axis=0)
		archiveF=np.delete(archiveF,-1,axis=0)
	return archiveX,archiveF
def SelectLeader(archiveX,archiveF):
	Nx=np.floor(len(archiveX)*0.1)
	h=random.randint(0,Nx)
	return h


def mopsocd(costFunction,lb,ub,Nobj,swarmsize=200, maxit=60,nRep=100,w=0.4,wdamp=0.99,c1=2.8,c2=1.4,alpha=0.1,mu=0.05,bait=[]):

    vmax = (ub-lb)*0.1

    D=len(ub)
    X = np.random.rand(swarmsize, D)  # particle positions
    V = np.zeros_like(X)  # particle velocities
    F= np.zeros((swarmsize,Nobj))  # current particle function values
    PXbest = np.zeros_like(X)  # best particle positions
    PFbest = np.ones(swarmsize)*np.inf  # best particle function values
    
    # Initialize the particle's position
    X = lb + X*(ub - lb)
#    X[0,:]=np.array([2.7109,2.7109,2.7109,2.7109,-4.3580,-4.7055,-4.7055,-4.3580,5.3454,5.8778,5.8778,5.3454])
#    Vhigh = np.abs(ub - lb)
#    Vlow = -Vhigh
#    V = Vlow + np.random.rand(S, D)*(Vhigh - Vlow)
    if bait !=[]:
        idxx = np.random.randint(0,high=swarmsize,size=len(bait))
        print idxx,swarmsize
        for idx,elem in zip(idxx,bait):
            
            X[idx,:]=np.array(elem)
            

    processes=1
    if processes > 1:
        import multiprocessing
        mp_pool = multiprocessing.Pool(processes)
        F= np.array(mp_pool.map(costFunction, X))
    else:
        for i in range(swarmsize):
            F[i]= costFunction(X[i, :])
#            print F[i]
    PXbest=X.copy()
    PFbest=F.copy()
    nondomX,nondomF=nondomSolutions(X,F)
    
    archiveX=nondomX.copy() 
    archiveF=nondomF.copy()
    
    if len(nondomX)>2:
        archiveX,archiveF=crowding(archiveX,archiveF)
    
    
    it=0
    while it<maxit:
        it=it+1
        h=SelectLeader(archiveX,archiveF)    
        w=wdamp-(wdamp-w)/maxit*it
        for i in xrange(0, swarmsize):
    	#:::: Rep[h] is a value that is taken from the repository ::::#
    #        h=SelectLeader(archiveX,archiveF)  
            
            V[i,:] = w*V[i,:] + random.random()*c1*(PXbest[i,:] - X[i,:]) + random.random()*c2*(archiveX[h,:] - X[i,:]);
            X[i,:] = X[i,:] + V[i,:]
    
    		#:::: Check that the particles do not fly out of the search space ::::#
            X[i,X[i,:] < lb] = lb[X[i,:] < lb]
            X[i,X[i,:] > ub] = ub[X[i,:] > ub]
    		#:::: Change the velocity direction ::::#
            V[i,X[i,:] < lb] = -V[i,X[i,:] < lb]
            V[i,X[i,:] > ub] = -V[i,X[i,:] > ub]
    		#:::: Constrain the velocity ::::#
            V[i, V[i,:] > vmax] = vmax[V[i,:] > vmax]
            V[i, V[i,:] < -vmax] = vmax[V[i,:] < -vmax]
        for i in xrange(0,swarmsize):
            F[i,:] = costFunction(X[i,:])
            pm=(1-(it-1)/maxit)**(1/mu)
    
            if np.random.rand() <pm :
                NewX = X[i,:].copy()
                NewX = Mutate(X[i],pm,mu,lb,ub)
                NewF = costFunction(NewX)
            else:
                NewX = X[i,:].copy()
                NewF = F[i,:].copy()
            if Dominates(NewF,F[i,:]):
                X[i,:] = NewX.copy()
                F[i,:] = NewF.copy()
            elif Dominates(F[i,:],NewF):
                pass
            elif np.random.rand() < 0.5:
                X[i,:] = NewX.copy()
                F[i,:] = NewF.copy()
            
            if Dominates(F[i,:],PFbest[i,:]):
                PFbest[i,:] = F[i,:].copy()
                PXbest[i,:] = X[i,:].copy()    
            elif Dominates(PFbest[i,:],F[i,:]):
                pass
            elif (random.randint(0,1) == 0):
                PFbest[i,:] = F[i,:].copy()
                PXbest[i,:] = X[i,:].copy()
        archiveX = np.concatenate((archiveX,X),axis = 0)    
        archiveF = np.concatenate((archiveF,F),axis = 0) 
        archiveX,archiveF = nondomSolutions(archiveX,archiveF)
        if len(archiveX)>2:
            archiveX,archiveF=crowding(archiveX,archiveF)
        if len(archiveX)>nRep:
            print len(archiveX)-nRep,'particles removed'
            archiveX,archiveF=removeparticles(archiveX,archiveF,len(archiveX)-nRep)
            
        print 'the',it,'th iteration finished'
    
    
    np.savetxt('archiveF.txt',archiveF)
    np.savetxt('archiveX.txt',archiveX)
    plt.plot(archiveF[:,0], archiveF[:,1],'ro')
    plt.xlabel('$f_{1}$')
    plt.ylabel('$f_{2}$')    
