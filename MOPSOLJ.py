# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 21:35:16 2016

@author: lj
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import time,random

def CostFunction(Xin):
    #:::: ZDT1 ::::#
    X=Xin[0]
    n = X.shape[0]
    G  = 1 + (9*(np.sum(X) - X[0])/(n-1)) # or G  = 1 + 9*(np.sum(X[2:n]))/(n-1)
    F1 = X[0]
    F2 = G*(1 - np.sqrt(np.divide(X[0],G)))
    F = np.array([F1, F2])
    return F    
def DetermineDomination(X):
    nX=len(X)
    idx=[]
    for i in range(nX):
        X[i]['IsDominated']=False
    for i in range(nX-1):
        for j in range(i+1,nX):
            if Dominates(X[i]['Cost'],X[j]['Cost'])  :
                X[j]['IsDominated']=True
                idx.append(j)
            elif Dominates(X[j]['Cost'],X[i]['Cost']) :
                X[i]['IsDominated']=True
                idx.append(i)
    return X,idx
def Dominates(X,Y):
    if isinstance(X,dict) and isinstance(Y,dict):
        X = X['Cost']
        Y = Y['Cost']
    else:
        pass
    b = False
    if all(X<=Y) and any(X<Y):
        b = True
    return b
        
        
def CreateGrid(pop,nGrid,alpha):
    C=np.zeros((len(pop),len(pop[0]['Cost'])))
    for i in range(len(pop)):
        C[i,:]=pop[i]['Cost']
    cmax=np.amax(C,axis=0)
    cmin=np.amin(C,axis=0)
    dc=cmax-cmin
    cmin=cmin-alpha*dc
    cmax=cmax+alpha*dc
#    print cmax
    nObj=len(pop[0]['Cost'])
    Grid=[{} for x in range(nObj)]
    for i in range(nObj):
        cj=np.linspace(cmin[i],cmax[i],num=nGrid+1)
        Grid[i]['LB']=np.array([-np.inf])
        Grid[i]['LB']=np.append(Grid[i]['LB'],cj)
        Grid[i]['UB']=cj
        Grid[i]['UB']=np.append(Grid[i]['UB'],np.array([np.inf]))
    return Grid

def FindGridIndex(particle,Grid):
    nObj=len(particle['Cost'])
    nGrid = len(Grid[0]['LB'])
    particle['GridSubIndex'] = np.zeros(nObj)
#    print particle['GridSubIndex'] 
    for j in range(nObj):
        buf=np.zeros_like(Grid[j]['UB'])
        for k in range(len(buf)):
            if particle['Cost'][j]<Grid[j]['UB'][k]:
                buf[k]=1
        particle['GridSubIndex'][j] = np.nonzero(buf)[0][0]
#        print particle['GridSubIndex'][j]
    particle['GridIndex']=particle['GridSubIndex'][0]
    for j in range(nObj-1): 
        particle['GridIndex'] = particle['GridIndex']-1
        particle['GridIndex']  = nGrid * particle['GridIndex'] 
        particle['GridIndex']  = particle['GridIndex'] + particle['GridSubIndex'] [j+1]
#        print particle['GridIndex']

    return particle
    
def SelectLeader(rep,beta):
    ### Grid Index of all repository Members
    GI =np.zeros(len(rep))
    for i in range(len(rep)):
        GI[i] = rep[i]['GridIndex']     
    ##### Occupied Cells
    OC = np.unique(GI)
    ##### Numbers of particles in occupied cells
    N = np.zeros_like(OC) 
    for k in range(len(N)): 
        N[k]=len(np.nonzero(OC[k]==GI))
    ##### Seletion probilities
    P = np.exp(-beta*N)
    P =np.divide(P,np.sum(P))    
    ##### Selected Cell Index
    sci = RouletteWheelSelection(P)    
    ##### Selected Cell Index
    sc = OC[sci]
#    print sc   
    ##### Selected Cell Members
    SCM = np.nonzero(GI==sc)    
    ##### Selected Member Index
    smi = np.random.randint(0,high = len(SCM))
    ##### Selected Member
    sm = SCM[smi][0]
    ##### Leader
    leader = rep[sm]   
    return leader
        
def RouletteWheelSelection(P):

    r = np.random.rand()
    C = np.cumsum(P)
    id = np.nonzero(r<=C)[0]
    i = sorted(set(id))[0]
    return i 
    
    
#def Mutate(x,pm,lb,ub):
#    D = len(x)
#    j = np.random.randint(0,high = D)
#    dx = pm * (ub - lb)
#    lb1 = x[j]-dx[j]  
#    maskl1 = lb1 < lb[j]
#    maskl2 = lb1 >= lb[j]
#    lb1 = maskl1 * lb[j] +maskl2 * lb1
#    ub1 = x[j] + dx[j]
#    masku1 = ub1 > ub[j]
#    masku2 = ub1 <= ub[j]
#    ub1 = masku1 * ub[j] +masku2 * ub1    
#    xnew = x
#    xnew[j] =np.random.uniform(lb1,ub1)
#    return xnew
    
def Mutate(x,pm,lb,ub):
    x = list(x[0])
    mu = 0.02
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
    
#def Mutate(x,p,lb,ub):
#    D = len(x)
#    j = np.random.randint(0,high = D)
#    dx = pm * (ub - lb)
#    lb1 = x[j]-dx[j]  
#    maskl1 = lb1 < lb
#    maskl2 = lb1 >= lb
#    lb1 = maskl1 * lb +maskl2 * lb1
#    ub1 = x[j] + dx[j]
#    masku1 = ub1 > ub
#    masku2 = ub1 <= ub
#    ub1 = masku1 * ub +masku2 * ub1    
#    xnew =lb1 + np.random.rand(1,D)*(ub1 -lb1)
#    return xnew

def DeleteOneRepMember(rep,gamma):
    ####### Grid index of all repository members
    GI =np.zeros(len(rep))
    for i in range(len(rep)):
        GI[i] = rep[i]['GridIndex']     
    ##### Occupied Cells
    OC = np.unique(GI)
    ##### Numbers of particles in occupied cells
    N = np.zeros_like(OC) 
    for k in range(len(N)): 
        N[k]=len(np.nonzero(OC[k]==GI))
    ##### Seletion probilities
    P = np.exp(gamma*N)
    P =np.divide(P,np.sum(P))    
    ##### Selected Cell Index
    sci = RouletteWheelSelection(P)    
    ##### Selected Cell Index
    sc = OC[sci]
#    print sc   
    ##### Selected Cell Members
    SCM = np.nonzero(GI==sc)    
    ##### Selected Member Index
    smi = np.random.randint(0,high = len(SCM))   
    ##### Selected Member
    sm = SCM[smi][0]
    ##### Leader
    del  rep[sm]
    return rep

    
        

           

Nobj = 2
ub = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]) #Lower Bound of Variables#
lb = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) #Upper Bound of Variables
D=len(ub)       # Number of Decision Variables
DSize = [1,D]   # Size of Decision Variables Matrix
maxit=200           #Maximum Number of Iterations
SwarmSize=200               # Swarm size
nRep=200            # Repository Size
w=0.5               # Inertia Weight
wdamp=0.99          # Intertia Weight Damping Rate
c1=2.8                # Personal Learning Coefficient
c2=1.3             # Global Learning Coefficient
nGrid=20            # Number of Grids per Dimension
alpha=0.1           # Inflation Rate
beta=0.1             # Leader Selection Pressure
gamma=2             # Deletion Selection Pressure
mu=0.05             # Mutation Rate


pop=[{} for x in range(SwarmSize)]

for x in pop:
    x['Position'] = lb + np.random.rand(1,D)*(ub - lb)
    x['Velocity'] = np.zeros(D,dtype=float)
    x['Cost'] = CostFunction(x['Position'])
    #update Personal Best
    x['Best']={}
    x['Best']['Position'] = x['Position'].copy()
    x['Best']['Cost']=x['Cost'].copy()
    x['IsDominated'] = False
    x['GridIndex'] = []
    x['GridSubIndex'] = []
pop,idx=DetermineDomination(pop)
#print idx
del_list=sorted(set(idx))
#print del_list
rep=[v for i, v in enumerate(pop) if i not in del_list]
#print rep[1]['Cost']
Grid = CreateGrid(rep,nGrid,alpha)

for i in range(len(rep)):
    rep[i]=FindGridIndex(rep[i],Grid)


#############  MOPSO main loop   ################
plt.ion()
fig = plt.figure()
axes = fig.add_subplot(111)
axes.set_autoscale_on(True) # enable autoscale
axes.autoscale_view(True,True,True)
l, = plt.plot([x['Cost'][0] for x in rep],[x['Cost'][1] for x in rep], 'ro') # Plot blank data
plt.xlabel('x')         # Set up axes
plt.title('test')



for it in range(maxit):
    for i in range(SwarmSize):
        w = 0.9 - ((0.9 - 0.4)/maxit)*(it+1)
        leader = SelectLeader(rep,beta)
        pop[i]['Velocity'] = w * pop[i]['Velocity']+c1*np.random.rand(1,D)*(pop[i]['Best']['Position']-pop[i]['Position'])+c2*np.random.rand(1,D)*(leader['Position']-pop[i]['Position'])
        pop[i]['Position'] = pop[i]['Position']+pop[i]['Velocity']
        pop[i]['Position'] = np.maximum(pop[i]['Position'],lb)
        pop[i]['Position'] = np.minimum(pop[i]['Position'],ub)
        pop[i]['Cost'] = CostFunction(pop[i]['Position'])
        
        ######## Apply Mutation  ########
        pm = (1-(it-1)/(maxit-1))**(1/mu)
        if np.random.rand() <pm :
            NewSol = pop[i]
            NewSol['Position'] = Mutate(pop[i]['Position'],pm,lb,ub)
            if Dominates(NewSol,pop[i]):
                pop[i]['Position'] = NewSol['Position']
                pop[i]['Cost'] = NewSol['Cost']
            elif Dominates(pop[i],NewSol):
                pass
            elif np.random.rand() < 0.5:
                pop[i]['Position'] = NewSol['Position']
                pop[i]['Cost'] = NewSol['Cost']
                
        if Dominates(pop[i],pop[i]['Best']):
            pop[i]['Best']['Position'] = pop[i]['Position']
            pop[i]['Best']['Cost'] = pop[i]['Cost']
        elif Dominates(pop[i]['Best'],pop[i]):
            pass
        elif np.random.rand() < 0.5 :
            pop[i]['Best']['Position'] = pop[i]['Position']
            pop[i]['Best']['Cost'] = pop[i]['Cost']
    
    ####### Add Non-Dominated Particle to Repository  ##########
    
#    pop,idx = DetermineDomination(pop)
#    del_list = sorted(set(idx))
#    rep = rep + [v for i, v in enumerate(pop) if i not in del_list]
    rep = rep + [v for v in pop]
    
    rep,idx = DetermineDomination(rep)
    del_list = sorted(set(idx))
    rep = [v for i, v in enumerate(rep) if i not in del_list]  
    Grid = CreateGrid(rep,nGrid,alpha)
    for i in range(len(rep)):
        rep[i]=FindGridIndex(rep[i],Grid)
    ####### Check if repository is Full
    
    if len(rep) > nRep:
        Extra = len(rep) - nRep
        for e in range(Extra):
            rep = DeleteOneRepMember(rep,gamma)
    
#    w = w*wdamp
    print 'the '+str(it)+'iteration is finished'
    print len(rep)
    l.set_data([x['Cost'][0] for x in rep],[x['Cost'][1] for x in rep])
    axes.relim()
    axes.autoscale_view(True,True,True)
    plt.draw()
    plt.pause(0.1)

#plt.plot([x['Cost'][0] for x in rep],[x['Cost'][1] for x in rep], 'ro')
#plt.xlabel('$f_{1}$')
#plt.ylabel('$f_{2}$')  
