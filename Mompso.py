# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 21:35:16 2016

@author: lj
"""


import numpy as np


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
            if all(X[i]['Cost']<=X[j]['Cost']) and any(X[i]['Cost']<X[j]['Cost']) :
                X[j]['IsDominated']=True
                idx.append(j)
            elif all(X[j]['Cost']<=X[i]['Cost']) and any(X[j]['Cost']<X[i]['Cost']) :
                X[i]['IsDominated']=True
                idx.append(i)
    return X,idx
                
def CreateGrid(pop,nGrid,alpha):
    C=np.zeros((len(pop),len(pop[0]['Cost'])))
    for i in range(len(pop)):
        C[i,:]=pop[i]['Cost']
    cmax=np.amax(C,axis=0)
    cmin=np.amin(C,axis=0)
    dc=cmax-cmin
    cmin=cmin-alpha*dc
    cmax=cmax+alpha*dc
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
    

        

           

Nobj = 2
ub = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]) #Lower Bound of Variables#
lb = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) #Upper Bound of Variables
D=len(ub)       # Number of Decision Variables
DSize = [1,D]   # Size of Decision Variables Matrix
maxit=200           #Maximum Number of Iterations
SwarmSize=100               # Swarm size
nRep=100            # Repository Size
w=0.5               # Inertia Weight
wdamp=0.99          # Intertia Weight Damping Rate
c1=1                # Personal Learning Coefficient
c2=2                # Global Learning Coefficient
nGrid=20            # Number of Grids per Dimension
alpha=0.1           # Inflation Rate
beta=2              # Leader Selection Pressure
gamma=2             # Deletion Selection Pressure
mu=0.1              # Mutation Rate

#empty_particle = {}
#empty_particle['Position'] = np.array([])
#empty_particle['Velocity'] = np.array([])
#empty_particle['Cost'] = np.array([])
#empty_particle['Best'] = {}
#empty_particle['Best']['Postion'] = np.zeros(D)
#empty_particle['Best']['Cost'] = np.zeros(Nobj)

pop=[{} for x in range(SwarmSize)]

for x in pop:
    x['Position'] = lb + np.random.rand(1,D)*(ub - lb)
    x['Velocity'] = np.zeros(SwarmSize,dtype=float)
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

#print rep[i]