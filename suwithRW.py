# -*- coding: utf-8 -*-
"""
Created on Sat Dec 31 22:18:42 2016

@author: lj
"""
from __future__ import division
from atelements import *
import numpy as np
from scipy import constants
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import Rbf
import matplotlib.pyplot as plt
#element.E = 0.629
frev=constants.c/48
A=element(629E-3,80*frev)  #ini#
A.set_global()     #ini #



Xin=np.array([ 2.459051012,2.459051012,3.180633754,3.180633754,-5.088073189,-3.724263674,-3.724263674,-5.088073189,5.248167791,4.864659493,4.864659493,5.248167791])
Xin=np.array([2.4746,2.9624,2.9624,2.9624, -4.50667, -4.17499, -4.17499, -4.50667,5.06574,5.20415,5.20415,5.06574])
Sin=np.array([0,0,0])

#A=np.loadtxt('fie12.map')
#
#
#Bx=A[:,5]
#By=A[:,4]
#Bz=A[:,3]
#
#
#Px=np.unique(A[:,2])
#Py=np.unique(A[:,1])
#Pz=np.unique(A[:,0])
#
#Nx=int((np.max(Px)-np.min(Px))/(A[1,2]-A[0,2]))+1
#Ny=int((np.max(Py)-np.min(Py))/(Py[1]-Py[0]))+1
#Nz=int((np.max(Pz)-np.min(Pz))/(Pz[1]-Pz[0]))+1
#
#Px=np.linspace(np.min(Px),np.max(Px),Nx)
#
#U = np.zeros((Nz,Ny,Nx))
#V = np.zeros((Nz,Ny,Nx))
#W = np.zeros((Nz,Ny,Nx))
#tt=0
#for i in range(Nz):
#    for j in range(Ny):
#        for k in range(Nx):
#            U[i,j,k] = Bx[tt]
#            V[i,j,k] = By[tt]
#            W[i,j,k] = Bz[tt]
#            tt=tt+1
#
#f1=RegularGridInterpolator((Pz, Py, Px), U)
#f2=RegularGridInterpolator((Pz, Py, Px), V)
#f3=RegularGridInterpolator((Pz, Py, Px), W)


#RW=Fieldmap('RW',2.3)

RW=Drift('RW',2.3,Nslices=23)
Q1P2K1=Xin[0]
Q1P1L2=Xin[1]
Q1P2L2=Xin[2]
Q1P1K3=Xin[3]
Q1P2K3=Xin[3]
Q1P1L4=Xin[2]
Q1P2L4=Xin[1]
Q1P1K1=Xin[0]

Q2P2K1=Xin[4]
Q2P1L2=Xin[5]
Q2P2L2=Xin[6]
Q2P1K3=Xin[7]
Q2P2K3=Xin[7]
Q2P1L4=Xin[6]
Q2P2L4=Xin[5]
Q2P1K1=Xin[4]   


Q3P2K1=Xin[8] 
Q3P1L2=Xin[9] 
Q3P2L2=Xin[10] 
Q3P1K3=Xin[11]  
Q3P2K3=Xin[11]  
Q3P1L4=Xin[10] 
Q3P2L4=Xin[9] 
Q3P1K1=Xin[8]   
rho =1.2*4/3.1415926
theta = 3.1415926/4

S3P2K1=Sin[2]
S2P2K1=Sin[1]
S1P2K1=Sin[0]
S1P1L2=Sin[0]
S2P1L2=Sin[1]
S3P1L2=Sin[2]

S3P2L2=Sin[2]
S2P2L2=Sin[1]
S1P2L2=Sin[0]

S1P1K3=Sin[0]
S2P1K3=Sin[1]
S3P1K3=Sin[2]

S3P2K3=Sin[2]
S3P2K3=Sin[1]
S3P2K3=Sin[0]
S1P1L4= Sin[0]
S2P1L4= Sin[1]
S3P1L4= Sin[2]

S3P2L4=Sin[2]
S2P2L4=Sin[1]
S1P2L4=Sin[0]
S1P1K1= Sin[0]
S2P1K1= Sin[1]
S3P1K1= Sin[2]

MLS_start= Drift('D0',0)
L1=Drift('D1',1.25,Nslices=25) #drift 1.25m#
S3M2K1=Sextupole('S3P2K1',S3P2K1,0.1,Nslices=4)
L101=Drift('D101',0.15,Nslices=3) #drift 1.25m#
Q3M2K1=Quadrupole('Q3P2K1',Q3P2K1,0.2,Nslices=10)
L2=Drift('D2',0.15,Nslices =3)
Q2M2K1=Quadrupole('Q2P2K1',Q2P2K1,0.2,Nslices=10)
L3= Drift('D3',0.425,Nslices=17)
BB11=Bedge('Bedge11',1/rho,theta/2,hgap=0.025,fint=0.5)
BB10=Sbend('Bend1',rho,theta,Nslices=48)
BB12=Bedge('Bedge12',1/rho,theta/2,hgap=0.025,fint=0.5)
L4=Drift('D4',0.425,Nslices=17)
S2M2K1=Sextupole('S2P2K1',S2P2K1,0.1,Nslices=4)
L401=Drift('D401',0.3,Nslices=12)
S1M2K1=Sextupole('S1P2K1',S1P2K1,0.1,Nslices=4)
L402=Drift('D402',0.15,Nslices=6)
Q1M2K1=Quadrupole('Q1P2K1',Q1P2K1,0.2,Nslices=10) 
L5=Drift('D5',0.35,Nslices=7)
Q1M1L2=Quadrupole('Q1P1L2',Q1P1L2,0.2,Nslices=10) 
L6 = Drift('D6',0.15,Nslices=6)
S1M1L2=Sextupole('S1P1L2',S1P1L2,0.1,Nslices=4)
L601 = Drift('D601',0.3,Nslices=12)
S2M1L2=Sextupole('S2P1L2',S2P1L2,0.1,Nslices=4)
L602 = Drift('D602',0.425,Nslices=17)
BB21=Bedge('Bedge21',1/rho,theta/2,hgap=0.025,fint=0.5)
BB20=Sbend('Bend2',rho,theta,Nslices=48)
BB22=Bedge('Bedge22',1/rho,theta/2,hgap=0.025,fint=0.5)
L7=Drift('D7',0.425,Nslices=85)
Q2M1L2=Quadrupole('Q2P1L2',Q2P1L2,0.2,Nslices=10) 
L8=Drift('D8',0.15,Nslices=3)
Q3M1L2= Quadrupole('Q3P1L2',Q3P1L2,0.2,Nslices=10)
L9=Drift('D9',0.15,Nslices=3)
S3M1L2=Sextupole('S3P1L2',S3P1L2,0.1,Nslices=4)
L901=Drift('D901',3,Nslices=60)
section1=[MLS_start,L1,S3M2K1,L101,Q3M2K1,L2,Q2M2K1,L3,BB11,BB10,BB12,L4,S2M2K1,L401,S1M2K1,L402,Q1M2K1,L5,Q1M1L2,L6,S1M1L2,L601,S2M1L2,L602,\
          BB21,BB20,BB22,L7,Q2M1L2,L8,Q3M1L2,L9,S3M1L2,L901]

L10=Drift('D10',3,Nslices=60)
S3M2L2=Sextupole('S3P2L2',S3P2L2,0.1,Nslices=4)
L1001=Drift('D1001',0.15,Nslices=3)
Q3M2L2=Quadrupole('Q3P2L2',Q3P2L2,0.2,Nslices=10)
L11=Drift('D11',0.15,Nslices=3)
Q2M2L2=Quadrupole('Q2P2L2',Q2P2L2,0.2,Nslices=10)
L12=Drift('D12',0.425,Nslices=85)
BB31=Bedge('Bedge31',1/rho,theta/2,hgap=0.025,fint=0.5)
BB30=Sbend('Bend3',rho,theta,Nslices=48)
BB32=Bedge('Bedge32',1/rho,theta/2,hgap=0.025,fint=0.5)
L13=Drift('D13',0.425,Nslices=85)
S2M2L2=Sextupole('S2P2L2',S2P2L2,0.1,Nslices=4)
L1301=Drift('D1301',0.3,Nslices=60)
S1M2L2=Sextupole('S1M2L2',S1P2L2,0.1,Nslices=4)
L1302=Drift('D1302',0.15,Nslices=85)
Q1M2L2=Quadrupole('Q1P2L2',Q1P2L2,0.2,Nslices=10)
L14=Drift('D14',0.35,Nslices=7)
Q1M1K3=Quadrupole('Q1P1K3',Q1P1K3,0.2,Nslices=10)
L15=Drift('D15',0.15,Nslices=3)
S1M1K3=Sextupole('S1P1K3',S1P1K3,0.1,Nslices=4)
L1501=Drift('D1501',0.3,Nslices=6)
S2M1K3=Sextupole('S2P1K3',S2P1K3,0.1,Nslices=4)
L1502=Drift('D1502',0.425,Nslices=17)
BB41=Bedge('Bedge41',1/rho,theta/2,hgap=0.025,fint=0.5)
BB40=Sbend('Bend4',rho,theta,Nslices=48)
BB42=Bedge('Bedge42',1/rho,theta/2,hgap=0.025,fint=0.5)
L16=Drift('D16',0.425,Nslices=85)
Q2M1K3=Quadrupole('Q2P1K3',Q2P1K3,0.2,Nslices=10)
L17=Drift('D17',0.15,Nslices=3)
Q3M1K3=Quadrupole('Q3P1K3',Q3P1K3,0.2,Nslices=10)
L18=Drift('D18',0.15,Nslices=3)
S3M1K3=Sextupole('S3P1K3',S3P1K3,0.1,Nslices=4)
L1801=Drift('D1801',0.1,Nslices=2)


section2=[L10,S3M2L2,L1001,Q3M2L2,L11,Q2M2L2,L12,BB31,BB30,BB32,L13,S2M2L2,L1301,S1M2L2,L1302,Q1M2L2,L14,Q1M1K3,L15,S1M1K3,L1501,S2M1K3,L1502,\
          BB41,BB40,BB42,L16,Q2M1K3,L17,Q3M1K3,L18,S3M1K3,L1801,RW]

L19=Drift('D19',0.1,Nslices=2) 
S3M2K3=Sextupole('S3P1K3',S3P1K3,0.1,Nslices=4)
L1901=Drift('D1901',0.15,Nslices=3) 
Q3M2K3=Quadrupole('Q3P2K3',Q3P2K3,0.2,Nslices=30)
L20=Drift('D20',0.15,Nslices=3)
Q2M2K3=Quadrupole('Q2P2K3',Q2P2K3,0.2,Nslices=10)
L21=Drift('D21',0.425,Nslices=17)
BB51=Bedge('Bedge51',1/rho,theta/2,hgap=0.025,fint=0.5)
BB50=Sbend('Bend5',rho,theta,Nslices=48)
BB52=Bedge('Bedge52',1/rho,theta/2,hgap=0.025,fint=0.5)

L22=Drift('D22',0.425,Nslices=17)
S2M2K3=Sextupole('S2P1K3',S2P1K3,0.1,Nslices=4)
L2201=Drift('D2201',0.3,Nslices=6)
S1M2K3=Sextupole('S1P1K3',S1P1K3,0.1,Nslices=4)
L2202=Drift('D2202',0.15,Nslices=15)
Q1M2K3=Quadrupole('Q1P2K3',Q1P2K3,0.2,Nslices=10)
L23=Drift('D23',0.35,Nslices=7)
Q1M1L4=Quadrupole('Q1P1L4',Q1P1L4,0.2,Nslices=10)
L24=Drift('D24',0.15,Nslices=3)
S1M1L4=Sextupole('S1P1L4',S1P1L4,0.1,Nslices=4)
L2401=Drift('D2401',0.3,Nslices=6)
S2M1L4=Sextupole('S2P1L4',S2P1L4,0.1,Nslices=4)
L2402=Drift('D2402',0.425,Nslices=17)
BB61=Bedge('Bedge61',1/rho,theta/2,hgap=0.025,fint=0.5)
BB60=Sbend('Bend6',rho,theta,Nslices=48)
BB62=Bedge('Bedge62',1/rho,theta/2,hgap=0.025,fint=0.5)

L25=Drift('D25',0.425,Nslices=85)
Q2M1L4=Quadrupole('Q2P1L4',Q2P1L4,0.2,Nslices=10)
L26=Drift('D26',0.15,Nslices=3)
Q3M1L4=Quadrupole('Q3P1L4',Q3P1L4,0.2,Nslices=10)
L27=Drift('D27',3,Nslices=60)
S3M1L4=Sextupole('S3P1L4',S3P1L4,0.1,Nslices=4)
L2701=Drift('D2701',0.15,Nslices=5)
section3=[L19,S3M2K3,L1901,Q3M2K3,L20,Q2M2K3,L21,BB51,BB50,BB52,L22,S2M2K3,L2201,S1M2K3,L2202,Q1M2K3,L23,Q1M1L4,L24,S1M1L4,L2401,S2M1L4,L2402,BB61,BB60,BB62,L25,Q2M1L4,L26,Q3M1L4,L27,S3M1L4,L2701]

L28=Drift('D28',3,Nslices=60)
S3M2L4=Sextupole('S3P2L4',S3P2L4,0.1,Nslices=4)
L2801=Drift('D2801',0.15,Nslices=5)
Q3M2L4=Quadrupole('Q3P2L4',Q3P2L4,0.2,Nslices=10)
L29=Drift('D29',0.15,Nslices=3)
Q2M2L4=Quadrupole('Q2P2L4',Q2P2L4,0.2,Nslices=10)
L30=Drift('D30',0.425,Nslices=85)
BB71=Bedge('Bedge71',1/rho,theta/2,hgap=0.025,fint=0.5)
BB70=Sbend('Bend7',rho,theta,Nslices=48)
BB72=Bedge('Bedge71',1/rho,theta/2,hgap=0.025,fint=0.5)
L31=Drift('D31',0.425,Nslices=215)
S2M2L4=Sextupole('S2P2L4',S2P2L4,0.1,Nslices=4)
L3101=Drift('D3101',0.3,Nslices=215)
S1M2L4=Sextupole('S1P2L4',S1P2L4,0.1,Nslices=4)
L3102=Drift('D3102',0.15,Nslices=215)
Q1M2L4=Quadrupole('Q1P2L4',Q1P2L4,0.2,Nslices=10)
L32=Drift('D32',0.35,Nslices=7)
Q1M1K1=Quadrupole('Q1P1K1',Q1P1K1,0.2,Nslices=10)
L33=Drift('D33',0.15,Nslices=3)
S1M1K1=Sextupole('S1P1K1',S1P1K1,0.1,Nslices=4)
L3301=Drift('D3301',0.3,Nslices=6)
S2M1K1=Sextupole('S2P1K1',S2P1K1,0.1,Nslices=4)
L3302=Drift('D3302',0.425,Nslices=17)
BB81=Bedge('Bedge81',1/rho,theta/2,hgap=0.025,fint=0.5)
BB80=Sbend('Bend8',rho,theta,Nslices=48) 
BB82=Bedge('Bedge81',1/rho,theta/2,hgap=0.025,fint=0.5)
L34=Drift('D34',0.425,Nslices=85)
Q2M1K1=Quadrupole('Q2P1K1',Q2P1K1,0.2,Nslices=10)
L35=Drift('D35',0.15,Nslices=3)
Q3M1K1=Quadrupole('Q3P1K1',Q3P1K1,0.2,Nslices=10)
L36=Drift('D36',1.25,Nslices=25)
S3M1K1=Sextupole('S3P1K1',S3P1K1,0.1,Nslices=4)
L3601=Drift('D3601',0.15,Nslices=3)
section4=[L28,S3M2L4,L2801,Q3M2L4,L29,Q2M2L4,L30,BB71,BB70,BB72,L31,S2M2L4,L3101,S1M2L4,L3102,Q1M2L4,L32,Q1M1K1,L33,S1M1K1,L3301,S2M1K1,L3302,BB81,BB80,BB82,L34,Q2M1K1,L35,Q3M1K1,L36,S3M1K1,L3601]

seq=section1+section2+section3+section4



ringtwisspara,[Qx,Qy]=Ringtwiss(seq)

print Qx,Qy
U0,emit,sigmaE,ap,[I1,I2,I3,I4,I5],[Jx,Jy,Js]=Cal_RadiationIntegrals(ringtwisspara)
print U0,emit,sigmaE,ap,[I1,I2,I3,I4,I5],[Jx,Jy,Js]
#Beamsize(ringtwisspara,U0,emit,sigmaE,ap,315E-6,0.01)
#
plt.plot(ringtwisspara[:,0],ringtwisspara[:,4],'b',ringtwisspara[:,0],ringtwisspara[:,7],'g',ringtwisspara[:,0],ringtwisspara[:,10]*10,'r')