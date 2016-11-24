# -*- coding: utf-8 -*-
"""
Created on Wed Nov 02 22:36:29 2016

@author: local_admin
"""
from __future__ import division
import numpy as np
from scipy import constants
import random
from scipy import integrate
from matplotlib.pylab import plt

from scipy.optimize import fsolve,broyden1

E = 629E-3
Emass = constants.physical_constants['electron mass energy equivalent in MeV'][0]*1E-3
gamma= E/Emass
beta = np.sqrt(1-1/gamma**2)
couplingfactor=0.01
h=80
Isingle=2
theta =constants.pi/4
rho =1.2/theta
lbend=1.2
Cr=8.846E-5 ##unit m/GeV3
Cq=3.832E-13
def Fquad(kk,L):
    A=np.zeros((6,6))
    K=np.sqrt(kk)
    a = K * L
    A[0,0]= np.cos(a)
    A[0,1] = np.sin(a)/K
    A[1,0] = -K*np.sin(a)
    A[1,1] = np.cos(a)
    A[2,2] = np.cosh(a)
    A[2,3] = np.sinh(a)/K
    A[3,2] = np.sinh(a)*K
    A[3,3] = np.cosh(a)
    A[4,4] = 1
    A[4,5] = L/gamma**2/beta**2
    A[5,5] = 1
    return [A,L,111]
def Dquad(kk,L):
    A=np.zeros((6,6))
    K=np.sqrt(-kk)
    a = K * L
    A[0,0]= np.cosh(a)
    A[0,1] = np.sinh(a)/K
    A[1,0] = np.sinh(a)*K
    A[1,1] = np.cosh(a)
    A[2,2] = np.cos(a)
    A[2,3] = np.sin(a)/K
    A[3,2] = -K*np.sin(a)
    A[3,3] = np.cos(a)
    A[4,4] = 1
    A[4,5] = L/gamma**2/beta**2
    A[5,5] = 1
    return [A,L,112]
def Drift(L):
    A=np.identity(6)
    A[0,1]= L
    A[2,3]= L
    A[4,5] = L/gamma**2/beta**2   
    return [A,L,0]
def Bedge(h,e1,hgap,fint):
    A=np.identity(6)
    psi =2*h*hgap*fint*(1+np.sin(e1)*np.sin(e1))
    A[1,0]= np.tan(e1)*h
    A[3,2] =-np.tan(e1-psi)*h
    L=0
    return [A,L,220]

def Sbend(rho,theta):
    A=np.identity(6)
    A[0,0]= np.cos(theta)
    A[0,1] = rho*np.sin(theta)
    A[0,5] = rho*(1-np.cos(theta))/beta
    A[1,0] = -np.sin(theta)/rho
    A[1,1] = np.cos(theta)
    A[1,5] = np.sin(theta)/beta
    A[2,2] = 1
    A[2,3] = rho*theta
    A[3,3] = 1
    A[4,0] = -np.sin(theta)/beta
    A[4,1] = -rho*(1-np.cos(theta))/beta
    A[4,4] = 1
    A[4,5] = rho*(np.sin(theta)-beta**2*theta)
    A[5,5] = 1
    return [A,rho*theta,221]
def Twiss(T):
    M = np.zeros((3,3))
    M[0,0] = T[0,0]**2
    M[0,1] = -2*T[0,0]*T[0,1]
    M[0,2] = T[0,1]**2
    M[1,0] = -T[0,0]*T[1,0]
    M[1,1] = T[0,0]*T[1,1]+T[0,1]*T[1,0]
    M[1,2] = -T[0,1]*T[1,1]
    M[2,0] = T[1,0]**2
    M[2,1] = -2*T[1,0]*T[1,1]
    M[2,2] = T[1,1]**2
    return M
def CheckStability(Mcell):
    if (np.abs(Mcell[0,0]+Mcell[1,1])>=2) or (np.abs(Mcell[2,2]+Mcell[3,3])>=2):
#        print 'Mcell is unstable.\n'
        return 0
    else:
        return 1
def Cal_Tunes(Mcell):
    mucellh = np.arccos((Mcell[0,0]+Mcell[1,1])/2)
    mucellv = np.arccos((Mcell[2,2]+Mcell[3,3])/2)
    Qh = mucellh/constants.pi/2
    Qv = mucellv/constants.pi/2
    return Qh, Qv,mucellh,mucellv
def Cal_Q(Mcell,s,betax,betay):
    Qx=0
    Qy=0
    for i in range(len(s)):
        if i>0:
            Qx= Qx + 1/betax[i]*(s[i]-s[i-1])
            Qy= Qy + 1/betay[i]*(s[i]-s[i-1])   
    mucellh = np.arccos((Mcell[0,0]+Mcell[1,1])/2)
    mucellv = np.arccos((Mcell[2,2]+Mcell[3,3])/2)
    Qh = mucellh/constants.pi/2+np.floor(Qx/2/constants.pi)
    Qv = mucellv/constants.pi/2+np.floor(Qy/2/constants.pi)
    return Qh,Qv
    
def Cal_alpha():
    alpha=0
    for x in twisspara:
        if x[-1]==220 or x[-1]==221:
            alpha=alpha + x[3]*theta/49
    alpha=alpha/twisspara[-1,0]
    return alpha
def Cal_emittance(twisspara):
    emit=0
    cnt=0
    for x in twisspara:
        if x[-1]==220 or x[-1]==221:
            emit=emit + x[4]
            cnt=cnt+1
    emit=emit/cnt*3.832E-13/rho*gamma**2
    return emit*1E9
def Cal_RadiationIntegrals(twisspara):
    I1=0
    I2=0
    I3=0
    I4=0
    I5=0
    cnt=0
    for i in range(len(twisspara)):
        x=twisspara[i,:]
        
        if x[-1]==220 or x[-1]==221:
            y=twisspara[i+1,:]
            if y[-1]==220 or y[-1]==221:
                ds=y[0]-x[0]
            else:
                ds=0
            I1=I1+(x[7]+y[7])/2/rho*ds  #x[7] ,y[7] represnt eta#
            I2=I2+1/rho**2*ds
            I3=I3+1/abs(rho)**3*ds
            I4=I4+(x[7]+y[7])/2/rho**3*ds
            I5=I5+(x[9]+y[9])/2/rho**3*ds  #x[9] ,y[9] represnt eta#
            cnt=cnt+1
    jx=1-I4/I2
#    print I2,jx
    jy=1
    jz=2+I4/I2
    U0=Cr/constants.pi/2*(E-0)**4*I2
    emit=Cq*gamma**2*I5/jx/I2
    sigmaE=np.sqrt(Cq*gamma**2*I3/jz/I2)
    ap=I1/twisspara[-1,0]
    return U0,emit,sigmaE,ap
def dfunction(zeta):
    f2=lambda x: np.log(x)*np.exp(-x)/x
    Term2=integrate.quad(f2,zeta,np.inf)[0]*0.5*zeta
    f3=lambda x: np.exp(-x)/x
    Term3=integrate.quad(f3,zeta,np.inf)[0]*0.5*(3*zeta-zeta*np.log(zeta)+2)
    Term1=-1.5*np.exp(-zeta)
    Dfunction = np.sqrt(zeta)*(Term1+Term2+Term3)
    return Dfunction
def get_optimumV(delta_acc,h,alphap,U0):
    Fq=delta_acc**2*constants.pi*h*abs(alphap)*E/U0  
    func=lambda q: np.sqrt(q**2-1)-np.arccos(1/q)-0.5*Fq  
    solutions=fsolve(func, 1)
    return solutions
def Cal_touschek(sigmaX,sigmaY,sigmaS,RFacc,Dvalue):
    N=Isingle*1E9
    re=2.8179403227E-15 #m
    tau_touschek=8*constants.pi*sigmaX*sigmaY*sigmaS*gamma**2*RFacc**3/Dvalue/N/re**2/constants.c
    return tau_touschek
    
def costFunction(Xin):
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
    #FINT=0.60614
    FINT=0.60614
    h=80
    MLS_start= [Drift(0)]
    L1=[Drift(1.5/30) for i in range(30)] #drift 1.5m#
    Q3M2K1=[ Fquad(Q3P2K1,0.2/10) for i in range(10)]
    L2=[Drift(0.15/3) for i in range(3)]
    Q2M2K1=[Dquad(Q2P2K1,0.2/10) for i in range(10)]
    L3=[Drift(0.4/10) for i in range(10)]+[Drift(0.025/5) for i in range(5)]
    BB1=[Bedge(1/rho,theta/2,0.025,FINT)]+[Sbend(rho,theta/48) for i in range(48)]+[Bedge(1/rho,theta/2,0.025,FINT)]
    
    L4=[Drift(1/20) for i in range(20)]+[Drift(0.075/3) for i in range(3)]
    Q1M2K1=[Fquad(Q1P2K1,0.2/10) for i in range(10)]
    L5=[Drift(0.35/7) for i in range(7)]
    Q1M1L2=[Fquad(Q1P1L2,0.2/10) for i in range(10)]
    L6=[Drift(1/20) for i in range(20)]+[Drift(0.075/3) for i in range(3)]
    BB2=[Bedge(1/rho,theta/2,0.025,FINT)]+[Sbend(rho,theta/48) for i in range(48)]+[Bedge(1/rho,theta/2,0.025,FINT)]

    L7=[Drift(0.4/10) for i in range(10)]+[Drift(0.025/5) for i in range(5)]
    Q2M1L2=[Dquad(Q2P1L2,0.2/10) for i in range(10)]
    L8=[Drift(0.15/3) for i in range(3)]
    Q3M1L2=[Fquad(Q3P1L2,0.2/10) for i in range(10)]
    L9=[Drift(3.25/65) for i in range(65)]
    section1=MLS_start+L1+Q3M2K1+L2+Q2M2K1+L3+BB1+L4+Q1M2K1+L5+Q1M1L2+L6+BB2+L7+Q2M1L2+L8+Q3M1L2+L9
    
    L10=[Drift(3.25/65) for i in range(65)]
    Q3M2L2=[Fquad(Q3P2L2,0.2/10) for i in range(10)]
    L11=[Drift(0.15/3) for i in range(3)]
    Q2M2L2=[Dquad(Q2P2L2,0.2/10) for i in range(10)]
    L12=[Drift(0.4/10) for i in range(10)]+[Drift(0.025/5) for i in range(5)]
    BB3=[Bedge(1/rho,theta/2,0.025,FINT)]+[Sbend(rho,theta/48) for i in range(48)]+[Bedge(1/rho,theta/2,0.025,FINT)]
    L13=[Drift(1/20) for i in range(20)]+[Drift(0.075/3) for i in range(3)]
    Q1M2L2=[Fquad(Q1P2L2,0.2/10) for i in range(10)]
    L14=[Drift(0.35/7) for i in range(7)]
    Q1M1K3=[Fquad(Q1P1K3,0.2/10) for i in range(10)]
    L15=[Drift(1/20) for i in range(20)]+[Drift(0.075/3) for i in range(3)]
    BB4=[Bedge(1/rho,theta/2,0.025,FINT)]+[Sbend(rho,theta/48) for i in range(48)]+[Bedge(1/rho,theta/2,0.025,FINT)]
    L16=[Drift(0.4/10) for i in range(10)]+[Drift(0.025/5) for i in range(5)]
    Q2M1K3=[Dquad(Q2P1K3,0.2/10) for i in range(10)]
    L17=[Drift(0.15/3) for i in range(3)]
    Q3M1K3=[ Fquad(Q3P1K3,0.2/10) for i in range(10)]
    L18=[Drift(1.5/30) for i in range(30)]
    section2=L10+Q3M2L2+L11+Q2M2L2+L12+BB3+L13+Q1M2L2+L14+Q1M1K3+L15+BB4+L16+Q2M1K3+L17+Q3M1K3+L18
    
    
    L19=[Drift(1.5/30) for i in range(30)] #drift 1.5m#
    Q3M2K3=[ Fquad(Q3P2K3,0.2/10) for i in range(10)]
    L20=[Drift(0.15/3) for i in range(3)]
    Q2M2K3=[Dquad(Q2P2K3,0.2/10) for i in range(10)]
    L21=[Drift(0.4/10) for i in range(10)]+[Drift(0.025/5) for i in range(5)]
    BB5=[Bedge(1/rho,theta/2,0.025,FINT)]+[Sbend(rho,theta/48) for i in range(48)]+[Bedge(1/rho,theta/2,0.025,FINT)]

    L22=[Drift(1/20) for i in range(20)]+[Drift(0.075/3) for i in range(3)]
    Q1M2K3=[Fquad(Q1P2K3,0.2/10) for i in range(10)]
    L23=[Drift(0.35/7) for i in range(7)]
    Q1M1L4=[Fquad(Q1P1L4,0.2/10) for i in range(10)]
    L24=[Drift(1/20) for i in range(20)]+[Drift(0.075/3) for i in range(3)]
    BB6=[Bedge(1/rho,theta/2,0.025,FINT)]+[Sbend(rho,theta/48) for i in range(48)]+[Bedge(1/rho,theta/2,0.025,FINT)]

    L25=[Drift(0.4/10) for i in range(10)]+[Drift(0.025/5) for i in range(5)]
    Q2M1L4=[Dquad(Q2P1L4,0.2/10) for i in range(10)]
    L26=[Drift(0.15/3) for i in range(3)]
    Q3M1L4=[Fquad(Q3P1L4,0.2/10) for i in range(10)]
    L27=[Drift(3.25/65) for i in range(65)]
    section3=L19+Q3M2K3+L20+Q2M2K3+L21+BB5+L22+Q1M2K3+L23+Q1M1L4+L24+BB6+L25+Q2M1L4+L26+Q3M1L4+L27
    
    
    L28=[Drift(3.25/65) for i in range(65)]
    Q3M2L4=[Fquad(Q3P2L4,0.2/10) for i in range(10)]
    L29=[Drift(0.15/3) for i in range(3)]
    Q2M2L4=[Dquad(Q2P2L4,0.2/10) for i in range(10)]
    L30=[Drift(0.4/10) for i in range(10)]+[Drift(0.025/5) for i in range(5)]
    BB7=[Bedge(1/rho,theta/2,0.025,FINT)]+[Sbend(rho,theta/48) for i in range(48)]+[Bedge(1/rho,theta/2,0.025,FINT)]
    L31=[Drift(1/20) for i in range(20)]+[Drift(0.075/3) for i in range(3)]
    Q1M2L4=[Fquad(Q1P2L4,0.2/10) for i in range(10)]
    L32=[Drift(0.35/7) for i in range(7)]
    Q1M1K1=[Fquad(Q1P1K1,0.2/10) for i in range(10)]
    L33=[Drift(1/20) for i in range(20)]+[Drift(0.075/3) for i in range(3)]
    BB8=[Bedge(1/rho,theta/2,0.025,FINT)]+[Sbend(rho,theta/48) for i in range(48)]+[Bedge(1/rho,theta/2,0.025,FINT)]
    L34=[Drift(0.4/10) for i in range(10)]+[Drift(0.025/5) for i in range(5)]
    Q2M1K1=[Dquad(Q2P1K1,0.2/10) for i in range(10)]
    L35=[Drift(0.15/3) for i in range(3)]
    Q3M1K1=[ Fquad(Q3P1K1,0.2/10) for i in range(10)]
    L36=[Drift(1.5/30) for i in range(30)]
    
    section4=L28+Q3M2L4+L29+Q2M2L4+L30+BB7+L31+Q1M2L4+L32+Q1M1K1+L33+BB8+L34+Q2M1K1+L35+Q3M1K1+L36
    
    seq=section1+section2+section3+section4
    
    Mcell =np.identity(6)
    s=0
    for x in reversed(seq):
        Mcell = np.dot(Mcell,x[0])
        s=s+x[1]
#    print Mcell,s
    
    
    cr_num = CheckStability(Mcell)
    
    if cr_num==1:
        Qh, Qv,mucellh,mucellv=Cal_Tunes(Mcell)
        Qh=Qh
        Qv=Qv
        betax0 = Mcell[0,1]/np.sin(mucellh)
        alphax0 = (Mcell[0,0]-Mcell[1,1])/2/np.sin(mucellh)
        betay0 = Mcell[2,3]/np.sin(mucellv)
        alphay0 = (Mcell[2,2]-Mcell[3,3])/2/np.sin(mucellv)
        etax0 = ((1 -Mcell[1,1])*Mcell[0,5] + Mcell[0,1]*Mcell[1,5])/2/(1-np.cos(mucellh))
        etaxp0 = ((1 -Mcell[0,0])*Mcell[1,5] + Mcell[1,0]*Mcell[0,5])/2/(1-np.cos(mucellh))
    
        s=0
        bagx0=np.array([betax0,alphax0,(1+alphax0**2)/betax0])
        bagy0=np.array([betay0,alphay0,(1+alphay0**2)/betay0])
        disp0=np.array([etax0,etaxp0,0,0,0,1])
        bagx=bagx0.copy()
        bagy=bagy0.copy()
        disp = disp0.copy()
        
        curl_H=bagx[2]*disp[0]**2+2*bagx[1]*disp[0]*disp[1]+bagx[0]*disp[1]**2
        twisspara=np.array([s,bagx[0],bagx[1],bagx[2],bagy[0],bagy[1],bagy[2],disp[0],disp[1],curl_H,0,0])
        twisspara=np.reshape(twisspara,(1,12))
        for x in seq:
            disp=np.dot(x[0],disp)
            bagx=np.dot(Twiss(x[0][0:2,0:2]),bagx)
            bagy=np.dot(Twiss(x[0][2:4,2:4]),bagy)
            curl_H=bagx[2]*disp[0]**2+2*bagx[1]*disp[0]*disp[1]+bagx[0]*disp[1]**2
            s=s+x[1]
            #s, betax, alphx, gammax,betay,alphy,gammay, eta,etap,curl_H.,#
            twisspara=np.vstack((twisspara,np.array([s,bagx[0],bagx[1],bagx[2],bagy[0],bagy[1],bagy[2],disp[0],disp[1],curl_H,x[1],x[2]])))
        np.savetxt('twiss.txt',twisspara)
        Qx,Qy=Cal_Q(Mcell,twisspara[:,0],twisspara[:,1],twisspara[:,4]) #Tunes calculation,mind the inputs#

        
        if Qx <3 or Qx >4 or Qy>3 or Qy<2 or np.max(twisspara[:,1])>25 or np.max(twisspara[:,4])>25 or np.max(abs(twisspara[:,7]))>2:
            tau=np.inf
            emit = np.inf
            QPD01=np.inf
            QPD00=np.inf
            EUV=np.inf
            EUVp=np.inf
            VUV=np.inf
            VUVp=np.inf
            IR=np.inf
            IRp=np.inf
            IRy=np.inf
            IRyp=np.inf
            THZ=np.inf
    #        Dcav=np.inf
            Dsep0=np.inf
            Dsep0p=np.inf
        else:
            U0,emit,sigmaE,alphap=Cal_RadiationIntegrals(twisspara)
            if emit>4E-8:
                tau=np.inf
                emit = np.inf
                QPD01=np.inf
                QPD00=np.inf
                EUV=np.inf
                EUVp=np.inf
                VUV=np.inf
                VUVp=np.inf
                IR=np.inf
                IRp=np.inf
                IRy=np.inf
                IRyp=np.inf
                THZ=np.inf
        #        Dcav=np.inf
                Dsep0=np.inf
                Dsep0p=np.inf
            else:
                Gacc=[]
                for i in range(len(twisspara)):
                    if i<=21 or i>=len(twisspara)-20:
                        if twisspara[i,7]==0:
                            Gacc=np.append(Gacc,np.inf)
                        else:
                            Gacc=np.append(Gacc,0.02/2/abs(twisspara[i,7]))
                    else:
                        if twisspara[i,7]==0:
                            Gacc=np.append(Gacc,np.inf)
                        else:
                            Gacc=np.append(Gacc,0.035/2/abs(twisspara[i,7]))
                RFacc=np.min(Gacc) 
                q_opt=get_optimumV(RFacc,h,alphap,U0)[0]
                Vopt=q_opt*U0##unit GV
                T0=twisspara[-1,0]/constants.c
                fs=np.sqrt(Vopt/E*h*abs(alphap)*np.sqrt(1-1/q_opt**2)/2/constants.pi)/T0
#                print fs,h
                beamsize=np.zeros((len(twisspara[:,0]),7))
                for i in range(len(twisspara)):
                    sigmaXbeta=np.sqrt(twisspara[i,1]*emit)                     #twisspara[i,1] betax#
                    sigmaXD=np.sqrt(twisspara[i,7]**2*sigmaE**2)       #twisspara[i,7] etax#
                    sigmaX=np.sqrt(sigmaXbeta**2+sigmaXD**2)
                    sigmaXp=np.sqrt(emit*twisspara[i,3]+twisspara[i,8]**2*sigmaE**2)  #twisspara[i,3]:gammax,  twisspara[i,8]: etap#
                    sigmaY=np.sqrt(couplingfactor*emit*twisspara[i,4])
                    sigmaYp=np.sqrt(couplingfactor*emit*twisspara[i,6])
                    sigmaS=abs(alphap)*constants.c/fs*sigmaE
                    zeta=(RFacc*twisspara[i,1]/gamma/sigmaX)**2
                    Dvalue=dfunction(zeta)
                    Tau_t=Cal_touschek(sigmaX,sigmaY,sigmaS,RFacc,Dvalue)
                    beamsize[i,0]=twisspara[i,0]
                    beamsize[i,1]=sigmaX
                    beamsize[i,2]=sigmaXp
                    beamsize[i,3]=sigmaY
                    beamsize[i,4]=sigmaYp
                    beamsize[i,5]=sigmaS
    ##                print sigmaX,sigmaY,sigmaS
                    beamsize[i,6]=Tau_t
                np.savetxt('beamsize.txt',beamsize)
                tau=0
                for i in range(len(twisspara)-1):
                    tau=tau+(beamsize[i,6]+beamsize[i+1,6])*0.5*(beamsize[i+1,0]-beamsize[i,0])      
                tau=-tau/3600/twisspara[-1,0]
                QPD01=beamsize[461,1]
                QPD00=beamsize[905,1]
                EUV=beamsize[584,1]
                EUVp=beamsize[584,2]
                VUV=beamsize[596,1]
                VUVp=beamsize[596,2]
                IR=beamsize[894,1]
                IRp=beamsize[894,2]
                IRy=beamsize[894,3]
                IRyp=beamsize[894,4]
                THZ=beamsize[770,1]
    #            Dcav=abs(twisspara[1033,3])
                Dsep0=abs(twisspara[0,7])
                Dsep0p=abs(twisspara[0,8])
##            print Qx,Qy
#            plt.plot(twisspara[:,0],twisspara[:,1],'b',twisspara[:,0],twisspara[:,4],'g',twisspara[:,0],twisspara[:,7]*10,'r')
##            plt.show()
    else:
        tau=np.inf
        emit = np.inf
        QPD01=np.inf
        QPD00=np.inf
        EUV=np.inf
        EUVp=np.inf
        VUV=np.inf
        VUVp=np.inf
        IR=np.inf
        THZ=np.inf
        IRp=np.inf
        IRy=np.inf
        IRyp=np.inf
#        Dcav=np.inf
        Dsep0=np.inf
        Dsep0p=np.inf
    F=np.array([emit,tau,QPD00,QPD01,IR,IRp,IRy,IRyp,EUV,EUVp,VUV,VUVp,Dsep0,Dsep0p])
    
    
    return F
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
	
Nobj = 14
lb = np.array([2,2,2,2,-5.2,-5.2,-5.2,-5.2,4.5,4.5,4.5,4.5])
ub =np.array([3.5,3.5,3.5,3.5,-3.5,-3.5,-3.5,-3.5,6,6,6,6])

maxit= 60     #Maximum Number of Iterations
             # Swarm size
nRep=100            # Repository Size
w=0.4               # Inertia Weight
wdamp=0.99          # Intertia Weight Damping Rate
c1=2.8                # Personal Learning Coefficient
c2=1.4             # Global Learning Coefficient

alpha=0.1          # Inflation Rate

mu=0.05            # Mutation Rate
vmax = (ub-lb)*0.1
S=200
D=len(ub)
X = np.random.rand(S, D)  # particle positions
V = np.zeros_like(X)  # particle velocities
F= np.zeros((S,Nobj))  # current particle function values
PXbest = np.zeros_like(X)  # best particle positions
PFbest = np.ones(S)*np.inf  # best particle function values

# Initialize the particle's position
X = lb + X*(ub - lb)
X[0,:]=np.array([2.7109,2.7109,2.7109,2.7109,-4.3580,-4.7055,-4.7055,-4.3580,5.3454,5.8778,5.8778,5.3454])
#X[0,:]=np.array([2.7109,2.7109,2.7109,2.7109,-4.3580,-4.7055,-4.7055,-4.3580,5.3454,5.8778,5.8778,5.3454])
processes=1
if processes > 1:
    F= np.array(mp_pool.map(costFunction, X))
else:
    for i in range(S):
        F[i]= costFunction(X[i, :])
        print F[i]
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
    for i in xrange(0, S):
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
        archiveX,archiveF=removeparticles(archiveX,archiveF,len(archiveX)-nRep)
    print 'the',it,'th iteration finished'


np.savetxt('archiveF.txt',archiveF)
np.savetxt('archiveX.txt',archiveX)
plt.plot(archiveF[:,0], archiveF[:,1],'ro')
plt.xlabel('$f_{1}$')
plt.ylabel('$f_{2}$')    
