from __future__ import division
from scipy import constants
import numpy as np 
from scipy.interpolate import RegularGridInterpolator


Cr=8.846E-5 ##unit m/GeV3
Cq=3.832E-13
echarge = -1.602176487e-19
m_e = 9.10938215e-31
c_light = 299792458


RWdata=np.loadtxt('fie12.map')


Bx=RWdata[:,5]
By=RWdata[:,4]
Bz=RWdata[:,3]


Px=np.unique(RWdata[:,2])
Py=np.unique(RWdata[:,1])
Pz=np.unique(RWdata[:,0])

Nx=int((np.max(Px)-np.min(Px))/(RWdata[1,2]-RWdata[0,2]))+1
Ny=int((np.max(Py)-np.min(Py))/(Py[1]-Py[0]))+1
Nz=int((np.max(Pz)-np.min(Pz))/(Pz[1]-Pz[0]))+1

Px=np.linspace(np.min(Px),np.max(Px),Nx)

U = np.zeros((Nz,Ny,Nx))
V = np.zeros((Nz,Ny,Nx))
W = np.zeros((Nz,Ny,Nx))
tt=0
for i in range(Nz):
    for j in range(Ny):
        for k in range(Nx):
            U[i,j,k] = Bx[tt]
            V[i,j,k] = By[tt]
            W[i,j,k] = Bz[tt]
            tt=tt+1

f1=RegularGridInterpolator((Pz, Py, Px), U)
f2=RegularGridInterpolator((Pz, Py, Px), V)
f3=RegularGridInterpolator((Pz, Py, Px), W)

#newmap=[]
#for x in Px:
#    for y in Py:
#        for z in Pz:
#            Bx = f1((z, y, x))
#            By = f2((z, y, x))
#            Bz = f3((z, y, x))
#            newmap.append([x,y,z,Bx,By,Bz])
#np.savetxt('newmap.txt',newmap)
            


        
class element(object):
    """
    element in the ring
    Attribute:
    E: energy of particles
    """
    gamma=1
    beta =1
    particle=np.array([0,0,0,0,0,0])
    def __init__(self,energy,f):
        self.energy = energy
        self.f=f
    def set_global(self):
        Emass = constants.physical_constants['electron mass energy equivalent in MeV'][0]*1E-3
        element.gamma= self.energy/Emass
        element.E = self.energy
        element.frf= self.f
 
class Drift(element):
    """
    element in the ring
    Attribute:
    name: element name,string
    Length: element length in meters
    Nslices: 10 default
    aperture: 1*2 array of elliptical aperture half-axes  in metres
    """
    def __init__(self,name,length,Nslices=1,aperture=np.array([10,10])):
        self.name=name
        self.length=length
        self.Nslices=Nslices
        self.aperture=aperture
    def transfermatrix(self):
        A=np.identity(6)
        A[0,1]= self.length/self.Nslices
        A[2,3]= self.length/self.Nslices
        A[4,5] =self.length/self.Nslices/self.gamma**2/self.beta**2  
        return [[A,self.length/self.Nslices,0] for i in range(self.Nslices)]
    def track(self):
        
        [x0, px0, y0, py0, ct0, dp0] = element.particle
        
        ds = self.length
        		
        d1  = np.sqrt(1  - px0 * px0 - py0 * py0 + 2 * dp0 / self.beta + dp0 * dp0)
        
        x1  = x0  + ds * px0 /d1 
        y1  = y0  + ds * py0 /d1 
        ct1 = ct0 + ds * (1 - (1 + self.beta * dp0) /d1)/self.beta
        
        element.particle = [x1, px0, y1, py0, ct1, dp0]

        
class Quadrupole(element):
    """
    element in the ring
    Attribute:
    name: element name,string
    strength: 
    Length: element length in meters
    Nslices: 10 default
    aperture: 1*2 array of elliptical aperture half-axes in metres
    """
    def __init__(self,name,strength,length,Nslices=10,aperture=np.array([10,10])):
        self.name = name
        self.strength = strength
        self.length = length
        self.Nslices = Nslices
        self.aperture = aperture
    def transfermatrix(self):
        A=np.identity(6)       
        if self.strength>0:
            K=np.sqrt(self.strength)
            a = K * self.length/self.Nslices
            A[0,0]= np.cos(a)
            A[0,1] = np.sin(a)/K
            A[1,0] = -K*np.sin(a)
            A[1,1] = np.cos(a)
            A[2,2] = np.cosh(a)
            A[2,3] = np.sinh(a)/K
            A[3,2] = np.sinh(a)*K
            A[3,3] = np.cosh(a)
            self.symbol=101
        if self.strength<0:
            K=np.sqrt(-self.strength)
            a = K * self.length/self.Nslices
            A[0,0]= np.cosh(a)
            A[0,1] = np.sinh(a)/K
            A[1,0] = np.sinh(a)*K
            A[1,1] = np.cosh(a)
            A[2,2] = np.cos(a)
            A[2,3] = np.sin(a)/K
            A[3,2] = -K*np.sin(a)
            A[3,3] = np.cos(a)    
            self.symbol=102
        if self.strength==0:
            A[0,1]= self.length/self.Nslices
            A[2,3]= self.length/self.Nslices
            self.symbol=0
        A[4,5] = self.length/self.gamma**2/self.beta**2
        return [[A,self.length/self.Nslices,self.symbol] for i in range(self.Nslices)]
    def track(self):
        
        [x0, px0, y0, py0, ct0, dp0] = element.particle
        if self.strength > 0: 
            [x0, px0, y0, py0, ct0, dp0] = element.particle

            ds = self.length
            k1 = self.strength
            
            d1 = np.sqrt(1 + 2 * dp0 /self.beta + dp0 * dp0)
            w  = np.sqrt(k1 / d1)
            
            xs  = np.sin( w * ds)
            xc  = np.cos( w * ds)
            ys  = np.sinh(w * ds)
            yc  = np.cosh(w * ds)
            xs2 = np.sin(2 * w * ds)
            ys2 = np.sinh(2 * w * ds)
            
            
            x1  =  x0 * xc       + px0 * xs * w / k1
            px1 = -k1 * x0 * xs / w + px0 * xc 
            y1  =  y0 * yc       + py0 * ys * w / k1 
            py1 =  k1 * y0 * ys / w + py0 * yc
            
            d0  = 1 / self.beta + dp0 
            d2  =-d0 / d1 / d1 / d1/2 
            
            c0  = (1 / self.beta - d0 / d1) * ds
            c11 = k1 * k1 * d2 *(xs2 /w - 2*ds) /w /w /4 
            c12 =-k1 * d2 * xs *xs /w /w
            c22 = d2 *(xs2 /w + 2 * ds) / 4 
            c33 = k1 * k1 * d2 *(ys2 /w - 2*ds) /w /w/4 
            c34 = k1 * d2 * ys * ys /w /w 
            c44 = d2 *(ys2 / w + 2 * ds)/4 
            
            ct1 = ct0 + c0 \
                      + c11 * x0 * x0 \
                      + c12 * x0 * px0 \
                      + c22 * px0 * px0 \
                      + c33 * y0 * y0 \
                      + c34 * y0 * py0 \
                      + c44 * py0 *py0 
                 
            element.particle= [x1, px1, y1, py1, ct1, dp0]

        elif self.strength < 0:

            ds = self.length;
            k1 = self.strength
            
            d1 = np.sqrt(1 + 2 * dp0 / self.beta + dp0 * dp0)
            w  = np.sqrt(abs(k1)  / d1) 
            
            xs  = np.sinh( w * ds) 
            xc  = np.cosh( w * ds) 
            ys  = np.sin( w * ds) 
            yc  = np.cos( w * ds) 
            xs2 = np.sinh(2 * w * ds) 
            ys2 = np.sin(2 * w * ds) 
            
            
            x1  =  x0 * xc       + px0 * xs * w / abs(k1) 
            px1 = -k1 * x0 * xs /w + px0 * xc 
            y1  =  y0 * yc       + py0 *ys * w/ abs(k1) 
            py1 =  k1 * y0 * ys /w + py0 *yc 
            
            d0  = 1/self.beta + dp0
            d2  =-d0 /d1 /d1 / d1/ 2 
            
            c0  = (1/self.beta - d0 /d1)*ds 
            c11 = k1 * k1 * d2 *(xs2 / w - 2 * ds) /w / w /4 
            c12 =-k1 * d2 * xs * xs /w /w 
            c22 = d2 *(xs2 /w + 2 * ds)/ 4 
            c33 = k1 * k1 * d2 *(ys2 /w - 2 * ds) /w / w /4 
            c34 = k1 * d2 * ys * ys /w /w 
            c44 = d2 *(ys2 /w + 2 * ds)/4 
            
            ct1 = ct0 + c0 \
                      + c11 * x0 * x0   \
                      + c12 * x0 * px0  \
                      + c22 * px0 * px0 \
                      + c33 * y0 * y0   \
                      + c34 * y0 * py0  \
                      + c44 * py0 *py0 
                 
            element.particle= [x1, px1, y1, py1, ct1, dp0]

        elif self.strength ==0: 
            [x0, px0, y0, py0, ct0, dp0] = element.particle
            		        
            d1  = np.sqrt(1 - px0 * px0 - py0 * py0 + 2 * dp0 / self.beta + dp0 * dp0)
            	 
            x1  = x0  + ds * px0 /d1 
            y1  = y0  + ds * py0 /d1 
            ct1 = ct0 + ds * (1 - (1 + self.beta * dp0) /d1)/self.beta
                    
            element.particle = [x1, px0, y1, py0, ct1, dp0]
    
class Sextupole(element):
    """
    element in the ring
    Attribute:
    name: element name,string
    strength: 
    Length: element length in meters
    Nslices: 10 default
    aperture: 1*2 array of elliptical aperture half-axes in metres
    """
    def __init__(self,name,strength,length,Nslices=10,aperture=np.array([10,10])):
        self.name=name
        self.strength = strength
        self.length = length
        self.Nslices = Nslices
        self.aperture = aperture
    def transfermatrix(self):
        A=np.identity(6)
        A[0,1] = self.length/self.Nslices
        A[2,3] = self.length/self.Nslices
        A[4,5] = self.length/self.Nslices/self.gamma**2/self.beta**2  
        return [[A,self.length/self.Nslices,0] for i in range(self.Nslices)]
    def track(self):
#        [x0,xp,y0,yp,ct0,dp0] = element.particle
        [x0,px0,y0,py0,ct0,dp0] = element.particle
        
        beta0  = self.beta;     
        ds = self.length;
        
        k2 = self.strength # normalised gradient

        # First apply a drift through 0.5 ds
        d1  = np.sqrt(1 - px0*px0 - py0*py0 + 2*dp0/beta0 + dp0*dp0)
        


        x1  = x0  + ds * px0/d1/2
        y1  = y0  + ds*py0/d1/2
        ct1 = ct0 + ds*(1 - (1 + beta0*dp0)/d1)/beta0/2
        
        # Next, apply a sextupole 'kick'

        px1 = px0 - (x1*x1 - y1*y1)*k2*ds/2
        py1 = py0 + x1*y1*k2*ds;
        
        # Finally, apply a second drift through ds/2
        d1  = np.sqrt(1 - px1*px1 - py1*py1 + 2*dp0/beta0 + dp0*dp0)

        x2  = x1  + ds*px1/d1/2
        y2  = y1  + ds*py1/d1/2
        ct2 = ct1 + ds*(1 - (1 + beta0*dp0)/d1)/beta0/2
        

        element.particle = [x2,px1,y2,py1,ct2,dp0]
        
        
class Bedge(element):
    """
    element in the ring
    Attribute:
    name: element name,string
    h: curvature
    e1: entrance angle of the dipole, usually half of the bending angle
    fint: fringel field factor
    Length: no length input, 0
    aperture: 1*2 array of elliptical aperture half-axes in metres
    """
    def __init__(self,name,h,e1,hgap=0.025,fint=0.5,aperture=np.array([10,10])):
        self.name = name
        self.h = h
        self.rho = 1/h
        self.e1 = e1
        self.hgap = hgap
        self.fint = fint
        self.aperture=aperture
        self.symbol = 201
        self.length = 0
    def transfermatrix(self):
        A=np.identity(6)       
        _psi =2*self.h*self.hgap*self.fint*(1+np.sin(self.e1)*np.sin(self.e1))/np.cos(self.e1)
        A[1,0]= np.tan(self.e1)*self.h
        A[3,2] =-np.tan(self.e1-_psi)*self.h
        return [[A,0,self.symbol]]
    def track(self):
        
        [x0, px0, y0, py0, ct0, dp0] =  element.particle

        k0     = self.h
            
        # First, apply a map for the entrance fringe field
        sine1 = np.sin(self.e1)
        _psi   = 2*self.fint*self.hgap*k0*(1+sine1*sine1)/np.cos(self.e1)
        r10   = k0*np.tan(self.e1)
        r32   =-k0*np.tan(self.e1 - _psi)
        
        px1   = px0 + r10*x0
        py1   = py0 + r32*y0
        
        element.particle = [x0,px1,y0,py1,ct0,dp0]


class Sbend(element):
    """
    bending magnet element in the ring, without gradient
    Attribute:
    name: element name,string
    h: curvature
    e1: entrance angle of the dipole, usually half of the bending angle
    fint: fringel field factor
    Length: no length input, 0
    aperture: 1*2 array of elliptical aperture half-axes in metres
    """
    def __init__(self,name,rho,theta,Nslices=10,gradient=0,aperture=np.array([10,10])):
        self.name = name
        self.rho = rho
        self.theta = theta
        self.symbol = 202
        self.length = rho*theta
        self.Nslices = Nslices
        self.gradient = gradient
        self.aperture=aperture
    def transfermatrix(self):
        A=np.identity(6)
        A[0,0]= np.cos(self.theta/self.Nslices)
        A[0,1] = self.rho*np.sin(self.theta/self.Nslices)
        A[0,5] = self.rho*(1-np.cos(self.theta/self.Nslices))/self.beta
        A[1,0] = -np.sin(self.theta/self.Nslices)/self.rho
        A[1,1] = np.cos(self.theta/self.Nslices)
        A[1,5] = np.sin(self.theta/self.Nslices)/self.beta
        A[2,2] = 1
        A[2,3] = self.length/self.Nslices
        A[3,3] = 1
        A[4,0] = -np.sin(self.theta/self.Nslices)/self.beta
        A[4,1] = -self.rho*(1-np.cos(self.theta/self.Nslices))/self.beta
        A[4,4] = 1
        A[4,5] = self.rho*(np.sin(self.theta/self.Nslices)-self.beta**2*self.theta/self.Nslices)
        A[5,5] = 1
        return [[A,self.length/self.Nslices,self.symbol] for i in range(self.Nslices)]
    def track(self):

        [x0,px0,y0,py0,ct0,dp0] = element.particle
        beta0  = self.beta   
        ds = self.length
        k0 = 1/self.rho
        d1     = np.sqrt(1 + 2*dp0/beta0 + dp0*dp0)
        
        h   = 1/self.rho
        k1  = self.gradient
        a1  = h - k0/d1
        
        wx  = np.sqrt((h*k0 + k1)/d1)
        xc  = np.cos(wx*ds)
        xs  = np.sin(wx*ds)/wx
        xs2 = np.sin(2*wx*ds)/wx
        
        wy  = np.sqrt(k1/d1)
        yc  = np.cosh(wy*ds)
        ys  = ds
        ys2 = 2*ds

        if wy!=0:
            ys  = np.sinh(wy*ds)/wy
            ys2 = np.sinh(2*wy*ds)/wy

        
        x1  =             x0 * xc + px0 * xs/d1 + a1*(1-xc)/wx/wx
        px1 = -d1 * wx * wx * x0 * xs + px0 * xc     + a1 * xs * d1
        
        y1  =             y0*yc + py0*ys/d1
        py1 = d1 * wy * wy * y0 * ys + py0 * yc

        d0  = 1/beta0 + dp0
        
        c0  = (1/beta0 - d0/d1)*ds - d0*a1*(h*(ds-xs) + a1*(2*ds-xs2)/8)/wx/wx/d1
                   
        c1  =-d0*(h*xs - a1*(2*ds-xs2)/4)/d1
               
        c2  =-d0*(h*(1-xc)/wx/wx + a1*xs*xs/2)/d1/d1
               
        c11 =-d0*wx*wx*(2*ds-xs2)/d1/8
        c12 = d0*wx*wx*xs*xs/d1/d1/2
        c22 =-d0*(2*ds+xs2)/d1/d1/d1/8
        
        c33 =-d0*wy*wy*(2*ds-ys2)/d1/8
        c34 =-d0*wy*wy*ys*ys/d1/d1/2
        c44 =-d0*(2*ds+ys2)/d1/d1/d1/8

        ct1 = ct0 + c0 + \
                    c1 * x0      + c2 * px0 + \
                    c11 * x0 * x0 + c12 * x0 * px0 + c22 * px0 * px0 + \
                    c33 * y0 * y0 + c34 * y0 * py0 + c44 * py0 *py0
                    
        element.particle=[x1,px1,y1,py1,ct1,dp0]
        
                    
                    
class Fieldmap(element):
    """
    element in the ring
    Attribute:
    name: element name,string
    strength: 
    Length: element length in meters
    Nslices: 10 default
    aperture: 1*2 array of elliptical aperture half-axes in metres
    """
    def __init__(self,name,length,Nslices=10,Nsteps=10000,aperture=np.array([10,10])):
        self.name = name

        self.length = np.max(RWdata[:,0])-np.min(RWdata[:,0])
        self.Nsteps = Nsteps
        self.aperture = aperture
        self.Nslices = Nslices 
        self.symbol = 300

    def transfermatrix(self):
        A = np.loadtxt('RWmatrix1.txt')
        
        return [[A,self.length,self.symbol]]              
        
    def track(self):
        
        [x0, px0, y0, py0, ct0, dp0] = element.particle
        
        par = []
        
        z0      = 0
        
        
        P0      = self.beta * self.gamma
        beta0   = self.beta 
        gamma0  = self.gamma 
        
        gamma   = gamma0 * (1 + beta0 * dp0 ) 
        
        b0      = np.sqrt(1 - gamma**(-2)) 
        bx0     = P0 * px0 / gamma 
        by0     = P0 * py0 / gamma 
        bz0     = np.sqrt(b0 * b0 - bx0 * bx0 - by0 * by0)
        
        k       = -(echarge / m_e / c_light ) / gamma
        
        for n in range(self.Nsteps):
       
            cdt   = (self.length - z0)  / (self.Nsteps  - n)  / bz0 
            if x0 > np.max(self.fdata[:,2]) or x0 < np.min(self.fdata[:,2]) or y0 > np.max(self.fdata[:,1]) or y0 < np.min(self.fdata[:,1]) :
                x0 =np.inf
                y0 =np.inf
                print 'break'
                break
                
            else:
                Bx0   = k* f1((z0,y0,x0))
                By0   = k* f2((z0,y0,x0))
                Bz0   = k* f3((z0,y0,x0))
            
#                Bx0   = k* self.f1(z0,y0,x0)
#                By0   = k* self.f2(z0,y0,x0)
#                Bz0   = k* self.f3(z0,y0,x0)

                bmag  = np.sqrt(Bx0 * Bx0 + By0 * By0 + Bz0 * Bz0) 
    
                bdotv = Bx0 * bx0 + By0 * by0 + Bz0 * bz0

                s1    = np.sin(bmag * cdt) / bmag 
                s2    = bdotv * ( bmag * cdt - np.sin(bmag * cdt)) *bmag **(-3) 
                c1    = np.cos(bmag * cdt) 
                c2    = (1 - c1) * bmag**(-2) 
                
    
                x1  = x0 + bx0 *s1 + (by0 *Bz0 - bz0 *By0) *c2 + Bx0 *s2 
                y1  = y0 + by0 *s1 + (bz0 *Bx0 - bx0 *Bz0) *c2 + By0 *s2 
                z1  = z0 + bz0 *s1 + (bx0 *By0 - by0 *Bx0) *c2 + Bz0 *s2 
                ct0 = ct0 + (bz0/beta0 - 1)  * cdt 
    
                bx1 = bx0 *c1 + (by0 *Bz0 - bz0 *By0) *s1 + Bx0 *bdotv *c2 
                by1 = by0 *c1 + (bz0 *Bx0 - bx0 *Bz0) *s1 + By0 *bdotv *c2 
                bz1 = bz0 *c1 + (bx0 *By0 - by0 *Bx0) *s1 + Bz0 *bdotv *c2 
    
                x0  = x1
                y0  = y1
                z0  = z1
    
                bx0 = bx1
                by0 = by1
                bz0 = bz1
                
                par.append([x0,y0,z0])
        
            
        
        px0 = bx0 * gamma / P0
        py0 = by0 * gamma / P0

        element.particle=[x0,px0,y0,py0,ct0,dp0]
        return par
        
class rfcsavity(element):
    """
    element in the ring
    Attribute:
    name: element name,string
    strength: 
    Length: element length in meters
    Nslices: 10 default
    aperture: 1*2 array of elliptical aperture half-axes in metres
    """
    def __init__(self,name,fdata,f1,f2,f3,length,Nslices=10,aperture=np.array([10,10])):
        self.name = name
        self.fdata= fdata

        self.length = np.max(fdata[:,0])-np.min(fdata[:,0])
        self.Nsteps = Nsteps
        self.aperture = aperture
        self.Nslices = Nslices 

        self.f1 = f1
        self.f2 = f2
        self.f3 = f3
    def transfermatrix(self):
        A=np.identity(6)
        A[0,1] = self.length/self.Nslices
        A[2,3] = self.length/self.Nslices
        A[4,5] = self.length/self.Nslices/self.gamma**2/self.beta**2  
        return [[A,self.length/self.Nslices,0] for i in range(self.Nslices)]                    
        
    def track(self):
        
        [x0, px0, y0, py0, ct0, dp0] = element.particle
#        
#        beta0 = self.beta
#        
#        ds  = self.length;
#        f   = rfcavity.harmonic*mofreq
#
#        # First apply a drift through ds/2
#        d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp0/beta0 + dp0.*dp0);
#
#        x1  = x0  + ds*px0./d1/2;
#        y1  = y0  + ds*py0./d1/2;
#        ct1 = ct0 + ds*(1 - (1 + beta0*dp0)./d1)/beta0/2;
#        
#        # Next, apply an rf 'kick'
#        p = floor(beam.globaltime*mofreq);
#        beam.globaltime = beam.globaltime - p/mofreq;
#        t = beam.globaltime - ct1/(beta0*PhysicalConstants.SpeedOfLight);
#        ft = f*t - floor(f*t);
#
#        vnorm = rfcavity.voltage/beam.rigidity/PhysicalConstants.SpeedOfLight;
#        dp1   = dp0 + vnorm*sin(2*pi*ft + rfcavity.phase);
#        
#        # Finally, apply a second drift through ds/2
#        d1  = sqrt(1 - px0.*px0 - py0.*py0 + 2*dp1/beta0 + dp1.*dp1);
#
#        x2  = x1  + ds*px0./d1/2;
#        y2  = y1  + ds*py0./d1/2;
#        ct2 = ct1 + ds*(1 - (1 + beta0*dp1)./d1)/beta0/2;
#        
#        element.particle=[x0,px0,y0,py0,ct0,dp0]
        

def CheckStability(Mcell):
    if (np.abs(Mcell[0,0]+Mcell[1,1])>=2) or (np.abs(Mcell[2,2]+Mcell[3,3])>=2):
        print 'Mcell is unstable.\n'
        return 0
    else:
        return 1

def Cal_QFractions(Mcell):
    
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
    
def Cal_RadiationIntegrals(twisspara):
    """
    U0: one turn radiation loss
    emit: emittance
    sigmaE: energy spread
    ap: momentum compaction factor
    Intergrals: I1,....I5
    Jxyz: Damping portion numbers [Jx,Jy,Jz]
    """
    I1=0
    I2=0
    I3=0
    I4=0
    I5=0
    cnt=0
    for i in range(len(twisspara)):
        x=twisspara[i,:]
        
        if x[2]==201 or x[2]==202:
            y=twisspara[i+1,:]
            if y[2]==201 or y[2]==202:
                ds=y[0]-x[0]
            else:
                ds=0
            I1=I1+(x[10]+y[10])/2*twisspara[i,3]*ds  #x[9] ,y[9] represnt eta#
            I2=I2+twisspara[i,3]**2*ds
            I3=I3+1/abs(twisspara[i,3])**3*ds
            I4=I4+(x[10]+y[10])/2*twisspara[i,3]**3*ds
            I5=I5+(x[12]+y[12])/2*twisspara[i,3]**3*ds  #x[12] ,y[12] represnt curl_H#
            cnt=cnt+1
    Jx=1-I4/I2
    Jy=1
    Js=2+I4/I2
    U0=Cr/constants.pi/2*(element.E)**4*I2
    emit=Cq*element.gamma**2*I5/Jx/I2
    sigmaE=np.sqrt(Cq*element.gamma**2*I3/Js/I2)
    ap=I1/twisspara[-1,0]
    print sigmaE
    return U0,emit,sigmaE,ap,[I1,I2,I3,I4,I5],[Jx,Jy,Js]
    
def Ringtwiss(seq):
    Mcell = np.identity(6)
    SS=0
    for x in reversed(seq):
        SS+=x.length
        for y in reversed(x.transfermatrix()):
            Mcell = np.dot(Mcell,y[0])
    print SS
    cr_num = CheckStability(Mcell)
    
    
    if cr_num==1:
        Qh, Qv,mucellh,mucellv=Cal_QFractions(Mcell)
        print Qh,Qv
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
    #    
        curl_H=bagx[2]*disp[0]**2+2*bagx[1]*disp[0]*disp[1]+bagx[0]*disp[1]**2
        twisspara=np.array([s,0,0,0,bagx[0],bagx[1],bagx[2],bagy[0],bagy[1],bagy[2],disp[0],disp[1],curl_H])
        twisspara=np.reshape(twisspara,(1,13))
        for y in seq:
            for x in y.transfermatrix():
                disp=np.dot(x[0],disp)
                bagx=np.dot(Twiss(x[0][0:2,0:2]),bagx)
                bagy=np.dot(Twiss(x[0][2:4,2:4]),bagy)
                curl_H=bagx[2]*disp[0]**2+2*bagx[1]*disp[0]*disp[1]+bagx[0]*disp[1]**2
                s=s+x[1]
                if x[2]==201 or x[2]==202 :
                    curvature=1/y.rho
                else:
                    curvature=0
                #s,length,element symbol number,curvature, betax, alphx, gammax,betay,alphy,gammay, eta,etap,curl_H,#
                twisspara=np.vstack((twisspara,np.array([s,x[1],x[2],curvature,bagx[0],bagx[1],bagx[2],bagy[0],bagy[1],bagy[2],disp[0],disp[1],curl_H])))
        np.savetxt('twiss.txt',twisspara)
        Qx,Qy=Cal_Q(Mcell,twisspara[:,0],twisspara[:,4],twisspara[:,7]) #Tunes calculation,mind the inputs#
        return twisspara,[Qx,Qy]
    else:
        print 'not stable!!!!!!!!'
def Beamsize(twisspara,U0,emit,sigmaE,alphap,V,couplingfactor):
    """
    tiwsspara: generated from Ringtwiss
    U0: one turn loss, from Cal_RadiationIntegrals
    emit: emittance, from Cal_RadiationIntegrals
    sigmaE: energy spread, from Cal_RadiationIntegrals
    alphap: momentum compaction factor, from Cal_RadiationIntegrals
    V: peak voltage in the cavity
    couplingfactor:
    """
    q=V/U0  #over voltage factor#
    T0=twisspara[-1,0]/constants.c  #single pass period#
    h=round(element.frf*T0)           #hamonic number#
    print h
    beamsize=np.zeros((len(twisspara[:,0]),7))
    fs=np.sqrt(V/element.E*h*abs(alphap)*np.sqrt(1-1/q**2)/2/constants.pi)/T0
    sigmaS=abs(alphap)*constants.c*element.beta/fs*sigmaE/2/constants.pi
    for i,x in enumerate(twisspara):
        sigmaXbeta=np.sqrt(x[4]*emit)                     #twisspara[i,4] betax#
        sigmaXD=np.sqrt(x[10]**2*sigmaE**2)       #twisspara[i,10] etax#
        sigmaX=np.sqrt(sigmaXbeta**2+sigmaXD**2)
        sigmaXp=np.sqrt(emit*x[6]+x[11]**2*sigmaE**2)  #twisspara[i,6]:gammax,  twisspara[i,11]: etap#
        sigmaY=np.sqrt(couplingfactor*emit*twisspara[i,7])
        sigmaYp=np.sqrt(couplingfactor*emit*twisspara[i,9])
        beamsize[i,0]=twisspara[i,0]
        beamsize[i,1]=sigmaX
        beamsize[i,2]=sigmaXp
        beamsize[i,3]=sigmaY
        beamsize[i,4]=sigmaYp
        beamsize[i,5]=sigmaS
    np.savetxt('beamsize.txt',beamsize)
    return beamsize
def track(seq,beam0,Nturn):

    element.particle=beam0
    Pdist=[]
    for i in range(Nturn):
        for x in seq:
            par = x.track()
            if element.particle[0]**2 +element.particle[3]**2 > x.aperture[0]**2+x.aperture[1]**2:
                return 0
            
        Pdist.append(element.particle)
        
        
#            print x.track()
#        np.append(Pdist,element.particle)
#        Pdist.append(element.particle)
    return [Pdist,np.array(par)]    
        
        
          