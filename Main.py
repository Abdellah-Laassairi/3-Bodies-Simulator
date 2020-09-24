import time
import math
import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from math import *
from matplotlib import animation, rc

#Animation function
def init():
    line1.set_data([], [])
    line2.set_data([], [])
    ttl.set_text('')
    
    return (line1,line2,ttl)
#forces
#Force of FSatellite1/sun
def force_es(r):
    F = np.zeros(2)
    Fmag = GG*Me*Ms/(np.linalg.norm(r)+1e-20)**2
    theta = math.atan(np.abs(r[1])/(np.abs(r[0])+1e-20))
    F[0] = Fmag * np.cos(theta)
    F[1] = Fmag * np.sin(theta)
    if r[0] > 0:
        F[0] = -F[0]
    if r[1] > 0:
        F[1] = -F[1]
        
    return F
#Force of FSatellite2/sun
def force_js(r):
    F = np.zeros(2)
    Fmag = GG*Mj*Ms/(np.linalg.norm(r)+1e-20)**2
    theta = math.atan(np.abs(r[1])/(np.abs(r[0])+1e-20))
    F[0] = Fmag * np.cos(theta)
    F[1] = Fmag * np.sin(theta)
    if r[0] > 0:
        F[0] = -F[0]
    if r[1] > 0:
        F[1] = -F[1]
        
    return F


#PFD
def force(r,satellite,ro,vo):
    if satellite == 'SatteliteA':
        return force_es(r) 
    if satellite == 'SatteliteB':
        return force_js(r) 

#Vitesse    
def dr_dt(t,r,v,satellite,ro,vo):
    return v
 
#Acceleration  
def dv_dt(t,r,v,satellite,ro,vo):
    F = force(r,satellite,ro,vo)
    if satellite == 'SatteliteA':
        y = F/Me
    if satellite == 'SatteliteB':
        y = F/Mj
    return y

#Differential equations solver : Resolution des equations differentiels

def EulerCromerSolver(t,r,v,h,satellite,ro,vo):
    z = np.zeros([2,2])
    r = r + h*dr_dt(t,r,v,satellite,ro,vo)
    v = v + h*dv_dt(t,r,v,satellite,ro,vo)
    z = [r, v]
    return z

def RK4Solver(t,r,v,h,satellite,ro,vo):
    k11 = dr_dt(t,r,v,satellite,ro,vo) 
    k21 = dv_dt(t,r,v,satellite,ro,vo)
    
    k12 = dr_dt(t + 0.5*h,r + 0.5*h*k11,v + 0.5*h*k21,satellite,ro,vo)
    k22 = dv_dt(t + 0.5*h,r + 0.5*h*k11,v + 0.5*h*k21,satellite,ro,vo)
    
    k13 = dr_dt(t + 0.5*h,r + 0.5*h*k12,v + 0.5*h*k22,satellite,ro,vo)
    k23 = dv_dt(t + 0.5*h,r + 0.5*h*k12,v + 0.5*h*k22,satellite,ro,vo)
    
    k14 = dr_dt(t + h,r + h*k13,v + h*k23,satellite,ro,vo)
    k24 = dv_dt(t + h,r + h*k13,v + h*k23,satellite,ro,vo)
    
    y0 = r + h * (k11 + 2.*k12 + 2.*k13 + k14) / 6.
    y1 = v + h * (k21 + 2.*k22 + 2.*k23 + k24) / 6.
    
    z = np.zeros([2,2])
    z = [y0, y1]
    return z
#Cinetic, Kinetic Energy and Momentum
def Energyc(v):
    vn = np.linalg.norm(v)
    return 0.5*Me*vn**2

def Energyp(r):
    fmag = np.linalg.norm(force_es(r))
    rmag = np.linalg.norm(r)
    return -fmag*rmag

def Momentc(r,v):
    rn = np.linalg.norm(r)
    vn = np.linalg.norm(v)
    r = r/rn
    v = v/vn
    rdotv = r[0]*v[0]+r[1]*v[1]
    theta = 0
    return theta

def mplot(fign,x,y,xl,yl,clr,lbl):
    py.figure(fign)
    py.xlabel(xl)    
    py.ylabel(yl)
    return py.plot(x,y,clr, linewidth =1.0,label = lbl)

#Units
th    = 60
rad = pi/180.*th
Me = 7.348 * 1e0              # Mass of satteliteA in kg
Ms = 2e30                     # Mass of Sun in kg                       
Mj = 7.348 * 1e0              # Mass of SatteliteB
G = 6.673e-11                 # Gravitational Constant
RR = 1.496e11                 # Normalizing distance in km (= 1 AU)
MM = 6e24                     # Normalizing mass
TT = 365*24*60*60.0           # Normalizing time (1 year)
FF = (G*MM**2)/RR**2          # Unit force
EE = FF*RR                    # Unit energy
GG = (MM*G*TT**2)/(RR**3)
Me = 1000                  # Normalized mass of Earth
Ms = Ms/MM                    # Normalized mass of Sun  
Mj = 1000              # Normalized mass of Jupiter/Super Jupiter

ti = 0                # initial time = 0
tf = 120              # final time = 120 years




N = 100*tf                   # 100 points per year
t = np.linspace(ti,tf,N)     # time array from ti to tf with N points 

h = t[2]-t[1]    # time step (uniform means that time is subdivided evenly)


# Initialization

EC = np.zeros(N)            # Energue cinetique
EP = np.zeros(N)            # Potential energy
MC = np.zeros(N)            # Moment cinetique
AreaVal = np.zeros(N)


r = np.zeros([N,2])         # position vector of SatteliteA
v = np.zeros([N,2])         # velocity vector of SatteliteA
rj = np.zeros([N,2])        # position vector of SatteliteB
vj = np.zeros([N,2])        # velocity vector of SatteliteB
ri = [-5,-4]      # initial position of SatteliteA
rji = [-4,-3]           # initial position of SAtteleteB


vi = [6, 6]                  # Initial velocity vector for Earth.Taken to be along y direction as ri is on x axis.
vji = [6, 6]                # Initial velocity vector for Jupiter


# initial values
t[0] = ti
r[0,:] = ri
v[0,:] = vi
rj[0,:] = rji
vj[0,:] = vji
EC[0] = Energyc(v[0,:])
EP[0] = Energyp(r[0,:])
MC[0] = Momentc(r[0,:],v[0,:])


#this calculates the speed and position using the damned runge kutta the 4th methode, you can change it to euler or the eulerCromer to see the difference
for i in range(0,N-1):
    [r[i+1,:],v[i+1,:]]=EulerCromerSolver(t[i],r[i,:],v[i,:],h,'SatteliteA',rj[i,:],vj[i,:])
    [rj[i+1,:],vj[i+1,:]]=RK4Solver(t[i],rj[i,:],vj[i,:],h,'SatteliteB',r[i,:],v[i,:])

#Animation function. Reads out the positon coordinates sequentially
def animate(i):
    Atrail = 400; #the lenght of the tail it leaves
    Btrail = 400; #same here
    tm_yr = 'Elapsed time = ' + str(round(t[i],1)) + ' years'
    ttl.set_text(tm_yr)
    line1.set_data(r[i:max(1,i-Atrail):-1,0], r[i:max(1,i-Atrail):-1,1])
    line2.set_data(rj[i:max(1,i-Btrail):-1,0], rj[i:max(1,i-Btrail):-1,1])
    

    return (line1,line2)

lbl = 'orbit'

mplot(3,t,v,r'Time, $t$ (years)',r'Vitesse satteliteA','blue',lbl)
mplot(3,t,vj,r'Time, $t$ (years)',r'Vitesse SatteliteB','red',lbl)
py.ylim([4, 8])

# Function for setting up the animation

fig, ax = py.subplots()
ax.axis('square')
ax.set_xlim(( -9, 9))
ax.set_ylim((-9, 9))
ax.get_xaxis().set_ticks([])    # enable this to hide x axis ticks
ax.get_yaxis().set_ticks([])    # enable this to hide y axis ticks

ax.plot(0,0,'o',markersize = 9, markerfacecolor = "#FDB813",markeredgecolor ="#FD7813" )
line1, = ax.plot([], [], 'o-',color = '#d2eeff',markevery=10000, markerfacecolor = '#0077BE',lw=2)   # line for sattelite1
line2, = ax.plot([], [], 'o-',color = '#e3dccb',markersize = 8, markerfacecolor = '#f66338',lw=2,markevery=10000)   # line for satellite2


ax.plot([-6,-5],[6.5,6.5],'r-')
ax.text(-4.5,6.3,r'1 AU = $1.496 \times 10^8$ km')

ax.plot(-6,-6.2,'o', color = '#d2eeff', markerfacecolor = '#0077BE')
ax.text(-5.5,-6.4,'sattelite')


ax.plot(-1.5,-6.2,'o', color = '#e3dccb',markersize = 8, markerfacecolor = '#f66338')
ax.text(-1,-6.4,'sattelite2')

ax.plot(5,-6.2,'o', markersize = 9, markerfacecolor = "#FDB813",markeredgecolor ="#FD7813")
ax.text(5.5,-6.4,'Sun')
ttl = ax.text(0.24, 1.05, '', transform = ax.transAxes, va='center')

#Plt.title('Elapsed time, T=%i years' %u)


#This last line just calls out the animation (like return f )

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=4000, interval=5, blit=True)
plt.show()
