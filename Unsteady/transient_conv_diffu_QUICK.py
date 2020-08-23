import numpy as np
import matplotlib.pyplot as plt
import math
import sys

def source_term(xgrid):
    source = np.zeros(nx)
    for i in range(nx):
        if xgrid[i] < x1:
            source[i] = aaa*xgrid[i] + bbb
        elif xgrid[i] < x1+x2:
            source[i] = 100*xgrid[i] - 80
    return source

def exactSolution0(xgrid):
    phi = np.zeros(np.size(xgrid))
    for i in range(np.size(xgrid)):
        temp = 0
        x = xgrid[i]
        for j in range(1,sn):
            coef = L / (j * math.pi)
            coef2 = ppp * math.sin(x / coef) + math.cos(x / coef) / coef
            coef3 = ppp**2 + 1/coef**2
            temp = temp + an(j) * coef * coef2 / coef3
        phi[i] = cal_C1() + cal_C2() * math.exp(ppp*x) - a0()/ppp**2 * (ppp*x + 1) - temp
    return phi

def cal_C2():
    coef = a0() / math.exp(ppp*L)
    temp = coef / ppp**2
    for i in range(1,sn):
        coef2 = i * math.pi
        temp = temp + coef * math.cos(coef2) / (ppp**2 + (coef2/L)**2)
    return temp

def cal_C1():
    coef = a0()/ppp**2
    temp = -cal_C2() + coef
    for i in range(1,sn):
        coef2 = (i * math.pi / L)**2
        temp = temp + an(i) / (ppp**2 + coef2)
    return temp

def a0():
    return ((x1 + x2) * (aaa*x1 + bbb) + bbb*x1) / (2*L)

def an(n):
    '''
    coef = 2*L/(n**2 * math.pi**2)
    coef2 = (aaa*(x1+x2)+bbb)/x2 * math.cos(n * math.pi * x1 / L)
    coef3 = aaa + (aaa*x1 + bbb)/x2 * math.cos( n*math.pi*(x1+x2) / L )'''
    ann = 2*L / (n*math.pi)**2 * ( (aaa*(x1+x2) + bbb)/x2 * math.cos(n*math.pi*x1/L) \
        - (aaa + (aaa*x1 + bbb)/x2) * math.cos(n*math.pi*(x1+x2)/L) )
    #return coef * ( coef2 - coef3)
    return ann

def exactSolution(xgrid):
    nx = np.size(xgrid)
    phi = np.zeros(nx)
    a0 = ( (x1+x2)*(aaa*x1+bbb) + bbb*x1 ) / (2*L)

    temp = 0.0
    for i in range(1,sn):
        temp = temp + an(i)/math.exp(ppp*L) * math.cos(i*math.pi) / (ppp**2 + (i*math.pi/L)**2)
        #print(ppp,ppp*L,math.exp(ppp*L))
    c2 = a0/(ppp**2 * math.exp(ppp*L)) + temp
    print(c2)
    temp = 0.0
    for i in range(1, sn):
        temp = temp + an(i)/ (ppp**2 + (i*math.pi/L)**2)
    c1 = a0/ppp**2 - c2 + temp
    for j in range(nx):
        temp = 0
        x = xgrid[j]
        for i in range(1,sn):
            temp = temp + an(i) * L/i/math.pi * (ppp*math.sin(i*math.pi*x)/L + i*math.pi/L * math.cos(i*math.pi*x/L)) \
                / (ppp**2 + (i*math.pi/L)**2) 
        phi[j] = c1 + c2*math.exp(ppp*x) - a0/ppp**2 * (ppp*x+1) - temp
    return phi


    


    

nx = 45
dt  = 0.01

L = 1.5
u = 2
rho = 1
gamma = 0.03
phi_a = 0
aaa = -200
bbb = 100
x1 = 0.6
x2 = 0.2
ppp = rho * u / gamma
sn = 32

dx = L / nx
x_grid = np.linspace(dx/2, L-dx/2, nx, endpoint=True)
x_exact = np.linspace(0, L, 50, endpoint=True)
#print(x_grid)

dd = gamma / dx
ff = rho * u
ap0 = rho * dx / dt
dt_limit  = rho * dx**2 / (2*gamma)
print('dt max = ', dt_limit)

phi = np.ones(nx)
phi0 = phi.copy()
source = source_term(x_grid) * dx
AA = np.zeros((nx, nx))
bb = np.zeros(nx)
time = []
t_temp = 0

for j in range(100):
    if j%100 == 0: print('j = ',j)
    for i in range(nx):
        if i == 0:
            AA[i,i+1] = -(dd + dd / 3)
            sp = -(8/3*dd + ff)
            bb[i] = (8/3*dd + ff)*phi_a + 1/8*ff*(phi0[i] - 3*phi0[i+1]) + ap0*phi0[i] + source[i]
            AA[i,i] = -AA[i,i+1] - sp + ap0
        elif i == 2:
            AA[i,i-1] = -(dd + ff)
            AA[i,i+1] = -dd
            bb[i] = 1/8 * ff * (3*phi0[i] - phi0[i-1]) + 1/8 * ff * (phi0[i-1] + 2*phi0[i] \
                -3*phi0[i+1]) + ap0*phi0[i] + source[i]
            AA[i,i] = -AA[i,i-1] - AA[i,i+1] + ap0
        elif i == nx-1:
            AA[i,i-1] = -(dd + ff)
            bb[i] = 1/8 * ff * (3*phi0[i] - 2*phi0[i-1] -phi0[i-2]) + ap0*phi0[i] + source[i]
            AA[i,i] = -AA[i,i-1] + ap0
        else:
            AA[i,i-1] = -(dd + ff)
            AA[i,i+1] = -dd
            bb[i] = 1/8 * ff * (3*phi0[i] - 2*phi0[i-1] - phi0[i-2]) + 1/8 * ff * (phi0[i-1] \
                + 2*phi0[i] -3*phi0[i+1]) + ap0*phi0[i] + source[i]
            AA[i,i] = -AA[i,i-1] - AA[i,i+1] + ap0
    
  #  print(AA[0],bb)
    phi = np.linalg.solve(AA, bb)
    phi0 = phi.copy()

    t_temp = t_temp + dt
    time.append(t_temp)
  #  if t_temp >= target_time :
  #      print('time = %d s'% t_temp)
  #      break
    

#print('T(' + str(time[-1]) + ') = ', phi)

#sys.exit()

plt.xlabel('Distance (m)')
plt.ylabel('T (C)')
plt.plot(x_grid, phi ,'bs--', label='Numerical')
plt.plot(x_grid, exactSolution(x_grid),'r-', label='Exact')
title = 'time = %d'%time[-1] 
plt.title(title + ' s')
plt.legend()
plt.show()