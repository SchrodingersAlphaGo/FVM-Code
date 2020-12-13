# -*- coding:utf-8 -*-
import sys
import math
import numpy as np
import matplotlib.pyplot as plt



'''
plt.xlabel('Distance (m)')
plt.ylabel('A')
#plt.plot(range(P_NODES_NUM), xu_grid,'bs', label='area_i')
#plt.plot(range(P_NODES_NUM), xp_grid, 'go-', label='area_I')
plt.plot(xu_grid[1:], area_i[1:] ,'bs', label='area_i')
plt.plot(xp_grid, area_I,'go-', label='area_I')
#plt.title('velocity')
plt.legend()
plt.show()
'''


#== initialization ==



for j in range(MAX_ITERATION):
    #== SIMPLE ALGORITHM ==
    #== STEP 1 : solving momentum equations ==
    #== inlet ==
    #print(uu[i], uu0[i])

    #== internal nodes ==
    for i in range(2,U_NODES_NUM):
        aw = DENSITY * area_I[i-1] * (uu[i]+uu0[i-1])/2

        if i == U_NODES_NUM-1:
            ai = DENSITY * area_i[i] * uu[i]
        else:
            ai = DENSITY * area_I[i] * (uu[i]+uu[i+1])/2

        dd[i] = area_i[i]/ai
        uu_residual[i] = abs(ai*uu[i] - aw*uu[i-1] - area_i[i]*(pp[i-1]-pp[i]) )
        uu[i] = aw/ai * uu[i-1] + dd[i] * (pp[i-1]-pp[i]) 

    #print(uu[1:])

    ##== STEP 2: solving pressure correction equation ==
    for I in range(1,P_NODES_NUM-1):
        Aap[I,I-1] = -DENSITY * dd[I] * area_i[I]
        Aap[I,I+1] = -DENSITY * dd[I+1] * area_i[I+1]
        bb[I] = DENSITY * (area_i[I]*uu[I] - area_i[I+1]*uu[I+1])
        Aap[I,I] = -(Aap[I,I-1] + Aap[I,I+1])

    for i in [0, P_NODES_NUM-1]:
        Aap[i,i] = 1
        bb[i] = 0

    pp_cor = np.linalg.solve(Aap, bb)
    #print(pp_cor)

    ##== STEP 3: correct presure and velocity ==
    for i in range(1,U_NODES_NUM):
        uu[i] = uu[i] + ALPHA*dd[i]*(pp_cor[i-1] - pp_cor[i])
        
    pp = pp + ALPHA*pp_cor
    pp[0] = P_INLET - 0.5*DENSITY*uu[1]**2 * (area_i[1]/AREA_A)

    uu0 = uu.copy()
    if j%200==0:
        print('i = ',j, '   u_residul = ', uu_residual.sum())

    if uu_residual.sum() < 1e-5:
        print('i = ',j, '   u_residul = ', uu_residual.sum())
        break
    if uu_residual.sum() > 1e50:
        print('i = ',j, '   u_residul = ', uu_residual.sum())
        sys.exit()
        

#print(uu)
#print(pp)

print('mass flow rate = ', DENSITY*uu[-1]*area_i[-1])

##== exact solution ==
nodes_excat = 50
AREA_Exact = np.linspace(AREA_A, AREA_E, nodes_excat, endpoint=True)
xe_grid = np.linspace(0, LENGTH, nodes_excat, endpoint=True)
pP_ENLETxact = np.zeros(nodes_excat)
uu_exact = np.zeros(nodes_excat)
for i in range(nodes_excat):
    pP_ENLETxact[i] = P_INLET - 0.5 * 0.44721**2 /DENSITY / AREA_Exact[i]**2
    uu_exact[i] = math.sqrt((P_INLET - pP_ENLETxact[i])*2/DENSITY)


plt.xlabel('Distance (m)')
plt.ylabel('Velocity (m/s)')
plt.plot(xu_grid[1:], uu[1:] ,'bs--', label='Numerical')
plt.plot(xe_grid, uu_exact,'k', label='Exact')
plt.title('velocity')
plt.legend()


plt.figure(2)
plt.xlabel('Distance (m)')
plt.ylabel('Pressure (Pa)')
plt.plot(xp_grid, pp ,'bs--', label='Numerical')
plt.plot(xe_grid, pP_ENLETxact,'k', label='Exact')
plt.title('Pressure')
plt.legend()

plt.show()




#=== parameters ===
U_NODES_NUM = 30
P_NODES_NUM = U_NODES_NUM
ALPHA = 0.01
MAX_ITERATION = 15001

LENGTH = 2
DENSITY = 1
AREA_A = 0.5
AREA_E = 0.1
P_INLET = 10
P_ENLET = 0


xp_grid = np.linspace(0,LENGTH,P_NODES_NUM, endpoint=True)
xu_grid = np.zeros(U_NODES_NUM)
area_I = np.linspace(AREA_A, AREA_E, P_NODES_NUM, endpoint=True)
area_i = np.zeros(U_NODES_NUM)

au = np.zeros(U_NODES_NUM)
dd = np.zeros(U_NODES_NUM)
bb = np.zeros(P_NODES_NUM)
Aap = np.zeros((P_NODES_NUM, P_NODES_NUM))
Bbp = np.zeros(P_NODES_NUM)

uu = np.zeros(U_NODES_NUM)
pp = np.linspace(P_INLET, P_ENLET, P_NODES_NUM, endpoint=True)
uu0 = uu.copy()
uu_residual = np.zeros(U_NODES_NUM)

print('xu_grid = ',xu_grid)
#==========================
def grid():
    for i in range(1,U_NODES_NUM):
        area_i[i] = (area_I[i] + area_I[i-1])/2
        xu_grid[i] = (xp_grid[i] + xp_grid[i-1])/2

def initialization():
    guessed_m = 1
    uu[1:] = guessed_m / DENSITY / area_i[1:]
    grid()
    
def BC_inlet():
    i = 1
    fw = DENSITY * uu[i] * area_i[i]
    fe = DENSITY * area_I[i] * (uu[i] + uu[i+1])/2
    ai = fe + 0.5*fw*(area_i[i]/AREA_A)**2
    su = (P_INLET - pp[i])*area_i[i] + fw * area_i[i]/AREA_A * uu[i]
    uu_residual[i] = abs(ai*uu[i] - su)
    uu[i] = su / ai
    dd[i] = area_i[i]/ai

def BC_outlet():

def sol_pseudo_velocity():

def sov_pressure_equ():

def sol_moment_equ():

def sol_pressure_cor_equ():

def correction_uvp():

def sol_transport_equ():

def simple():

def simpleR():

def exact_solution():

def show_result():

#def simpleC():
#def piso():


def main():
    simpleR()

if __name__ == '__main__':
    main()
    