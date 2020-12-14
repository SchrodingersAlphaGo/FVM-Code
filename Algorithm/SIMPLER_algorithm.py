# -*- coding:utf-8 -*-
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
    

#==========================
def grid():
    for i in range(1,U_NODES_NUM):
        area_i[i] = (area_I[i] + area_I[i-1])/2  # area of sections of 1,2,3,4 ( Velocity nodes )
        xu_grid[i] = (xp_grid[i] + xp_grid[i-1])/2


def initialization():
    guessedMassFlowRate = 1
    uu[1:] = guessedMassFlowRate / DENSITY / area_i[1:]
    global uu0
    uu0 = uu.copy()


def sol_pseudo_velocity():
    i = 1
    fw = DENSITY * uu[i] * area_i[i]
    fe = DENSITY * area_I[i] * (uu[i] + uu[i+1])/2
    ai = fe + 0.5 * fw * (area_i[i]/AREA_A)**2
    su = fw * area_i[i]/AREA_A * uu[i] # + (P_INLET - pp[i])*area_i[i]
    #psd_u[i] = area_i[i] * area_I[i] / AREA_A * uu[i]
    psd_u[i] = su / ai
    dd[i] = area_i[i] / ai

    for i in range(2, U_NODES_NUM):
        aw = DENSITY * area_I[i-1] * (uu[i] + uu0[i-1]) / 2

        if i == U_NODES_NUM-1:
            ai = DENSITY * area_i[i] * uu[i]
        else:
            ai = DENSITY * area_I[i] * (uu[i] + uu[i+1])/2

        psd_u[i] = aw/ai * uu[i-1]
        dd[i] = area_i[i] / ai


def sol_pressure_equ():
    for I in range(1,P_NODES_NUM-1):
        Aap[I,I-1] = -DENSITY * dd[I] * area_i[I]
        Aap[I,I+1] = -DENSITY * dd[I+1] * area_i[I+1]
        bb[I] = DENSITY * (area_i[I]*psd_u[I] - area_i[I+1]*psd_u[I+1])
        Aap[I,I] = -(Aap[I,I-1] + Aap[I,I+1])
    for i in [0, P_NODES_NUM-1]:
        Aap[i,i] = 1
        bb[i] = 0
    i = 0
    Aap[i, i] = 1
    bb[i] = P_INLET - 0.5 * DENSITY * psd_u[i+1]**2 * (area_i[i+1] / area_I[i])**2
    
    i = P_NODES_NUM - 1
    Aap[i, i] = 1
    bb[i] = 0

    return np.linalg.solve(Aap, bb)

def sol_moment_equ():
    i = 1
    fw = DENSITY * uu[i] * area_i[i]
    fe = DENSITY * area_I[i] * (uu[i] + uu[i+1])/2
    ai = fe + 0.5*fw*(area_i[i]/AREA_A)**2
    su = (P_INLET - pp[i]) * area_i[i] + fw * area_i[i] / AREA_A * uu[i]
    uu_residual[i] = abs(ai*uu[i] - su)
    uu[i] = su / ai
    dd[i] = area_i[i]/ai

    for i in range(2, U_NODES_NUM):
        aw = DENSITY * area_I[i-1] * (uu[i] + uu0[i-1]) / 2

        if i == U_NODES_NUM-1:
            ai = DENSITY * area_i[i] * uu[i]
        else:
            ai = DENSITY * area_I[i] * (uu[i] + uu[i+1])/2

        dd[i] = area_i[i] / ai
        uu_residual[i] = abs(ai * uu[i] - aw * uu[i-1] - area_i[i] * (pp[i-1] - pp[i]))
        uu[i] = aw/ai * uu[i-1] + dd[i] * (pp[i-1] - pp[i]) 


def sol_pressure_cor_equ():
    for I in range(1,P_NODES_NUM-1):
        Aap[I,I-1] = -DENSITY * dd[I] * area_i[I]
        Aap[I,I+1] = -DENSITY * dd[I+1] * area_i[I+1]
        bb[I] = DENSITY * (area_i[I]*uu[I] - area_i[I+1]*uu[I+1])
        Aap[I,I] = -(Aap[I,I-1] + Aap[I,I+1])
    for i in [0, P_NODES_NUM-1]:
        Aap[i,i] = 1
        bb[i] = 0
    return np.linalg.solve(Aap, bb)


def correction_uvp(pp_cor):
    for i in range(1,U_NODES_NUM):
        uu[i] = uu[i] + ALPHA*dd[i]*(pp_cor[i-1] - pp_cor[i])
    global pp, uu0
    pp = pp + ALPHA*pp_cor
    pp[0] = P_INLET - 0.5*DENSITY*uu[1]**2 * (area_i[1]/AREA_A)
    uu0 = uu.copy()


def residual(j):
    if j % N_PRINT == 0:
        print('i = ',j, '   u_residul = ', uu_residual.sum())
        return False
    if uu_residual.sum() < 1e-5:
        print('i = ',j, '   u_residul = ', uu_residual.sum())
        return True
    if uu_residual.sum() > 1e50:
        print('\ni = ',j, '   u_residul = ', uu_residual.sum(),' @(&_&)@  !!!!!!!!')
        print('Please reduce the under-relaxation factor (alpha) !\n')
        sys.exit()
        return False

#def sol_transport_equ():


def exact_solution():
    nodes_excat = 50
    area_exact = np.linspace(AREA_A, AREA_E, nodes_excat, endpoint=True)
    xe_grid = np.linspace(0, LENGTH, nodes_excat, endpoint=True)
    pp_exact = np.zeros(nodes_excat)
    uu_exact = np.zeros(nodes_excat)
    for i in range(nodes_excat):
        pp_exact[i] = P_INLET - 0.5 * 0.44721**2 /DENSITY / area_exact[i]**2
        uu_exact[i] = math.sqrt((P_INLET - pp_exact[i])*2/DENSITY)
    return uu_exact, pp_exact, xe_grid

def show_result():
    print('mass flow rate = ', DENSITY*uu[-1]*area_i[-1])
    exact_u, exact_p ,xe_grid = exact_solution()

    plt.figure(1)
    plt.xlabel('Distance (m)')
    plt.ylabel('Velocity (m/s)')   # ('Pressure (Pa)')
    plt.plot(xu_grid[1:], uu[1:] ,'bs--', label='Numerical')
    plt.plot(xe_grid, exact_u,'k', label='Exact')
    plt.title('Velocity')  # ('Pressure')
    plt.legend()
    
    plt.figure(2)
    plt.xlabel('Distance (m)')
    plt.ylabel('Pressure (Pa)')
    plt.plot(xp_grid, pp ,'bs--', label='Numerical')
    plt.plot(xe_grid, exact_p,'k', label='Exact')
    plt.title('Pressure')
    plt.legend()
    plt.show()


def simple():
    grid()
    initialization()

    for j in range(MAX_ITERATION):
        sol_moment_equ()
        p_correction = sol_pressure_cor_equ()
        correction_uvp(p_correction)
        if residual(j): break

    show_result()

    
def simpleR():
    grid()
    initialization()

    for j in range(300):
        sol_pseudo_velocity()
        global pp
       # print(pp)
        pp = sol_pressure_equ()
       # print(pp)
      #  print()
        sol_moment_equ()
        p_correction = sol_pressure_cor_equ()
        correction_uvp(p_correction)
        if residual(j): break

    show_result()

#def simpleC():
#def piso():


def main():
    #simple()
    simpleR()


#=== parameters ====================
U_NODES_NUM = 10
P_NODES_NUM = U_NODES_NUM
ALPHA = 0.01
MAX_ITERATION = 15001

LENGTH = 2
DENSITY = 1
AREA_A = 0.5
AREA_E = 0.1
P_INLET = 10
P_ENLET = 0

N_PRINT = 10

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
uu0 = uu.copy()
psd_u = np.zeros(U_NODES_NUM)
pp = np.linspace(P_INLET, P_ENLET, P_NODES_NUM, endpoint=True)
uu_residual = np.zeros(U_NODES_NUM)
#===================================================================

if __name__ == '__main__':
    main()
    