import numpy as np
import matplotlib.pyplot as plt
import math
import sys

#== parameters ===
u_nodes_num = p_nodes_num =  5  # number of mesh cells or nodes
alpha = 0.8  # under-relaxation factor
nmax = 2000  # max number of iteration
residual_criteria = 1e-6

length = 2
density = 1
area_A = 0.5
area_E = 0.1
p_0 = 10  # inlet pressure
p_e = 0   # outlet pressure
#===================
# ==== creating mesh ======
xp_grid = np.linspace(0,length,p_nodes_num, endpoint=True) # pressure mesh
print('xp_grid = ',xp_grid)

xu_grid = np.zeros(u_nodes_num) # velocity mesh
for i in range(1,u_nodes_num):
    xu_grid[i] = (xp_grid[i] + xp_grid[i-1])/2
print('xu_grid = ',xu_grid)
#===== preparing constant values ======
area_I = np.linspace(area_A, area_E, p_nodes_num, endpoint=True) 
# area of sections A,B,C,D,E  ( Pressure nodes )

area_i = np.zeros(u_nodes_num) 
for i in range(1,u_nodes_num):
    area_i[i] = (area_I[i] + area_I[i-1])/2
# area of sections of 1,2,3,4 ( Velocity nodes )

# === coefficients =====
au = np.zeros(u_nodes_num)
dd = np.zeros(u_nodes_num)
bb = np.zeros(p_nodes_num)
Aap = np.zeros((p_nodes_num, p_nodes_num))
Bbp = np.zeros(p_nodes_num)


'''
# === display area  ===
plt.xlabel('Distance (m)')
plt.ylabel('A')
plt.plot(xu_grid[1:], area_i[1:] ,'bs', label='area_i')
plt.plot(xp_grid, area_I,'go-', label='area_I')
#plt.title('')
plt.legend()
plt.show()
sys.exit()
'''


#== initialization ==
guessedMassFlowRate = 1  # mass flow rate

uu = np.zeros(u_nodes_num)   # velocity , the index starts with 1
uu[1:] = guessedMassFlowRate/density/area_i[1:]   # initial velocity

uu0 = uu.copy()  # 
uu_residual = np.zeros(u_nodes_num) # residual of velcity

pp = np.linspace(p_0, p_e, p_nodes_num, endpoint=True) # initial pressure

print('initial u = ',uu)
print('initial p = ',pp)
# ==================

#== SIMPLE ALGORITHM ==

for j in range(nmax):
    ##== STEP 1 : solving momentum equations ==
    #== inlet ==
    i = 1
    fw = density * uu[i] * area_i[i]  
    fe = density * area_I[i] * (uu[i] + uu[i+1])/2
    ai = fe + 0.5*fw*(area_i[i]/area_A)**2          #ap
    su = (p_0 - pp[i])*area_i[i] + fw * area_i[i]/area_A * uu[i]  # source term
    uu_residual[i] = abs(ai*uu[i] - su)
    uu[i] = su / ai
    dd[i] = area_i[i]/ai
    #print(uu[i], uu0[i])

    #== internal nodes ==
    for i in range(2,u_nodes_num):
        aw = density * area_I[i-1] * (uu[i]+uu0[i-1])/2   # aw = Fw = (rho * u * A)_w
                                                          # ae = 0
        if i == u_nodes_num-1:
            ai = density * area_i[i] * uu[i]  # outlet
        else:
            ai = density * area_I[i] * (uu[i]+uu[i+1])/2   # ap = aw + ae + (Fe - Fw) 
                                                           #    = Fw + Fe - Fw 
                                                           #    = Fe
                                                           #    = (rho * u * A)_e 

        dd[i] = area_i[i]/ai  # d = a_av / ap
        uu_residual[i] = abs(ai*uu[i] - aw*uu[i-1] - area_i[i]*(pp[i-1]-pp[i]) )
            # according to formula (3), assume ap * up = aw * uw + ae * ue + su + residual. 
            # so, residual = ap * up - aw * uw - ae * ue - su 

        uu[i] = aw/ai * uu[i-1] + dd[i] * (pp[i-1]-pp[i])  # up = aw/ap * uw + su/ap
                                                           #    = aw/ap * uw + d * Delta p

    #print(uu[1:])

    ##== STEP 2: solving pressure correction equation ==
    for I in range(1,p_nodes_num-1):
        Aap[I,I-1] = -density * dd[I] * area_i[I]
        Aap[I,I+1] = -density * dd[I+1] * area_i[I+1]
        bb[I] = density * (area_i[I]*uu[I] - area_i[I+1]*uu[I+1])
        Aap[I,I] = -(Aap[I,I-1] + Aap[I,I+1])

    for i in [0, p_nodes_num-1]:  # inlet and outlet
        Aap[i,i] = 1
        bb[i] = 0

    pp_correction = np.linalg.solve(Aap, bb)  # solve pressure correction equations
    #print(pp_correction)

    ##== STEP 3: correct presure and velocity ==
    for i in range(1,u_nodes_num):
        uu[i] = uu[i] + alpha *dd[i]*(pp_correction[i-1] - pp_correction[i]) 

    pp = pp + alpha *pp_correction
    pp[0] = p_0 - 0.5*density*uu[1]**2 * (area_i[1]/area_A)

    uu0 = uu.copy()
    if j%100==0:
        print('i = ',j, '   u_residul = ', uu_residual.sum())

    if uu_residual.sum() < residual_criteria:
        print('i = ',j, '   u_residul = ', uu_residual.sum())
        break
    if uu_residual.sum() > 1e10:
        print('\ni = ',j, '   u_residul = ', uu_residual.sum(),' @(&_&)@  !!!!!!!!')
        print('Please reduce the under-relaxation factor (alpha) !\n')
        sys.exit()
        

print('mass flow rate = ', density*uu[-1]*area_i[-1])

##== exact solution ==
nodes_excat = 50
area_exact = np.linspace(area_A, area_E, nodes_excat, endpoint=True)
xe_grid = np.linspace(0, length, nodes_excat, endpoint=True)
pp_exact = np.zeros(nodes_excat)
uu_exact = np.zeros(nodes_excat)
for i in range(nodes_excat):
    pp_exact[i] = p_0 - 0.5 * 0.44721**2 /density / area_exact[i]**2
    uu_exact[i] = math.sqrt((p_0 - pp_exact[i])*2/density)


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
plt.plot(xe_grid, pp_exact,'k', label='Exact')
plt.title('Pressure')
plt.legend()

plt.show()
