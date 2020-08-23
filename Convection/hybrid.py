import numpy as np 
import matplotlib.pyplot as plt
import math

#== parameters ===
nx = 25  # 网格单元数
nndoes = nx + 2 # 节点数，含边界

L = 1.0 # 长度，m
Area = 0.01 #截面面积 ,m2
gamma = 0.1  #扩散系数 , kg/m.s
phi_a = 1 # 边界A的温度值 
phi_b = 0 # 边界B的温度值 
rho = 1.0 # 密度， kg/m^3
u = 2.5 # 速度，m/s # 0.1 , 2.5
# =========================


#==  x grid ==
dx = L/nx  # 网格间距
print('dx = ',dx)
x = np.zeros(nndoes)  # x网格
x[1:nndoes-1] = np.linspace(dx/2, L-dx/2, nx) # 以边界A为原点创建网格点的坐标值
x[-1] = x[-2] + dx/2 #边界B的坐标值
print('x grid = ', x, '\n')

#==  solution array ==
phi = np.zeros(nndoes)  # 解向量
phi[0] = phi_a # 边界值
phi[-1] = phi_b

DD = gamma / dx  # D
FF = rho * u     # F
Pe = rho * u * dx / gamma      # Peclet number
#== matrix ==
A = np.zeros((nx, nx)) 
b = np.zeros(nx)

#### 内部网格节点  #########
su = 0.0
sp = 0.0
for i in range(1, nx-1):    
    A[i][i-1] = -max(FF, DD+FF/2, 0)
    A[i][i+1] = -max(-FF, DD-FF/2, 0)
    A[i][i] = -(A[i][i-1] + A[i][i+1]) - sp
    b[i] = su

# for boundary A
i = 0  
A[i][i+1] = -max(-FF, 2*DD-FF/2, 0)
su = (2*DD + FF) * phi_a
sp = -(2*DD + FF)
A[i][i] = -A[i][i+1] - sp
b[i] = su

# for boundary B
i = nx-1  
A[i][i-1] = -max(FF, DD+FF/2, 0)

if Pe>= 2:
    su = 2*DD*phi_b
    sp = -2*DD
else:
    su = (2*DD - FF) * phi_b
    sp = FF - 2*DD

A[i][i] = -A[i][i-1] - sp
b[i] = su
#==========================

print('A = \n', A, '\n')
print('b = \n', np.matrix(b).T ,'\n')

phi_temp = np.linalg.solve(A, b)
print('solution = \n', np.matrix(phi_temp).T, '\n')
phi[1:nndoes-1] = phi_temp

#===== for exact solution ======
xx = np.linspace(0, L, 50, endpoint=True)
exact_solution = np.zeros(50)
for i in range(50):
    exact_solution[i] = (math.exp(rho*u*xx[i] / gamma) -1) / (math.exp(rho*u*L / gamma) -1) * (phi_b - phi_a) + phi_a
    

#UD_solution = np.array([1., 0.99984252, 0.99874016, 0.99212598, 0.95244094, 0.71433071, 0.])
UD_solution = np.array([1.0, 0.9999999867545231, 0.9999999470180921, 0.9999998675452304, 0.999999708599507, 
0.9999993907080597, 0.9999987549251652, 0.9999974833593763, 0.9999949402277984, 0.9999898539646429, 
0.9999796814383317, 0.9999593363857093, 0.9999186462804649, 0.9998372660699758, 0.9996745056489976, 
0.9993489848070412, 0.9986979431231279, 0.9973958597553012, 0.9947916930196478, 0.9895833595483406, 
0.9791666926057266, 0.9583333587204984, 0.9166666909500418, 0.8333333554091289, 0.6666666843273031, 
0.3333333421636515, 0.0])
plt.xlabel('Distance (m)')
plt.ylabel('Phi')
plt.plot(x,phi ,'bs--', label='Numerical (hybrid)')
plt.plot(x,UD_solution ,'go--', label='Numerical (UD)')
plt.plot(xx,exact_solution,'k', label='Exact')
title = 'u= '+str(u)+',  Pe= %.3f'% Pe
plt.title(title.rstrip('0'))
plt.legend()
plt.show()

