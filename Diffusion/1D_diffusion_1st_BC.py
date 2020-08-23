'''
有限体积法求解一维扩散问题（第一类边界条件）
离散方程推导见博文：
https://blog.csdn.net/weixin_42562856/article/details/106799542
'''

import numpy as np 
import matplotlib.pyplot as plt

#== parameters ===
nx = 5  # 网格单元数
nndoes = nx + 2 # 节点数，含边界

#===  for example 1 ==========
''' 
L = 0.5 # 长度，m
Area = 0.01 #截面面积 ,m2
k = 1000  #导热系数 , W/m-k
Ta = 100 # 边界A的温度值 ,C
Tb = 500 # 边界B的温度值 ,C
q = 0 # kW/m3, 体积热源
''' 
# =========================

#***** for example 2 *********
#'''
L = 0.02 # 长度，m
Area = 1 #截面面积 ,m2
k = 0.5  #导热系数 , W/m-k
Ta = 100 # 边界A的温度值 ,C
Tb = 200 # 边界B的温度值 ,C
q = 1e6 # W/m3, 体积热源
#'''
#****************************


#==  x grid ==
dx = L/nx  # 网格间距
print('dx = ',dx)
x = np.zeros(nndoes)  # x网格
x[1:nndoes-1] = np.linspace(dx/2, L-dx/2, nx) # 以边界A为原点创建网格点的坐标值
x[-1] = x[-2] + dx/2 #边界B的坐标值
print('x grid = ', x, '\n')

#==  solution array ==
tt = np.zeros(nndoes)  # 解向量
tt[0] = Ta # 边界值
tt[-1] = Tb

#== matrix ==
A = np.zeros((nx, nx)) 
b = np.zeros(nx)

su = q*dx
for i in range(1, nx-1):    # 内部网格节点
    A[i][i-1] = -k*Area/dx
    A[i][i+1] = -k*Area/dx
    A[i][i] = -(A[i][i-1] + A[i][i+1])
    b[i] = su

# for boundary A
i = 0  
A[i][i+1] = -k*Area/dx
su = 2*k*Area*Ta/dx + q*dx
sp = -2*k*Area/dx
A[i][i] = -A[i][i+1] - sp
b[i] = su

# for boundary B
i = nx-1  
A[i][i-1] = -k*Area/dx
su = 2*k*Area*Tb/dx + q*dx
sp = -2*k*Area/dx
A[i][i] = -A[i][i-1] - sp
b[i] = su

print('A = \n', A, '\n')
print('b = \n', np.matrix(b).T ,'\n')

t_temp = np.linalg.solve(A, b)
print('solution = \n', np.matrix(t_temp).T, '\n')
tt[1:nndoes-1] = t_temp

xx = np.linspace(0, L, 50, endpoint=True)
exact_tt = np.zeros(50)
for i in range(50):
    #exact_tt[i] = 800*xx[i] + 100  # 例1 解析解
    exact_tt[i] = ((Tb-Ta)/L + q*(L-xx[i])/(2*k)) * xx[i] + Ta # 例2 解析解
    
plt.xlabel('Distance (cm)')
plt.ylabel('Temperature (C)')
plt.plot(x*100,tt ,'bs', label='Numerical')
plt.plot(xx*100,exact_tt,'k', label='Exact')
plt.legend()
plt.show()
