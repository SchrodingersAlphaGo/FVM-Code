import numpy as np
import matplotlib.pyplot as plt
import math

#=== parameters ===
x_nodes_number = 10
y_nodes_number = 10
x_length = 1.0
y_length = 1.0
BC_W = 100.0
BC_E = 0.0
BC_N = 100.0
BC_S = 0.0

xx_length = x_length * math.sqrt(2.0)
xy_nodes_number = x_nodes_number * y_nodes_number

#=== grid ===
dx = x_length / x_nodes_number
dy = y_length / y_nodes_number
dxx = xx_length / x_nodes_number
xx = np.linspace(dxx/2, xx_length-dxx/2, x_nodes_number) 

#=== initial ===
phi = np.zeros((x_nodes_number, y_nodes_number))

#=== calculate ===

#=== excat solution ===

#=== show result ===
