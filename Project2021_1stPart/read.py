import numpy as np
from matplotlib import pyplot as plt


it = 30
file_vort = "Results/w_"+str(it)+".txt"
file_problem = "problem.txt"
problem = np.loadtxt(file_problem)
[h, H, L, Nx, Ny, h_hill, d_hill, sigma_hill, u_hill, u_tau, C] = problem

Nx = int(Nx)
Ny = int(Ny)

vort = np.loadtxt(file_vort)

print(Nx)
#X = np.linspace(0, L, Nx, endpoint = True)
#Y = np.linspace(0, H, Ny, endpoint = True)

X = np.arange(0, L, h+0.001)
Y = np.arange(0, H, h+0.001)

XX, YY = np.meshgrid(X, Y)


plt.figure()
plt.contourf(XX, YY, vort)