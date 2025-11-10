import numpy as np
from numpy.random import *
from copy import deepcopy

seed(10)

N = 1e5
Np = int(1e7) # 1e12
P = N/Np

dt = 1e-6 # 7.3e-9
N_det = int(1/dt)
N_t= int(Np*dt)

foto_init = np.array(rand(Np)<=P,dtype=int)
G = deepcopy(foto_init)
T = deepcopy(foto_init)
R = deepcopy(foto_init)

indx = foto_init.nonzero()[0]

BS= (rand(len(indx))<=0.5)


T[indx] = BS
R[indx] = ~BS

G_det = []
T_det = []
R_det = []

for n in range(N_det):
	if np.any(G[n*N_t:n*N_t+N_t]>0):
		G_det.append(1)
	else:
		G_det.append(0)


	if  np.any(T[n*N_t:n*N_t+N_t]>0):
		T_det.append(1)
	else:
		T_det.append(0)


	if  np.any(R[n*N_t:n*N_t+N_t]>0):
		R_det.append(1) 
	else:
		R_det.append(0)

indx_G = np.array(G_det).nonzero()[0]
indx_T = np.array(T_det).nonzero()[0]
indx_R = np.array(R_det).nonzero()[0]

N_G = np.sum(G_det) 
N_T = np.sum(T_det)
N_R = np.sum(R_det)

N_TR = len(set(indx_T) & set(indx_R))

N_GT = len(set(indx_G) & set(indx_T))
N_GR =len(set(indx_G) & set(indx_R))
N_GTR = len((set(indx_G) & set(indx_R))& set(indx_T))

g2_3D = (N_GTR*N_G)/(N_GT*N_GR)

g2_2D = (N_TR*N_det)/(N_T*N_R)

#print(set(indx_T) & set(indx_R),'\n \n \n ',(set(indx_G) & set(indx_R))& set(indx_T))
print(g2_2D ,g2_3D)


