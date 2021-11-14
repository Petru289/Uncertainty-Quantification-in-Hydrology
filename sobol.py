import numpy as np
import chaospy as cp
from numpy.core.fromnumeric import size
from functions import c

import parameters as par

numberOfRealizations = 1000
distance = par.x[0]
timestep = par.t[0]

# Get Matrix A 
n = cp.Uniform(par.nRange[0], par.nRange[1])
D = cp.Uniform(par.DRange[0], par.DRange[1])
q = cp.Uniform(par.qRange[0], par.qRange[1])
Lambda = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
jointdist = cp.J(n, D, q, Lambda)
A = jointdist.sample(size = numberOfRealizations)

# Get matrix B
n2 = cp.Uniform(par.nRange[0], par.nRange[1])
D2 = cp.Uniform(par.DRange[0], par.DRange[1])
q2 = cp.Uniform(par.qRange[0], par.qRange[1])
Lambda2 = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
jointdist2 = cp.J(n2, D2, q2, Lambda2)
B = jointdist2.sample(size = numberOfRealizations)

print(type(B))

# Get matrices BA
BAn = B.copy()
BAD = B.copy()
BAq = B.copy()
BALambda = B.copy()

BAn[0,:] = A[0,:]
BAD[1,:] = A[1,:]
BAq[2,:] = A[2,:]
BALambda[3,:] = A[3,:]

# Get the matrices AB
ABn = A.copy()
ABD = A.copy()
ABq = A.copy()
ABLambda = A.copy()

ABn[0,:] = B[0,:]
ABD[1,:] = B[1,:]
ABq[2,:] = B[2,:]
ABLambda[3,:] = B[3,:]

# Get matrix C
C = np.append(A.copy(),B.copy(),axis=1)

# Get concentrations
yA = c(distance, timestep, par.M, A[0,:], A[1,:], A[2,:], A[3,:])
yB = c(distance, timestep, par.M, B[0,:], B[1,:], B[2,:], B[3,:])

yABn = c(distance, timestep, par.M, ABn[0,:], ABn[1,:], ABn[2,:], ABn[3,:])
yABD = c(distance, timestep, par.M, ABD[0,:], ABD[1,:], ABD[2,:], ABD[3,:])
yABq = c(distance, timestep, par.M, ABq[0,:], ABq[1,:], ABq[2,:], ABq[3,:])
yABLambda = c(distance, timestep, par.M, ABLambda[0,:], ABLambda[1,:], ABLambda[2,:], ABLambda[3,:])

yBAn = c(distance, timestep, par.M, BAn[0,:], BAn[1,:], BAn[2,:], BAn[3,:])
yBAD = c(distance, timestep, par.M, BAD[0,:], BAD[1,:], BAD[2,:], BAD[3,:])
yBAq = c(distance, timestep, par.M, BAq[0,:], BAq[1,:], BAq[2,:], BAq[3,:])
yBALambda = c(distance, timestep, par.M, BALambda[0,:], BALambda[1,:], BALambda[2,:], BALambda[3,:])

yC = c(distance, timestep, par.M, C[0,:], C[1,:], C[2,:], C[3,:])

# First order sobol
Sn = ((np.sum(np.multiply(yA,yBAn) - np.multiply(yA, yB))) / numberOfRealizations)/((np.sum(np.multiply(yC, yC) - np.mean(yC)**2))/ 2 * numberOfRealizations)
print(Sn)