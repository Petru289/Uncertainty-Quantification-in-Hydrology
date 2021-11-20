import numpy as np
import chaospy as cp
from numpy.core.fromnumeric import size
from functions import c

import parameters as par


x = par.x
t = par.t

file = open('output.txt', 'w')

for N in [1000, 10000, 100000]: #number of realizations

    print("Number of realizations =", N)
    print()
    file.write("Number of realizations = ")
    file.write(str(N))
    file.write('\n')
    file.write('\n')
    file.write('\n')

    for distIndex, distance in enumerate(x):
        for timeIndex, timestep in enumerate(t):
            
            # Get Matrix A 
            n = cp.Uniform(par.nRange[0], par.nRange[1])
            D = cp.Uniform(par.DRange[0], par.DRange[1])
            q = cp.Uniform(par.qRange[0], par.qRange[1])
            Lambda = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
            jointdist = cp.J(n, D, q, Lambda)
            A = np.transpose(jointdist.sample(size = N))
            #transpose, so that it's the same as in the algorithm description

            #print(A.shape)
            #print(np.max(A, axis = 0))
            #print(np.min(A, axis = 0))


            # Get matrix B
            n2 = cp.Uniform(par.nRange[0], par.nRange[1])
            D2 = cp.Uniform(par.DRange[0], par.DRange[1])
            q2 = cp.Uniform(par.qRange[0], par.qRange[1])
            Lambda2 = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
            jointdist2 = cp.J(n2, D2, q2, Lambda2)
            B = np.transpose(jointdist2.sample(size = N))

            #print(B.shape)
            #print(np.max(B, axis = 0))
            #print(np.min(B, axis = 0))

            # Get matrices BA
            #https://stackoverflow.com/questions/3059395/numpy-array-assignment-problem

            BAn = B.copy()
            BAD = B.copy()
            BAq = B.copy()
            BALambda = B.copy()

            BAn[:,0] = A[:,0]
            BAD[:,1] = A[:,1]
            BAq[:,2] = A[:,2]
            BALambda[:,3] = A[:,3]

            #print(B)
            #print(BAn)

            #print(B)
            #print(BAD)

            # Get the matrices AB
            ABn = A.copy()
            ABD = A.copy()
            ABq = A.copy()
            ABLambda = A.copy()

            ABn[:,0] = B[:,0]
            ABD[:,1] = B[:,1]
            ABq[:,2] = B[:,2]
            ABLambda[:,3] = B[:,3]
            
            #print(A)
            #print(ABn)

            # Get matrix C
            C = np.append(A, B, axis = 0)

            # Get concentrations
            yA = c(distance, timestep, par.M, A[:,0], A[:,1], A[:,2], A[:,3])
            yB = c(distance, timestep, par.M, B[:,0], B[:,1], B[:,2], B[:,3])
        

            yABn = c(distance, timestep, par.M, ABn[:,0], ABn[:,1], ABn[:,2], ABn[:,3])
            yABD = c(distance, timestep, par.M, ABD[:,0], ABD[:,1], ABD[:,2], ABD[:,3])
            yABq = c(distance, timestep, par.M, ABq[:,0], ABq[:,1], ABq[:,2], ABq[:,3])
            yABLambda = c(distance, timestep, par.M, ABLambda[:,0], ABLambda[:,1], ABLambda[:,2], ABLambda[:,3])


            yBAn = c(distance, timestep, par.M, BAn[:,0], BAn[:,1], BAn[:,2], BAn[:,3])
            yBAD = c(distance, timestep, par.M, BAD[:,0], BAD[:,1], BAD[:,2], BAD[:,3])
            yBAq = c(distance, timestep, par.M, BAq[:,0], BAq[:,1], BAq[:,2], BAq[:,3])
            yBALambda = c(distance, timestep, par.M, BALambda[:,0], BALambda[:,1], BALambda[:,2], BALambda[:,3])

            yC = c(distance, timestep, par.M, C[:,0], C[:,1], C[:,2], C[:,3])
            

            #We need nicer notations


            print("x = ", distance, "time = ", timestep)
            file.write("x = ")
            file.write(str(distance))
            file.write(" m")
            file.write(" time = ")
            file.write(str(timestep))
            file.write(" days")
            file.write('\n')
            file.write('\n')

            # First order sobol
            file.write("First order Sobol indices")
            file.write('\n')
            print("First order Sobol indices")
            Sn =  ((np.inner(yA, yBAn) - np.inner(yA, yB)) / N ) / (np.inner(yC, yC) / (2 * N) - np.mean(yC) ** 2)
            file.write("n: ")
            file.write(str(Sn))
            file.write('\n')
            print("n: ", Sn)

            SD = ((np.inner(yA, yBAD) - np.inner(yA, yB)) / N ) / (np.inner(yC, yC) / (2 * N) - np.mean(yC) ** 2)
            file.write("D: ")
            file.write(str(SD))
            file.write('\n')
            print("D: ", SD)

            Sq = ((np.inner(yA, yBAq) - np.inner(yA, yB)) / N ) / (np.inner(yC, yC) / (2 * N) - np.mean(yC) ** 2)
            file.write("q: ")
            file.write(str(Sq))
            file.write('\n')
            print("q: ", Sq)

            SLambda = ((np.inner(yA, yBALambda) - np.inner(yA, yB)) / N ) / (np.inner(yC, yC) / (2 * N) - np.mean(yC) ** 2)
            file.write("Lambda: ")
            file.write(str(SLambda))
            file.write('\n')
            file.write('\n')
            print("Lambda: ", SLambda)

            file.write("Total Sobol indices")
            file.write('\n')
            print("Total Sobol indices")
            STn = (np.sum((yA - yABn) ** 2) / (2 * N)) / (np.inner(yC, yC) / (2 * N) - np.mean(yC) ** 2)
            file.write("n: ")
            file.write(str(STn))
            file.write("\n")
            print("n: ", STn)

            STD = (np.sum((yA - yABD) ** 2) / (2 * N)) / (np.inner(yC, yC) / (2 * N) - np.mean(yC) ** 2)
            file.write("D: ")
            file.write(str(STD))
            file.write("\n")
            print("D: ", STD)

            STq = (np.sum((yA - yABq) ** 2) / (2 * N)) / (np.inner(yC, yC) / (2 * N) - np.mean(yC) ** 2)
            file.write("q: ")
            file.write(str(STq))
            file.write("\n")
            print("q: ", STq)

            STLambda = (np.sum((yA - yABLambda) ** 2) / (2 * N)) / (np.inner(yC, yC) / (2 * N) - np.mean(yC) ** 2)
            file.write("LAmbda: ")
            file.write(str(STLambda))
            file.write("\n")
            file.write("\n")
            print("Lambda: ", STLambda)

            print()

            #Reference: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis


file.close()

