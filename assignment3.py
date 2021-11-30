import pandas as pd
from grid import grid
from montecarlo import montecarlo
from neldermead import neldermead
from gridbrute import gridbrute

#For the loss function you want to use, choose between 'least squares', 'least log-squares', 'relative error'
#Additionaly, as first argument, 
#grid() expects the number of discretization points on each axis
#montecarlo() expects the sample size
#gridbrute uses the function 'brute' from scipy.optimize and performs grid search as well

file = open('optimals.txt', 'w')

objectives = ['least squares', 'least log-squares', 'relative error']

#This loop takes some time
for method in objectives:
    #Too many lines, file.write doesn't accept many arguments... maybe we find another way to write this

    file.write(method)
    file.write('\n')
    file.write('\n')

    res1 = grid(30, method)
    res2 = montecarlo(50, method)
    res3 = neldermead(method)
    res4 = gridbrute(method)

    file.write('Grid Search:')
    file.write('\n')
    file.write('[n D q Lambda] =')
    file.write(str(res1[0]))
    file.write('\n')
    file.write('Minimum value of the objective function: ')
    file.write('\n')
    file.write(str(res1[1]))
    file.write('\n')
     
    file.write('Monte Carlo:')
    file.write('\n')
    file.write('[n D q Lambda] =')
    file.write(str(res2[0]))
    file.write('\n')
    file.write('Minimum value of the objective function: ')
    file.write('\n')
    file.write(str(res2[1]))
    file.write('\n')
    

    file.write('Nelder Mead:')
    file.write('\n')
    file.write('[n D q Lambda] =')
    file.write(str(res3[0]))
    file.write('\n')
    file.write('Minimum value of the objective function: ')
    file.write('\n')
    file.write(str(res3[1]))
    file.write('\n')


    file.write('Built-in Grid Search:')
    file.write('\n')
    file.write('[n D q Lambda] =')
    file.write(str(res4[0]))
    file.write('\n')
    file.write('Minimum value of the objective function: ')
    file.write('\n')
    file.write(str(res4[1]))
    file.write('\n')
    file.write('\n')

file.close()