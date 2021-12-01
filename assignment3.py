import numpy as np
import pandas as pd
from grid import grid
from montecarlo import montecarlo
from neldermead import neldermead
from gridbrute import gridbrute
import matplotlib.pyplot as plt

#For the loss function you want to use, choose between 'least squares', 'least log-squares', 'relative error'
#Additionaly, as first argument, 
#grid() expects the number of discretization points on each axis
#montecarlo() expects the sample size
#gridbrute uses the function 'brute' from scipy.optimize and performs grid search as well


fig, axes = plt.subplots(nrows = 2, ncols = 2)
axisLabelFontSize = 10
labels = ['n', 'D', 'q', 'Lambda']
l = np.arange(len(labels))

file = open('optimals.txt', 'w')

objectives = ['least squares', 'least log-squares', 'relative error']

#This loop takes some time
for j, objective in enumerate(objectives):
    #Too many lines, file.write accepts only one argument... maybe we find another way to write this

    file.write(objective)
    file.write('\n')
    file.write('\n')

    res = [grid(30, objective), montecarlo(50, objective), neldermead(objective), gridbrute(objective)]

    fig, axes = plt.subplots(nrows = 2, ncols = 2)

    for i, method in enumerate(['Grid Search', 'Monte Carlo', 'Nelder Mead', 'Built-in Grid Search']):
    
        file.write(method)
        file.write('\n')
        file.write('[n D q Lambda] =')
        file.write(str(res[i][0]))
        file.write('\n')
        file.write('Minimum value of the objective function: ')
        file.write('\n')
        file.write(str(res[i][1]))
        file.write('\n')

        axes[int(i / 2), i % 2].bar(l, res[i][0])
        axes[int(i / 2), i % 2].set_title(method, fontsize = axisLabelFontSize)
        axes[int(i / 2), i % 2].set_xticks(l)
        axes[int(i / 2), i % 2].set_xticklabels(labels)
        axes[int(i / 2), i % 2].legend()
        fig.suptitle(objective)

    plt.show()

    file.write('\n')

file.close()



