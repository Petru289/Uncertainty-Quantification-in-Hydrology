import numpy as np
from grid import grid
from montecarlo import montecarlo
from neldermead import neldermead
import matplotlib.pyplot as plt

#For the loss function you want to use, choose between 'least squares', 'least log-squares', 'relative error'
#Additionaly, as first argument, 
#grid() expects the number of discretization points on each axis
#montecarlo() expects the sample size
#gridbrute uses the function 'brute' from scipy.optimize and performs grid search as well

# Definitions of plotting parameters and labels
axisLabelFontSize = 12
titleFontSize = 20
methodLabels = ['Grid-Search', 'Monte Carlo', 'Nelder-Mead']
parameterLabels = ['Porosity [$n$]', 'Dispersion Coefficient [$D$]', 'Specific Discharge [$q$]', 'Decay Rate [$\lambda$]']
yAxisLabels = ['$[-]$', '$m^2/d$', '$m/d$', '$1/d$']
l = np.arange(len(methodLabels))
width = 0.2
objectives = ['least squares', 'least log-squares', 'relative error']
fig, axes = plt.subplots(nrows = 2, ncols = 2)
optimals = np.zeros((4,len(methodLabels),len(objectives)))

file = open('optimals.txt', 'w')

# Calculate optimal values with different algorithms and objective functions
for j, objective in enumerate(objectives):

    file.write(objective + '\n\n')
    res = [grid(30, objective), montecarlo(50, objective), neldermead(objective)]

    for i, method in enumerate(['Grid Search', 'Monte Carlo', 'Nelder Mead']):
        
        string = method + '\n[n D q Lambda] =' + str(res[i][0])
        string += '\nMinimum value of the objective function:\n' + str(res[i][1]) + '\n'
        file.write(string)

        optimals[0][i][j] = res[i][0][0]
        optimals[1][i][j] = res[i][0][1]
        optimals[2][i][j] = res[i][0][2]
        optimals[3][i][j] = res[i][0][3]

    file.write('\n')

file.close()

# Create bar plots  
for index, ax in enumerate(fig.axes):  
    bar1 = ax.bar(l- 1.3*width, optimals[index, :, 0], width, label=objectives[0])
    bar2 = ax.bar(l, optimals[index, :, 1], width, label=objectives[1])
    bar3 = ax.bar(l + 1.3*width, optimals[index, :, 2], width, label=objectives[2])
    ax.set_title(parameterLabels[index], fontsize = titleFontSize)
    ax.set_xticks(l)
    ax.set_xticklabels(methodLabels, fontsize=axisLabelFontSize)
    ax.set_ylabel(yAxisLabels[index])
    ax.legend()

    ax.bar_label(bar1, padding=3, fmt='%.3f')
    ax.bar_label(bar2, padding=3, fmt='%.3f')
    ax.bar_label(bar3, padding=3, fmt='%.3f')

axes[0,0].set_ylim([0, 0.6])
axes[0,1].set_ylim([0, 1.1])
axes[1,0].set_ylim([0, 0.75])
axes[1,1].set_ylim([0, 0.03])

plt.tight_layout()
plt.show()