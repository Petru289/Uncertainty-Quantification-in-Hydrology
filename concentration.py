import matplotlib.pyplot as plt
import numpy as np
import csv
from functions import c
import parameters as par

filePath = './optimals.csv'

def readOptimalValues(filePath):
    
    # 3 methods, 4 parameters, 3 objective functions
    optimalValues = np.zeros([3,4,3])
    objIndex = -1
    methodIndex = 0

    with open(filePath) as csvFile:
        reader = csv.reader(csvFile)

        for line in reader:

            if line[0].startswith('Objective'):
                objIndex += 1
            elif line[0].startswith('method'):
                pass
            else:
                optimalValues[methodIndex%3,:,objIndex] = line[1:]
                methodIndex += 1

    return optimalValues

#fig, axes = plt.subplots(nrows=1, ncols=3)

optimalValues = readOptimalValues(filePath)

time = np.arange(1,101)

obj = ['Least Squares', 'Least Log-Squares', 'Relative Error']
methods = ['Grid Search', 'Monte Carlo', 'Nelder-Mead']

for j in range(3):
    plt.figure(j+1)
    ax = plt.gca()
    ax.set_title('Prediction of the concentration at 100 m\nobjective function: ' + obj[j])
    ax.set_ylabel('concentration $[kg/m^3]$')
    ax.set_xlabel('days')

    peaks = []
    for i in range(3):
        concentration = c(100, time, par.M, optimalValues[i,0,j],optimalValues[i,1,j],optimalValues[i,2,j],optimalValues[i,3,j])
        peaks.append([concentration.max(), concentration.argmax()])
        ax.plot(time, concentration, label=methods[i])

    textString = '\n'.join((
        '$\t{Peaks:}$',
        'Grid-Search: %.4f at day %.0f' %(peaks[0][0], peaks[0][1]),
        'Monte Carlo: %.4f at day %.0f' %(peaks[1][0], peaks[1][1]),
        'Nelder-Mead: %.4f at day %.0f' %(peaks[2][0], peaks[2][1])
    ))

    props = dict(boxstyle='round', facecolor='white')

    # place a text box in upper left in axes coords
    ax.text(0.02, 0.97, textString, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)
    
    ax.legend()

plt.tight_layout()
plt.show()


