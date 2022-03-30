import py_compile
from turtle import color
import numpy as np
from functions import c
import matplotlib.pyplot as plt
import spotpy.analyser


#plt.cm.summer(np.linspace(0, 1, no_plots))
def getXfromConcentration(c, t, M, n, D, q ,Lambda):
    A = M / (np.sqrt(4 * np.pi * D * t))
    x1 = q * t / n - np.sqrt(- 4 * D * t * (np.log(c / A) + Lambda * t))
    x2 = q * t / n + np.sqrt(- 4 * D * t * (np.log(c / A) + Lambda * t))
    return [x1, x2]



results = spotpy.analyser.load_csv_results('HydrologyDREAM')
param_all = results[['parn', 'parD', 'parq', 'parLambda']]
time = np.arange(1,101)
no_plots = len(time)
param = param_all[-1000:]
param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))
param_means = np.array([np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])])
count = 0
tcrit = []
plt.figure(1)
for i in range(param.shape[0]):
    concentration = c(100, time, 200, param[i,0], param[i,1], param[i,2], param[i,3])
    plt.plot(time, concentration, alpha = 0.2, color = 'navy')
    if np.max(concentration) >= 2.5:
        count += 1
    indices = np.where(concentration >= 2.5)
    if len(indices[0]) != 0:
        tcrit.append(np.min(indices))

tcrit = np.array(tcrit)
print(tcrit)

plt.ylabel('concentration $[kg/m^3]$')
plt.xlabel('time [days]')

plt.figure(2)
concentration = c(100, time, 200, param_means[0], param_means[1], param_means[2], param_means[3])
plt.plot(time, concentration, color = 'navy')
plt.title("Mean prediction of concentration at 100 meters using only measurements at well x = 5 meters")

plt.show()

print("The probability that the concentration exceeds the critical level at 100 meters is: ", count/param.shape[0])
print(" and this event is expected on day ", np.floor(np.mean(tcrit)))
print("The earliest expected time that the contaminant concentration can exceed the critical level is: ", np.min(tcrit))