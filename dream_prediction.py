from turtle import color
import numpy as np
from functions import c
from dream_results import param_all
import matplotlib.pyplot as plt


#plt.cm.summer(np.linspace(0, 1, no_plots))
time = np.arange(1,101)
no_plots = len(time)
param = param_all[-1000:]
param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))
for i in range(param.shape[0]):
    concentration = c(100, time, 200, param[i,0], param[i,1], param[i,2], param[i,3])
    plt.plot(time, concentration, alpha = 0.2, color = 'navy')

plt.ylabel('concentration $[kg/m^3]$')
plt.xlabel('time [days]')
plt.show()
