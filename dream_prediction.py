from cProfile import label
from turtle import color
import numpy as np
from functions import c
from dream_results import param_all, results
import matplotlib.pyplot as plt
from spotpy.analyser import get_best_parameterset


#plt.cm.summer(np.linspace(0, 1, no_plots))
time = np.arange(1,101)
no_plots = len(time)
param = param_all[-1000:]
param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))
mean_param = np.mean(param,0)

concentrations = []
concentration_plots = []
for i in range(param.shape[0]):
    concentration = c(100, time, 200, param[i,0], param[i,1], param[i,2], param[i,3])
    plot, = plt.plot(time, concentration, alpha = 0.2, color = 'navy')
    concentration_plots.append(plot)
    concentrations.append(concentration)
concentration_plots[-1].set_label('Realisations of posterior distribution')

# Plot mean concentration
mean_concentration = np.mean(concentrations,0)
mean_plot, = plt.plot(time,mean_concentration, color = 'firebrick',linewidth=2, label='Mean concentration')

# Plot concentration using the parameter means
concentration = c(100, time, 200, mean_param[0],mean_param[1],mean_param[2],mean_param[3])
plt.plot(time,concentration, color = 'green', linewidth=2, label='Concentration using parameter means')
plt.legend()

# Plot concentration using the optimal parameter set
optimal_parameters = get_best_parameterset(results)[0]
concentration = c(100, time, 200, optimal_parameters[0],optimal_parameters[1],optimal_parameters[2],optimal_parameters[3])
plt.plot(time,concentration, color = 'orange', linewidth=2, label='Concentration using the optimal parameter set')
plt.legend()

# Plot threshold concentration
plt.axhline(y=2.5, color='r', linestyle='--', alpha=0.5, linewidth=0.5)

print(get_best_parameterset(results)[0])

plt.ylabel('concentration $[kg/m^3]$', fontsize=14)
plt.xlabel('time [days]', fontsize=14)
plt.suptitle('Prediction of the concentration at x = 100 m using the posterior distribution of the Dream-Algorithm', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend()
plt.show()
