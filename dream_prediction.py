<<<<<<< HEAD
import py_compile
from turtle import color
import numpy as np
from functions import c
import matplotlib.pyplot as plt
import spotpy.analyser
=======
from cProfile import label
from turtle import color
import numpy as np
from functions import c
from dream_results import param_all, results
import matplotlib.pyplot as plt
from spotpy.analyser import get_best_parameterset
>>>>>>> 8c52e2142450822e8c1570bab70fe28a1e004289


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

print("The probability that the concentration exceeds the critical level at 100 meters is: ", count/param.shape[0])
print(" and this event is expected on day ", np.floor(np.mean(tcrit)))
print("The earliest expected time that the contaminant concentration can exceed the critical level is: ", np.min(tcrit))