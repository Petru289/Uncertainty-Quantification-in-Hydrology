import numpy as np
import spotpy
from dream_setup import spotpy_setup
from dream_run import dream_run
import matplotlib.pyplot as plt
from spotpy.analyser import plot_parameter_trace
from spotpy.analyser import plot_posterior_parameter_histogram



spot_setup = spotpy_setup()

#M, d, k = 4, 2, 4 #https://spotpy.readthedocs.io/en/latest/Sensitivity_analysis_with_FAST/
#rep = (1 + 4 * M ** 2 * (1 + (k - 2) * d)) * k
rep = 1000
sampler = spotpy.algorithms.fast(spot_setup,  dbname='FAST_hydrology',  dbformat='csv')
sampler.sample(repetitions = rep)

results_FAST = spotpy.analyser.load_csv_results('FAST_hydrology')
spotpy.analyser.plot_fast_sensitivity(results_FAST, number_of_sensitiv_pars=4)

r_hat = dream_run()
results = spotpy.analyser.load_csv_results('HydrologyDREAM')
fields = [word for word in results.dtype.names if word.startswith('sim')]
print(results[fields])

fig = plt.figure(figsize=(16,9))
ax = plt.subplot(1,1,1)
q5,q25,q75,q95=[],[],[],[]
for field in fields:
    q5.append(np.percentile(results[field][-100:-1],2.5))
    q95.append(np.percentile(results[field][-100:-1],97.5))
ax.plot(q5,color='dimgrey',linestyle='solid')
ax.plot(q95,color='dimgrey',linestyle='solid')
ax.fill_between(np.arange(0,len(q5),1),list(q5),list(q95),facecolor='dimgrey',zorder=0,
                linewidth=0,label='parameter uncertainty')  
ax.plot(spot_setup.evaluation(),'r.',label='data')
ax.set_ylim(-5,70)
ax.set_xlim(0,15)
ax.legend()
fig.savefig('python_hydrology.png',dpi=300)


spotpy.analyser.plot_gelman_rubin(results, r_hat, fig_name='python_hydrology_convergence.png')


parameters = spotpy.parameter.get_parameters_array(spot_setup)

fig, ax = plt.subplots(nrows=4, ncols=2)
for par_id in range(len(parameters)):
    plot_parameter_trace(ax[par_id][0], results, parameters[par_id])
    plot_posterior_parameter_histogram(ax[par_id][1], results, parameters[par_id])

ax[-1][0].set_xlabel('Iterations')
ax[-1][1].set_xlabel('Parameter range')

plt.show()
fig.savefig('hydrology_parameters.png',dpi=300)