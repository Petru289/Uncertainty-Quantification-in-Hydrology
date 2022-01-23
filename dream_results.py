import numpy as np
import spotpy
from dream_setup import spotpy_setup
from dream_run import dream_run
import matplotlib.pyplot as plt
from spotpy.analyser import plot_parameter_trace
from spotpy.analyser import plot_posterior_parameter_histogram


#obsindices = np.random.choice(13, size = 7, replace = False)
#obsindices = np.sort(obsindices)
#obsindices = np.array([0, 2, 5, 7, 10, 12, 13])
#obsindices1 = np.array([0, 2, 5])
#obsindices2 = np.array([7, 10, 12, 13])

#For all the observations, choose obsindices = np.arange(0,14,1)

obsindices = np.arange(0,14,1)
obsindices1 = obsindices[obsindices <= 6]
obsindices2 = obsindices[obsindices > 6]

ind1 = np.full(7, -1.0)
ind2 = np.full(7, -1.0)
for i in range(7):
    if i in obsindices1:
        ind1[i] = i
    if i in obsindices2 - 7:
        ind2[i] = i
if -1 in ind1:
    ind1[ind1 == -1] = np.nan
if -1 in ind2:
    ind2[ind2 == -1] = np.nan

spot_setup = spotpy_setup(obsindices)

#M, d, k = 4, 2, 4 #https://spotpy.readthedocs.io/en/latest/Sensitivity_analysis_with_FAST/
#rep = (1 + 4 * M ** 2 * (1 + (k - 2) *f d)) * k
#sampler = spotpy.algorithms.fast(spot_setup,  dbname='FAST_hydrology',  dbformat='csv')
#sampler.sample(repetitions = rep)
#results_FAST = spotpy.analyser.load_csv_results('FAST_hydrology')
#spotpy.analyser.plot_fast_sensitivity(results_FAST, number_of_sensitiv_pars=4)

time = [5, 10, 15, 20, 25, 30, 35]

r_hat = dream_run(obsindices)
results = spotpy.analyser.load_csv_results('HydrologyDREAM')
fields = [word for word in results.dtype.names if word.startswith('sim')]
Fields1 = fields[:len(obsindices1)]
Fields2 = fields[len(obsindices1):]

fig, axes = plt.subplots(nrows = 2, ncols = 1)

for i in range(2):
    Fields = Fields1 if i == 0 else Fields2
    q5, q25, q75, q95 = [], [], [], [] 
    for field in Fields:
        q5.append(np.percentile(results[field][-100:-1], 2.5))
        q95.append(np.percentile(results[field][-100:-1], 97.5))
    axes[i].plot(q5,color='dimgrey',linestyle='solid')
    axes[i].plot(q95,color='dimgrey',linestyle='solid')
    q5plot = []
    q95plot = []
    k = 0
    ind = ind1 if i==0 else ind2
    for j in range(7):
        if np.isnan(ind[j]):
            q5plot.append(0)
            q95plot.append(0.5)
        else:
            q5plot.append(q5[k])
            q95plot.append(q95[k])
            k += 1
    axes[i].fill_between(np.arange(0,7,1),q5plot,q95plot,facecolor='dimgrey',zorder=0,
                    linewidth=0,label='parameter uncertainty')
    
    if i == 0:
        axes[i].plot(obsindices1, spot_setup.evaluation()[:len(obsindices1)], 'r.', label='observations')
        axes[i].set_ylim(-0.3,40)
        axes[i].set_ylabel('C(5, . ) [kg/m^3]')

    else:
        axes[i].plot(obsindices2 - 7, spot_setup.evaluation()[len(obsindices1):], 'r.', label='observations')
        axes[i].set_ylim(-0.3,15)
        axes[i].set_ylabel('C(50, . ) [kg/m^3]')
    
    axes[i].set_xlabel('time [days]') 
    axes[i].set_xticklabels(labels = time) #I don t know why it misses the first tick. To be fixed.
    axes[i].set_xlim(-0.05,6.05)   
    axes[i].legend()
    fig.savefig('python_hydrology.png', dpi = 300)               
                   

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


# SCRAP PAPER

# def plotter(x):
#     if x == 5:
#         Fields = fields[:7]
#     else:
#         Fields = fields[7:]
#     fig = plt.figure(figsize=(16,9))
#     ax = plt.subplot(1,1,1)
#     q5, q25, q75, q95 = [], [], [], [] 
#     for field in Fields:
#         q5.append(np.percentile(results[field][-100:-1], 2.5))
#         q95.append(np.percentile(results[field][-100:-1], 97.5))
#     ax.plot(q5,color='dimgrey',linestyle='solid')
#     ax.plot(q95,color='dimgrey',linestyle='solid')
#     ax.fill_between(np.arange(0,len(q5),1),list(q5),list(q95),facecolor='dimgrey',zorder=0,
#                     linewidth=0,label='parameter uncertainty')
#     if x == 5:
#         ax.plot(spot_setup.evaluation()[:7],'r.', label='data')
#         ax.set_ylim(-0.1,40)
#         ax.set_xlim(0,7)
#         ax.legend()
#         fig.savefig('python_hydrology1.png', dpi = 300)
#     else:
#         ax.plot(spot_setup.evaluation()[7:],'r.', label='data')
#         ax.set_ylim(-0.1,15)
#         ax.set_xlim(0,6)
#         ax.legend()
#         fig.savefig('python_hydrology2.png', dpi = 300)

# plotter(x = 5)
# plotter(x = 50)


# fig = plt.figure(figsize=(16,9))
# ax = plt.subplot(1,1,1)
# q5,q25,q75,q95=[],[],[],[]
# for field in fields:
#     q5.append(np.percentile(results[field][-100:-1],2.5))
#     q95.append(np.percentile(results[field][-100:-1],97.5))
# ax.plot(q5,color='dimgrey',linestyle='solid')
# ax.plot(q95,color='dimgrey',linestyle='solid')
# ax.fill_between(np.arange(0,len(q5),1),list(q5),list(q95),facecolor='dimgrey',zorder=0,
#                 linewidth=0,label='parameter uncertainty')  
# ax.plot(spot_setup.evaluation(),'r.',label='data')
# ax.set_ylim(-5,70)
# ax.set_xlim(0,15)
# ax.legend()
# fig.savefig('python_hydrology.png',dpi=300)