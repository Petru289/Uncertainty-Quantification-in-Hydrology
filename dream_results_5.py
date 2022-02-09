import numpy as np
import spotpy
from dream_setup_5 import spotpy_setup
from dream_run_5 import dream_run_5
import matplotlib.pyplot as plt
from spotpy.analyser import plot_parameter_trace
from spotpy.analyser import plot_posterior_parameter_histogram
from functions import c


def dream_results_5(new_meas, plot = False):

    spot_setup = spotpy_setup(new_meas)

    time = [5, 10, 15, 20, 25, 30, 35]

    r_hat = dream_run_5(new_meas)
    results = spotpy.analyser.load_csv_results('HydrologyDREAM_5')


    param_all = results[['parn', 'parD', 'parq', 'parLambda']]
    param = param_all[-100:] # Only the last 100 iterations of the chain
    param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))
    
    if plot == True:

        new_well = new_meas[0]
        new_times = new_meas[1]
        #new_obs = new_meas[2]

        sim1 = []
        sim2 = []
        sim3 = []
        for t in time:
            sim1.append(c(5, t, 200, param[:,0], param[:,1], param[:,2], param[:,3]))
            sim2.append(c(50, t, 200, param[:,0], param[:,1], param[:,2], param[:,3]))
        for t in new_times:
            sim3.append(c(new_well, t, 200, param[:,0], param[:,1], param[:,2], param[:,3]))

        obsindices1 = np.arange(0,7,1)
        obsindices2 = np.arange(7,14,1)
        obsindices3 = np.arange(0,5,1)
        fig, axes = plt.subplots(nrows = 3, ncols = 1)
        for i in range(3):
            if i == 0:
                sim = sim1
            elif i == 1:
                sim = sim2
            else: 
                sim = sim3
            q5, q25, q75, q95 = [], [], [], [] # We don't use q25 and q75
            for j in range(len(sim)):
                q5.append(np.percentile(sim[j], 2.5))
                q95.append(np.percentile(sim[j], 97.5))
            axes[i].plot(q5,color='dimgrey',linestyle='solid')
            axes[i].plot(q95,color='dimgrey',linestyle='solid')
            axes[i].fill_between(np.arange(0,len(q5),1),q5,q95,facecolor='dimgrey',zorder=0,
                            linewidth=0,label='parameter uncertainty')
            
            if i == 0:
                axes[i].plot(obsindices1, spot_setup.evaluation()[:len(obsindices1)], 'r.', label='observations')
                axes[i].set_ylim(-0.3,40)
                axes[i].set_ylabel('C(5, . ) [kg/m^3]')

            elif i == 1:
                axes[i].plot(obsindices2 - 7, spot_setup.evaluation()[len(obsindices1):len(obsindices1) + len(obsindices2)], 'r.', label='observations')
                axes[i].set_ylim(-0.3,15)
                axes[i].set_ylabel('C(50, . ) [kg/m^3]')

            else:
                axes[i].plot(obsindices3, spot_setup.evaluation()[len(obsindices1) + len(obsindices2):len(obsindices1) + len(obsindices2) + len(obsindices3)], 'r.', label='observations')
                axes[i].set_ylim(-0.3,15)
                axes[i].set_ylabel('C(, . ) [kg/m^3]')
            
            axes[i].set_xlabel('time [days]') 
            if i == 0 or i == 1:
                axes[i].set_xticklabels(labels = ['5'] + time) # I know this is weird, but otherwise it doesn't put the first tick, Idk why
            else:
                axes[i].set_xticklabels(labels = new_times) #Labels to fix here
            axes[i].set_xlim(-0.05,6.05)   
            axes[i].legend()
            fig.savefig('python_hydrology_5.png', dpi = 300)                   
                                        

        spotpy.analyser.plot_gelman_rubin(results, r_hat, fig_name='python_hydrology_convergence_5.png')
        parameters = spotpy.parameter.get_parameters_array(spot_setup)

        fig, ax = plt.subplots(nrows=4, ncols=2)
        for par_id in range(len(parameters)):
            plot_parameter_trace(ax[par_id][0], results, parameters[par_id])
            plot_posterior_parameter_histogram(ax[par_id][1], results, parameters[par_id])

        ax[-1][0].set_xlabel('Iterations')
        ax[-1][1].set_xlabel('Parameter range')

        plt.show()
        fig.savefig('hydrology_parameters_5.png',dpi=300)

    return param

