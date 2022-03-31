import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from functions import c

# Load the parameters
import parameters as par

labelFontSize = 10

def histogramm():
    plt.rc('xtick', labelsize=8) 
    plt.rc('ytick', labelsize=8)

    # Initialize variables
    concentrations = np.zeros(shape = (len(par.x) * len(par.t), 1000))

    # Prepare the figure
    fig, axes = plt.subplots(nrows = len(par.x), ncols = len(par.t),figsize = (8,6))
    plt.subplots_adjust(left=0.1, bottom = 0.1, right = 0.9, top = 0.9, wspace = 0.7, hspace = 0.8)
    fig.suptitle('Histogramm of the $\it{concentration}$ for different $\it{x}$ and $\it{t}$', fontsize=16)

    for distIndex, distance in enumerate(par.x): # Loop over all distances
        for timeIndex, timestep in enumerate(par.t): # For all timesteps

            # This is just an index to store all concentration values in a matrix
            concentrationIndex = 3 * distIndex + timeIndex
            
            # Sample parameters from a joint distribution
            n = cp.Uniform(par.nRange[0], par.nRange[1])
            D = cp.Uniform(par.DRange[0], par.DRange[1])
            q = cp.Uniform(par.qRange[0], par.qRange[1])
            Lambda = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
            jointdist = cp.J(n, D, q, Lambda)
            jointsample = jointdist.sample(size = 1000)

            # Calculate concentration
            concentrations[concentrationIndex, :] = c(distance, timestep, par.M, jointsample[0,:], jointsample[1,:], jointsample[2,:], jointsample[3,:])
            
            # Statistical measures
            mean = np.round(np.mean(concentrations[concentrationIndex, :]), 2)
            std =  np.round(np.std(concentrations[concentrationIndex, :]), 2)
            
            # Plot subplots
            axes[distIndex, timeIndex].hist(concentrations[concentrationIndex,:], 20) #20 bars optimal, but with sample size = 1000
            axes[distIndex, timeIndex].set_title('x = {} m, t = {} days'.format(distance, timestep), fontsize = labelFontSize)
            #axes[distIndex, timeIndex].set_xlabel('concentration [$kg/m^3$]', fontsize = labelFontSize)
            #axes[distIndex, timeIndex].set_ylabel('Frequency', fontsize = labelFontSize)
            textstring = '\n'.join((
                '$\mu=%.2f$' % float(mean),
                '$\sigma=%.2f$' % std
            ))
            axes[distIndex, timeIndex].text(0.75, 0.95, textstring, transform=axes[distIndex, timeIndex].transAxes, fontsize=8,
            verticalalignment='top')

    axes[2, 0].set_xlabel('concentration [$kg/m^3$]', fontsize = labelFontSize)
    axes[2, 1].set_xlabel('concentration [$kg/m^3$]', fontsize = labelFontSize)
    axes[2, 2].set_xlabel('concentration [$kg/m^3$]', fontsize = labelFontSize)

    axes[0,0].set_ylabel('Frequency', fontsize = labelFontSize)
    axes[1,0].set_ylabel('Frequency', fontsize = labelFontSize)
    axes[2,0].set_ylabel('Frequency', fontsize = labelFontSize)
    plt.show()

if __name__ == "__main__":
    histogramm()