import numpy as np
import matplotlib.pyplot as plt
import chaospy as cp
from functions import c, dcdn, dcdq, dcdD, dcdLambda

# Load the input parameters
import parameters as par

def sensitivityAnalysis(parameter):
    '''Performs a one-at-a-time sensitivity analysis changing only the parameter given. Plots the result.'''

    # Define constant parameters
    n = np.mean(par.nRange)
    D = np.mean(par.DRange)
    q = np.mean(par.qRange)
    Lambda = np.mean(par.LambdaRange)
    #D - max, to illustrate the interaction between parameters - see report
    # n = np.mean(par.nRange)
    # D = np.max(par.DRange)
    # q = np.mean(par.qRange)
    # Lambda = np.mean(par.LambdaRange)
    x = par.x
    t = par.t
    M = par.M

    # Overwrite the parameter we want to change
    if parameter == 'n':
        unif = cp.Uniform(par.nRange[0], par.nRange[1]) 
        n = xstep = unif.sample(size = 1000) #sample size 1000 for a better picture of the distribution
        xLabel = 'porosity [-]'
        derivative = dcdn

    elif parameter == 'q':
        unif = cp.Uniform(par.qRange[0], par.qRange[1]) 
        q = xstep = unif.sample(size = 1000)
        xLabel = 'specific discharge [$m/d$]'
        derivative = dcdq 

    elif parameter == 'D':
        unif = cp.Uniform(par.DRange[0], par.DRange[1]) 
        D = xstep = unif.sample(size = 1000)
        xLabel = 'dispersion coefficient [$m^2/d$]'
        derivative = dcdD 

    elif parameter == 'Lambda':
        unif = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1]) 
        Lambda = xstep = unif.sample(size = 1000)
        xLabel = 'decay rate [$1/d$]'
        derivative = dcdLambda 
    
    concentrations = np.zeros(shape = (len(x) * len(t), 1000))
    dconcentrations = np.zeros(shape = (len(x) * len(t), 1000))

    # Prepare the figures
    fig1, axes1 = plt.subplots(nrows = len(x), ncols = len(t))
    fig1.tight_layout()
    fig2, axes2 = plt.subplots(nrows = len(x), ncols = len(t))
    fig2.tight_layout()
    fig3, axes3 = plt.subplots(nrows = len(x), ncols = len(t))
    fig3.tight_layout()
    axisLabelFontSize = 10

    for distIndex, distance in enumerate(x):
        for timeIndex, timestep in enumerate(t):

            rownumber = 3 * distIndex + timeIndex
            
            concentrations[rownumber, :] = c(distance, timestep, M, n, D, q, Lambda)
            dconcentrations[rownumber, :] = derivative(distance, timestep, M, n, D, q, Lambda)
            mean = np.round(np.mean(concentrations[rownumber, :]), 2)
            std =  np.round(np.std(concentrations[rownumber, :]), 2)

            # Clean overflowing skalars
            # This is needed because for i = 2 and j = 0 we encountered an overflow exeption during the np.hist function call
            concentrations[rownumber,:] = [0 if x<1e-300 else x for x in concentrations[rownumber,:]]
            
            # Plot histogramm
            axes1[distIndex, timeIndex].hist(concentrations[rownumber,:], 20) #20 -> more bars (but with samplesize = 1000), better picture of distribution
            axes1[distIndex, timeIndex].set_title('x = {}, t = {}'.format(distance, timestep), fontsize = axisLabelFontSize)
            axes1[distIndex, timeIndex].set_xlabel('concentration [$kg/m^3$], mean = {}, sd = {}'.format(mean, std), fontsize = axisLabelFontSize)
            print('mean = ', mean, 'sd = ', std)
            # Plot concentration
            axes2[distIndex, timeIndex].scatter(xstep, concentrations[rownumber,:], s = 5, c = "firebrick")
            axes2[distIndex, timeIndex].set_title('x = {}, t = {}'.format(distance, timestep), fontsize = axisLabelFontSize)
            axes2[distIndex, timeIndex].set_ylabel('concentration [$kg/m^3$]', fontsize = axisLabelFontSize)
            axes2[distIndex, timeIndex].set_xlabel(xLabel, fontsize = axisLabelFontSize)

            # Plot derivative
            axes3[distIndex, timeIndex].scatter(xstep, dconcentrations[rownumber,:], s = 5, c = "darkorange")
            axes3[distIndex, timeIndex].set_title('x = {}, t = {}'.format(distance, timestep), fontsize = axisLabelFontSize)
            axes3[distIndex, timeIndex].set_ylabel('change in concentration', fontsize = axisLabelFontSize)
            axes3[distIndex, timeIndex].set_xlabel(xLabel, fontsize = axisLabelFontSize)

    plt.show()

if __name__ == "__main__":
    sensitivityAnalysis('D')