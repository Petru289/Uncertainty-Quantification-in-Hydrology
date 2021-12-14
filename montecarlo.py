from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import chaospy as cp
import parameters as par
from functions import c
import matplotlib.pyplot as plt
import obj


def montecarlo(N, method, plot = False):  #N = sample size
 
    t = np.array([5, 10, 15, 20, 25, 30, 35])

    #run 'pip install openpyxl' before executing the following line
    file = pd.read_excel('./measurements.xlsx')
    obs = file.to_numpy()

    # Sample
    n = cp.Uniform(par.nRange[0], par.nRange[1])
    D = cp.Uniform(par.DRange[0], par.DRange[1])
    q = cp.Uniform(par.qRange[0], par.qRange[1])
    Lambda = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
    jointdist = cp.J(n, D, q, Lambda)
    jointsample = jointdist.sample(size = N)
    x = np.sort(jointsample[0,:])
    y = np.sort(jointsample[1,:])
    z = np.sort(jointsample[2,:])
    w = np.sort(jointsample[3,:])
    
    #This method is maybe inefficient because of sorting the arrays
    X, Y, Z, W = np.meshgrid(x, y, z, w)
    error = 0

    # Calculate the error depending on the method
    error = obj.getError(X, Y, Z, W, method)
    min = np.min(error)
    indices = np.where(error == min)

    #VERY CAREFUL HERE, MESHGRID TRICKS YOU. INDICES ARE CORRECT NOW.
    i1, i2, i3, i4 = int(indices[1]), int(indices[0]), int(indices[2]), int(indices[3])
    n_opt = x[i1]
    D_opt = y[i2]
    q_opt = z[i3]
    Lambda_opt = w[i4]

    #Plotting projections of the objective functions
    if plot:

        # Get X and Y values. If you want to plot different parameters change them here.
        # And don't forget to adjust the labels for the plot as well.
        xValues = X[:, :, i3, i4]
        yValues = Y[:, :, i3, i4]
        zValues = Z[:, :, i3, i4]
        wValues = W[:, :, i3, i4]
        
        # Evalueate the objective function
        f = obj.getError(xValues, yValues, zValues, wValues, method)

        # Define meaningful levels. This highly depends on the objective function.
        if method == 'least squares':
            vmin = 0.02
            levels = np.linspace(vmin, f.max()**(1/6), 10) ** 6
        
        if method == 'least log-squares':
            vmin = 0.02
            levels = np.linspace(vmin, f.max()**(1/6), 10) ** 6
        
        if method == 'relative error':
            vmin = f.min()
            levels = np.linspace(vmin, f.max(), 10)
            levels = np.linspace(vmin, f.max()**(1/10), 10) ** 10

        #levels = [0.02, 0.3, 4.0, 14.0, 38.0, 60.0, 200, 400, 800, 2000, 4000]

        cpl = plt.contour(xValues, yValues, f, levels, colors='black', linestyles='dashed', linewidths=1)
        plt.clabel(cpl, inline=1, fontsize=10)
        cpl = plt.contourf(xValues, yValues, f, levels, cmap='YlGnBu_r', vmin=vmin/2, vmax=f.max(), norm=LogNorm())
        plt.xlabel('Porosity $[-]$')
        plt.ylabel('Dispersion Coefficient $[m^2/d]$')
        plt.title('Value of the objective function')
        plt.show()

    return [np.array([n_opt, D_opt, q_opt, Lambda_opt]), min]



if __name__ == "__main__":
    opt_par = montecarlo(40, 'least log-squares', plot=True)
