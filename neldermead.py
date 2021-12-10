import numpy as np
import pandas as pd
import chaospy as cp
from scipy.optimize import minimize
from functions import c
import parameters as par
import obj

def neldermead(method):

    #run 'pip install openpyxl' before executing the following line
    file = pd.read_excel('./measurements.xlsx')
    obs = file.to_numpy()
    t = np.array([5, 10, 15, 20, 25, 30, 35])

    local_minimas = np.zeros((20, 4)) #Store the minimas here
    local_min_values = np.zeros(20) #Store the local minimal values here
    for i in range(20):
        #Nelder Mead algorithm can stumble upon local minimas that can be
        #far from being global minimas, so we randomly generate initial points
        #and pick the smallest value among the minimas obtained 
        n = cp.Uniform(par.nRange[0], par.nRange[1])
        D = cp.Uniform(par.DRange[0], par.DRange[1])
        q = cp.Uniform(par.qRange[0], par.qRange[1])
        Lambda = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
        jointdist = cp.J(n, D, q, Lambda)
        jointsample = jointdist.sample(size = 1)

        x0 = jointsample
        bnds = ((par.nRange[0], par.nRange[1]), (par.DRange[0], par.DRange[1]), (par.qRange[0], par.qRange[1]), (par.LambdaRange[0], par.LambdaRange[1]))

        if method == 'least squares':
            def f(param):
                return obj.ls(param[0], param[1], param[2], param[3])
        
        if method == 'least log-squares':
            def f(param):
                return obj.llogs(param[0], param[1], param[2], param[3])

        if method == 'relative error':
            def f(param):
                return obj.relerr(param[0], param[1], param[2], param[3])
        
        results = minimize(f, x0, method = 'Nelder-Mead', bounds = bnds)
        local_min_values[i] = results.fun
        local_minimas[i] = np.array(results.x)
    
    min = np.min(local_min_values)
    index = np.where(local_min_values == min)
    opt_par = np.reshape(local_minimas[index],(4,)) #Reshape, so that we get only 1 dimension
    
    return [opt_par, min]




if __name__ == "__main__":
    neldermead('least squares')

