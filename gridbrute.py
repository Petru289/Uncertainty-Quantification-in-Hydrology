import numpy as np
import pandas as pd
from scipy.optimize import optimize
from functions import c
import parameters as par
import obj


def gridbrute(method):

    #run 'pip install openpyxl' before executing the following line
    file = pd.read_excel(r'C:\Users\Gheorghe Pascu\OneDrive - tum.de\WiSe 21-22\Uncertainty Quantification in Hydrology\measurements.xlsx')
    obs = file.to_numpy()
    t = np.array([5, 10, 15, 20, 25, 30, 35])

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

    results = optimize.brute(f, ranges = bnds, full_output = 1)
        
    return [results[0], results[1]]



if __name__ == "__main__":
    gridbrute('least squares')
