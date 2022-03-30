import numpy as np
import spotpy
from dream_setup import spotpy_setup
from scipy.optimize import fsolve
from functions import c
import matplotlib.pyplot as plt
import scipy.optimize

# Read the results from Assignment 4
results = spotpy.analyser.load_csv_results('HydrologyDREAM')
param_all = results[['parn', 'parD', 'parq', 'parLambda']]
#param = param_all[-10:]
param = param_all[-1:]

cThreshold = 2.5
cDectTrheshold = 10e-5
sampleMethod = 'fixed'

def getX2(t, M, n, D, q ,Lambda):

    def cOfX(x):
        return c(x, t, M, n, D, q ,Lambda) - cDectTrheshold

    return scipy.optimize.broyden1(cOfX, [60], f_tol=1e-14)
    
def getXfromConcentration(c, t, M, n, D, q ,Lambda):
    A = M / (np.sqrt(4 * np.pi * D * t))
    #x1 = (q*t)/n + sqrt(4*D*t * (np.log((c * sqrt(4*np.pi*D*t))/M) + Lambda * t))
    x1 = q * t / n + np.sqrt(-(np.log(c/A) + Lambda*t)/(4*D*t))
    x2 = q * t / n - np.sqrt(-(np.log(c/A) + Lambda*t)/(4*D*t))
    return [x2[0], x1[0]]

        

def getT(c,x,M,n,D,q,Lambda, tEstimate):
    func = lambda t: (M / (4 * 3.1415 * D * t)**(1/2))*(np.exp(-(x - q * t / n) ** 2 / (4 * D * t)))*(np.exp(-Lambda * t)) - c
    return fsolve(func, tEstimate)

def getX(c, t, M, n, D, q ,Lambda, xEstimate):
    func = lambda x: (M / (4 * 3.1415 * D * t)**(1/2))*(np.exp(-(x - q * t / n) ** 2 / (4 * D * t)))*(np.exp(-Lambda * t)) - c
    return fsolve(func, xEstimate)

# Run dream one time

# Take the latest n parameter sets

# For each:
for paramSet in param:
    
    # Get the time we reach the critical concentration at the well
    tCrit = getT(cThreshold, 100, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], 64)
    
    tDetection = tCrit - 10

    # Get the location where we reach the detection threshold 10 times before tCrit
    # Solve for x c(x, t, M, n, D, q ,Lambda) and take the bigger sollution
    xDetection = getX2(tDetection, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3])
    
    sampleLocation = xDetection
    
    # Loop
        # Define sample time points (uniform, random...)

        # Create synthetic samples

        # Run Dream with the new samples

        # Evalueate the output variance and check the improvement that came from this points

        # Pick sollution with the smallest variance
    
    # Pick sollution with the smallest variance between all runs
    




if __name__ == '__main__':
    pass