import numpy as np
import spotpy
from dream_setup import spotpy_setup
from scipy.optimize import fsolve
from functions import c

# Read the results from Assignment 4
results = spotpy.analyser.load_csv_results('HydrologyDREAM')
param_all = results[['parn', 'parD', 'parq', 'parLambda']]
param = param_all[-10:]

cThreshold = 2.5
cDectTrheshold = 10e-5


def getXfromConcentration(c, t, M, n, D, q ,Lambda):
    A = M / (np.sqrt(4 * np.pi * D * t))
    #x1 = (q*t)/n + sqrt(4*D*t * (np.log((c * sqrt(4*np.pi*D*t))/M) + Lambda * t))
    x1 = q * t / n + np.sqrt(-(np.log(c/A) + Lambda*t)/(4*D*t))
    x2 = q * t / n - np.sqrt(-(np.log(c/A) + Lambda*t)/(4*D*t))
    return [x2[0], x1[0]]

def getTfromConcentration(c,x,M,n,D,q,Lambda):
    # Estimate the location of t using information from assignment 4
    tEstimate = 64 * x / 100
    
    # Define the function and return the sollution
    func = lambda t: (M / (4 * 3.1415 * D * t)**(1/2))*(np.exp(-(x - q * t / n) ** 2 / (4 * D * t)))*(np.exp(-Lambda * t)) - c
    return fsolve(func, tEstimate)


for paramSet in param:
    
    # Get the time we reach the critical concentration at the well
    tCrit = getTfromConcentration(cThreshold, 100, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3])
    
    # Get the location where we reach the detection threshold 10 times before tCrit
    xDetection = getXfromConcentration(cDectTrheshold, tCrit - 10, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3])

    # Loop
        # Define sampling time points

        # Create synthetic samples

        # Run Dream with the new samples

        # Evalueate the output variance and check the improvement that came from this points
    

    print(tCrit)
    print(xDetection)
    




if __name__ == '__main__':

    #getXfromConcentration(0.00001, 40, 200, 0.3213, 0.6855, 0.4678, 0.01567)
    #getTfromConcentration(0.00001, 100, 200, 0.3213, 0.6855, 0.4678, 0.01567)
    pass