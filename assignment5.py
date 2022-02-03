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



def getX(c, t, M, n, D, q ,Lambda, xEstimate):
    func = lambda x: (M / (4 * 3.1415 * D * t)**(1/2))*(np.exp(-(x - q * t / n) ** 2 / (4 * D * t)))*(np.exp(-Lambda * t)) - c
    return fsolve(func, xEstimate)

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

        

def getTfromConcentration(c,x,M,n,D,q,Lambda, tEstimate):
    func = lambda t: (M / (4 * 3.1415 * D * t)**(1/2))*(np.exp(-(x - q * t / n) ** 2 / (4 * D * t)))*(np.exp(-Lambda * t)) - c
    return fsolve(func, tEstimate)


for paramSet in param:
    
    # Get the time we reach the critical concentration at the well
    tCrit = getTfromConcentration(cThreshold, 100, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], 64)

    x = np.linspace(40,150,100)
    plt.plot(x, c(x, tCrit - 10, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]), 'c')
    plt.plot(x, c(x, tCrit, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]), 'm')
    
    tDetectionStart = tCrit - 10

    # Get the location where we reach the detection threshold 10 times before tCrit
    #xDetection = getX(cDectTrheshold, tCrit - 10, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], 50)
    
    xDetection = getX2(tDetectionStart, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3])
    print(xDetection)

    #tDetectionEnd = getTfromConcentration(cDectTrheshold, xDetection, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], tDetectionStart + 0.5)

    #print([tDetectionStart[0], tDetectionEnd[0]])

    # Test:
    #print(c(xDetection, np.array([tDetectionStart, tDetectionEnd]), 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]))

    if sampleMethod == 'fixed':
        samplePoints = np.array([tDetectionStart, tDetectionStart + 1, tDetectionStart + 2, tDetectionStart + 3, tDetectionStart + 4])

        # Create Samples
        #samples = c(xDetection, samplePoints, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3])
        #print(samples)
    
    
    
    
    
    # Loop
        # Define sampling time points

        # Create synthetic samples

        # Run Dream with the new samples

        # Evalueate the output variance and check the improvement that came from this points
plt.show()    
    




if __name__ == '__main__':

    #getXfromConcentration(0.00001, 40, 200, 0.3213, 0.6855, 0.4678, 0.01567)
    #getTfromConcentration(0.00001, 100, 200, 0.3213, 0.6855, 0.4678, 0.01567)
    pass