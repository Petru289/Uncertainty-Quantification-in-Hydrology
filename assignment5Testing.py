from turtle import color
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


def getX2(t, M, n, D, q ,Lambda, xEstimate, cThreshold):

    def cOfX(x):
        return c(x, t, M, n, D, q ,Lambda) - cThreshold

    return scipy.optimize.broyden1(cOfX, [xEstimate], f_tol=1e-10)
    

    

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

    fig, axes = plt.subplots(2,2)
    timeSteps = np.linspace(1,100,100)
    axes[0,0].plot(timeSteps, c(100, timeSteps, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]), 'c')
    axes[0,0].axvline(x=tCrit, color='r')
    axes[0,0].set_xticks([0,25,50,75,100, int(tCrit)])
    axes[0,0].set_xlabel('Days')
    axes[0,0].set_ylabel('C')
    axes[0,0].set_title('Concentration over time at x = 100m')
    

    xSteps = np.linspace(80,120,100)
    axes[0,1].plot(xSteps, c(xSteps, tCrit, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]), 'c')
    axes[0,1].axvline(x=100, color='r')
    axes[0,1].axhline(y=2.5, color='r')
    axes[0,1].set_xlabel('[m]')
    axes[0,1].set_ylabel('C')
    axes[0,1].set_title('Concentration over x at time = tCrit')


    x = np.linspace(40,150,100)
    axes[1,0].plot(x, c(x, tCrit - 10, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]), 'm')
    axes[1,0].plot(x, c(x, tCrit, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]), 'c')
    axes[1,0].set_xlabel('[m]')
    axes[1,0].set_ylabel('C')
    axes[1,0].set_title('Concentration over x at time = [tCrit - 10, tCrit]')
    axes[1,0].legend(['t = tCrit - 10', 't = tCrit'])
    
    tDetectionStart = tCrit - 10
    
    xDetection = [0,0]
    xDetection[0] = getX2(tDetectionStart, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], 120, cDectTrheshold)
    xDetection[1] = getX2(tDetectionStart, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], 40, cDectTrheshold)
    xThresholdLower = getX2(tDetectionStart, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], 70, cThreshold)
    xThresholdUpper = getX2(tDetectionStart, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], 85, cThreshold)
    print([xThresholdLower, xThresholdUpper])
    axes[1,0].axvline(x=xDetection[0], color='k')
    axes[1,0].axvline(x=xDetection[1], color='k')
    axes[1,0].axvline(x=xThresholdLower, color='r')

    timeSteps = np.linspace(tCrit - 30, tCrit, 100)
    axes[1,1].plot(timeSteps, c(xThresholdLower, timeSteps, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]))
    axes[1,1].axvline(x=tCrit - 10, color='r')
    axes[1,1].text((tCrit - 10) + 0.5, 0.1, 'tCrit - 10', rotation=90)
    axes[1,1].set_xlabel('Days')
    axes[1,1].set_ylabel('C')

    timeFirstSample = getTfromConcentration(cThreshold, xThresholdLower, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3], 43)
    sampleTimePoints = np.linspace(timeFirstSample, tCrit - 10, 5)
    axes[1,1].plot(sampleTimePoints, c(xThresholdLower, sampleTimePoints, 200, paramSet[0],paramSet[1],paramSet[2],paramSet[3]), 'ko')

    fig.tight_layout()
    plt.show()
    




if __name__ == '__main__':

    #getXfromConcentration(0.00001, 40, 200, 0.3213, 0.6855, 0.4678, 0.01567)
    #getTfromConcentration(0.00001, 100, 200, 0.3213, 0.6855, 0.4678, 0.01567)
    pass