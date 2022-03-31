from tkinter import E
import spotpy 
import sys
import numpy as np
sys.path.insert(0, './Assignment1')
sys.path.insert(0, './Assignment3')
import parameters as par
from functions import c
import obj
from spotpy.likelihoods import gaussianLikelihoodMeasErrorOut as GausianLike

#See documentation: https://spotpy.readthedocs.io/en/latest/Rosenbrock/


n_obs = 14 
time = np.array([5, 10, 15, 20, 25, 30, 35])
time2 = np.array([35, 38, 41, 44, 47, 50, 53])
space = np.array([5, 50])

class spotpy_setup(object):
    def __init__(self, obsindices):
        self.params = [spotpy.parameter.Uniform('n',par.nRange[0],par.nRange[1],1.5,3.0,par.nRange[0],par.nRange[1]),
                       spotpy.parameter.Uniform('D',par.DRange[0],par.DRange[1],1.5,3.0,par.DRange[0],par.DRange[1]),
                       spotpy.parameter.Uniform('q',par.qRange[0],par.qRange[1],1.5,3.0,par.qRange[0],par.qRange[1]),
                       spotpy.parameter.Uniform('Lambda',par.LambdaRange[0],par.LambdaRange[1],1.5,3.0,par.LambdaRange[0],par.LambdaRange[1])
                       ]
        self.obsindices1 = obsindices[obsindices <= 6]
        self.obsindices2 = obsindices[obsindices > 6]
    def parameters(self):
        return spotpy.parameter.generate(self.params)
    
    def simulation(self, x):
        n, D, q, Lambda = x
        sim = []
        for i in self.obsindices1:
            sim.append(c(5, time[i], 200, n, D, q, Lambda))
        for i in self.obsindices2:
            sim.append(c(50, time[i-7], 200, n, D, q, Lambda)) #For experiment in ass. 5, change time with time2 here
        return sim

    def evaluation(self):
        observations = []
        for i in self.obsindices1:
            observations.append(obj.obs[i,0])
        for i in self.obsindices2:
            observations.append(obj.obs[i-7,1]) #For experiment in assignment 5, change here 1 with 2
        return observations

    def objectivefunction(self, simulation, evaluation):    
        objective = GausianLike(data = evaluation, comparedata = simulation)
        return objective

