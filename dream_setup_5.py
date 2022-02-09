from tkinter import E
import spotpy 
import numpy as np
import parameters as par
from functions import c
import obj
from spotpy.likelihoods import gaussianLikelihoodMeasErrorOut as GausianLike

#See documentation: https://spotpy.readthedocs.io/en/latest/Rosenbrock/


n_obs = 14 
time = np.array([5, 10, 15, 20, 25, 30, 35])
space = np.array([5, 50])


class spotpy_setup(object):
    def __init__(self, new_meas):
        self.params = [spotpy.parameter.Uniform('n',par.nRange[0],par.nRange[1],1.5,3.0,par.nRange[0],par.nRange[1]),
                       spotpy.parameter.Uniform('D',par.DRange[0],par.DRange[1],1.5,3.0,par.DRange[0],par.DRange[1]),
                       spotpy.parameter.Uniform('q',par.qRange[0],par.qRange[1],1.5,3.0,par.qRange[0],par.qRange[1]),
                       spotpy.parameter.Uniform('Lambda',par.LambdaRange[0],par.LambdaRange[1],1.5,3.0,par.LambdaRange[0],par.LambdaRange[1])
                       ]
        self.new_well = new_meas[0]
        self.new_times = new_meas[1]
        self.new_obs = new_meas[2]
    def parameters(self):
        return spotpy.parameter.generate(self.params)
    
    def simulation(self, x):
        n, D, q, Lambda = x
        sim = []
        for i in range(len(time)):
            sim.append(c(5, time[i], 200, n, D, q, Lambda))
        for i in range(len(time)):
            sim.append(c(50, time[i], 200, n, D, q, Lambda))
        for i in range(len(self.new_times)):
            sim.append(c(self.new_well, self.new_times[i], 200, n, D, q, Lambda))
        return sim

    def evaluation(self):
        observations = []
        for i in range(obj.obs.shape[0]):
            observations.append(obj.obs[i,0])
        for i in range(obj.obs.shape[0]):
            observations.append(obj.obs[i,1])
        for i in range(len(self.new_times)):
            observations.append(self.new_obs[i])
        return observations

    def objectivefunction(self, simulation, evaluation):    
        objective = GausianLike(data = evaluation, comparedata = simulation)
        return objective

