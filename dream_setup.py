import spotpy 
import numpy as np
import parameters as par
from functions import c
import obj

#See documentation: https://spotpy.readthedocs.io/en/latest/Rosenbrock/

#Inspired from spotpy.spot_setupdtlz1.py

n_obs = 14 
t = np.array([5, 10, 15, 20, 25, 30, 35])

def conc(x):
    f = []
    for i in range(0, n_obs):
        if i <= 6:
            f.append(c(5, t[i], 200, x[0], x[1], x[2], x[3]))
        else:
            f.append(c(50, t[i - 7], 200, x[0], x[1], x[2], x[3]))
    return f


class spotpy_setup(object):
    def __init__(self):
        self.params = [spotpy.parameter.Uniform('n',par.nRange[0],par.nRange[1],1.5,3.0,par.nRange[0],par.nRange[1]),
                       spotpy.parameter.Uniform('D',par.DRange[0],par.DRange[1],1.5,3.0,par.DRange[0],par.DRange[1]),
                       spotpy.parameter.Uniform('q',par.qRange[0],par.qRange[1],1.5,3.0,par.qRange[0],par.qRange[1]),
                       spotpy.parameter.Uniform('Lambda',par.LambdaRange[0],par.LambdaRange[1],1.5,3.0,par.LambdaRange[0],par.LambdaRange[1])
                       ]
    def parameters(self):
        return spotpy.parameter.generate(self.params)


    def simulation(self,vector):      
            x = np.array(vector)
            sim = conc(x)
            return sim


    def evaluation(self):
        observations = obj.obs[:,0].tolist() + obj.obs[:,1].tolist()
        return observations


    def objectivefunction(self, simulation, evaluation):
        objective = []
        for i,f in enumerate(simulation):
            objective.append(spotpy.objectivefunctions.mae(evaluation = [evaluation[i]], simulation = [simulation[i]]))
        return objective
    