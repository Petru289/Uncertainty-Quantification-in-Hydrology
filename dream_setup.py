import spotpy 
import numpy as np
import parameters as par
from functions import c
import obj
from spotpy.likelihoods import gaussianLikelihoodMeasErrorOut as GausianLike

#See documentation: https://spotpy.readthedocs.io/en/latest/Rosenbrock/

#Inspired from spotpy.spot_setupdtlz1.py. The issue is that the model (our concentration function) changes with the observations.
#The following code is inspired from spot_setupdtlz1.py. Not sure if it's the correct approach; paramaters don't seem to converge to a distribution

n_obs = 14 
time = np.array([5, 10, 15, 20, 25, 30, 35])
space = np.array([5, 50])

# def conc(x):
#     f = []
#     for i in range(0, n_obs):
#         if i <= 6:
#             f.append(c(5, t[i], 200, x[0], x[1], x[2], x[3]))
#         else:
#             f.append(c(50, t[i - 7], 200, x[0], x[1], x[2], x[3]))
#     return f


class spotpy_setup(object):
    def __init__(self):
        self.params = [spotpy.parameter.Uniform('n',par.nRange[0],par.nRange[1],1.5,3.0,par.nRange[0],par.nRange[1]),
                       spotpy.parameter.Uniform('D',par.DRange[0],par.DRange[1],1.5,3.0,par.DRange[0],par.DRange[1]),
                       spotpy.parameter.Uniform('q',par.qRange[0],par.qRange[1],1.5,3.0,par.qRange[0],par.qRange[1]),
                       spotpy.parameter.Uniform('Lambda',par.LambdaRange[0],par.LambdaRange[1],1.5,3.0,par.LambdaRange[0],par.LambdaRange[1])
                       ]
    def parameters(self):
        return spotpy.parameter.generate(self.params)


    # def simulation(self,vector):      
    #         x = np.array(vector)
    #         sim = conc(x)
    #         return sim


    # def evaluation(self):
    #     observations = obj.obs[:,0].tolist() + obj.obs[:,1].tolist()
    #     return observations


    # def objectivefunction(self, simulation, evaluation):    
    #     objective = []
    #     for i,f in enumerate(simulation):
    #         objective.append(GausianLike(data = [evaluation[i]], comparedata = [simulation[i]]))
    #     return objective
    
    def simulation(self, x):
        n, D, q, Lambda = x
        sim = []
        for s in range(len(space)):
            for t in range(len(time)):
                sim.append(c(space[s], time[t], 200, n, D, q, Lambda))
        return sim

    def evaluation(self):
        observations = obj.obs[:,0].tolist() + obj.obs[:,1].tolist()
        return observations

    def objectivefunction(self, simulation, evaluation):    
        objective = GausianLike(data = evaluation, comparedata = simulation)
        return objective

