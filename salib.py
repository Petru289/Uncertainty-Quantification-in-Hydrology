from SALib.sample import saltelli
from SALib.analyze import sobol
import numpy as np
import parameters as par
from functions import c

problem = {
    'num_vars': 4,
    'names': ['n', 'D', 'q', 'Lambda'],
    'bounds': [[par.nRange[0], par.nRange[1]],
                [par.DRange[0], par.DRange[1]],
                [par.qRange[0], par.qRange[1]],
                [par.LambdaRange[0], par.LambdaRange[1]]]
}

for x in par.x:
    for t in par.t:
        param_values = saltelli.sample(problem, 1024)
        Y = np.zeros([param_values.shape[0]])

        for i, X in enumerate(param_values):
            Y[i] = c(x, t, par.M ,X[0],X[1],X[2],X[3])

        Si = sobol.analyze(problem, Y)

        print('----------------------------------------')
        print('x: ' + str(x) + ' t: ' + str(t))
        print('\n')
        print('First order sobol indices:')
        print(Si['S1'])
        print('\n')

        print('\n')
        print('Total sobol indices:')
        print(Si['ST'])
        print('\n')

        print('\n')
        print('Second order sobol indices: ')
        print("n-D:", Si['S2'][0,1])
        print("n-q:", Si['S2'][0,2])
        print("n-Lambda:", Si['S2'][0,3])
        print("D-q:", Si['S2'][1,2])
        print("D-Lambda:", Si['S2'][1,3])
        print("q-Lambda:", Si['S2'][2,3])
        print('\n')