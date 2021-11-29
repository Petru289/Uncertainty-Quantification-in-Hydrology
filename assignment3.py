import pandas as pd
from grid import grid
from montecarlo import montecarlo
from neldermead import neldermead

#For the loss function you want to use, choose between 'least squares', 'least log-squares', 'relative error'
#Additionaly, as first argument, 
#grid() expects the number of discretization points on each axis
#montecarlo() expects the sample size

grid(30, 'least squares') 

montecarlo(30, 'least squares')

neldermead('least squares')

