from tabnanny import check
import numpy as np
import spotpy
from scipy.optimize import fsolve
from functions import c
from dream_results_5 import dream_results_5
from assignment5_P import choose
import matplotlib.pyplot as plt



# Uncomment these 3 lines if you want to run the process one more time. It takes about 90 minutes on my PC. 
# Otherwise, just read the output from the files.
# choose('close')
# choose('halfway')
# choose('far')

# Close, halfway and far relative to the position of the drinking well

stds_close = np.loadtxt("stds_close.txt", dtype='f', delimiter=' ') #pos here is about 88 m
stds_halfway = np.loadtxt("stds_halfway.txt", dtype='f', delimiter=' ') #pos here is about 67 m
stds_far = np.loadtxt("stds_far.txt", dtype='f', delimiter=' ') #pos here is about 57 m
position_close = 88 # TO DO: write the positions in a file in assignment5_P and read them here
position_halfway = 67
position_far = 57
#position_close = np.loadtxt("position_close.txt", dtype='f', delimiter=' ')
#position_halfway = np.loadtxt("positions_halfway.txt", dtype='f', delimiter=' ')
#position_far = np.loadtxt("positions_far.txt", dtype='f', delimiter=' ')

# If we are also interested in the means
# means_close = np.loadtxt("means_close.txt", dtype='f', delimiter=' ')
# means_halfway = np.loadtxt("means_halfway.txt", dtype='f', delimiter=' ')
# means_far = np.loadtxt("means_far.txt", dtype='f', delimiter=' ')

meanstds_close = np.array([np.mean(stds_close[:,0]), np.mean(stds_close[:,1]), np.mean(stds_close[:,2]), np.mean(stds_close[:,3])])
meanstds_halfway = np.array([np.mean(stds_halfway[:,0]), np.mean(stds_halfway[:,1]), np.mean(stds_halfway[:,2]), np.mean(stds_halfway[:,3])])
meanstds_far = np.array([np.mean(stds_far[:,0]), np.mean(stds_far[:,1]), np.mean(stds_far[:,2]), np.mean(stds_far[:,3])])



def getTfromConcentration(c, x, M, n, D, q, Lambda):
    func = lambda t: M / np.sqrt(4 * np.pi * D * t) * np.exp(- (x - q * t / n) ** 2 / (4 * D * t) - Lambda * t) - c
    return fsolve(func, x0 = 60)

def visualisation(position):

    results = spotpy.analyser.load_csv_results('HydrologyDREAM')
    param_all = results[['parn', 'parD', 'parq', 'parLambda']]
    param = param_all[-100:]
    param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))
    n, D, q, Lambda = param[15,0], param[15,1], param[15,2], param[15,3] #we should randomize this, but for the moment consider 15 random:)
    tcrit = np.floor(getTfromConcentration(2.5, 100, 200, n, D, q, Lambda)[0])
    tdet = tcrit - 10
    h = (tdet - 36) / 4
    obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

    if position == 'close':
        pos = position_close
    if position == 'halfway':
        pos = position_halfway
    if position == 'far':
        pos = position_far

    obs = c(pos, obstimes, 200, n, D, q, Lambda)
    new_meas = [pos, obstimes, obs]
    dream_results_5(new_meas, plot = True)

visualisation('far') #choose position: 'close', 'halfway' or 'far'

labels = ['n', 'D', 'q', 'Lambda']
l = np.arange(len(labels))
            
width = 0.2
plt.bar(l - width, meanstds_close, width, label = 'Close')
plt.bar(l, meanstds_halfway, width, label = 'Halfway')
plt.bar(l + width, meanstds_far, width, label = 'Far')
plt.title('Standard deviations of parameters depending on the position of the new monitoring well relative to the drinking well')
plt.xticks(l, labels = labels)
plt.legend()
plt.show()