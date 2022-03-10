from tabnanny import check
import numpy as np
import spotpy
from scipy.optimize import fsolve
from functions import c
from dream_results_5 import dream_results_5
from assignment5_P import choose
import matplotlib.pyplot as plt
from scipy.stats import bootstrap, levene





# 'close', 'halfway' and 'far' relative to the position of the drinking well

def chooseposition(position, visualisation = False):

    (new_param, pos) = choose(position)
    stds = np.array([np.std(new_param[:,0]), np.std(new_param[:,1]), np.std(new_param[:,2]), np.std(new_param[:,3])])
    res_n = bootstrap((new_param[:,0], ), np.std, confidence_level=0.95)
    print(res_n.confidence_interval)
    res_D = bootstrap((new_param[:,1], ), np.std, confidence_level=0.95)
    print(res_D.confidence_interval)
    res_q = bootstrap((new_param[:,2], ), np.std, confidence_level=0.95)
    print(res_q.confidence_interval)
    res_Lambda = bootstrap((new_param[:,3], ), np.std, confidence_level=0.95)
    print(res_Lambda.confidence_interval)
    conf_int = [res_n.confidence_interval,res_D.confidence_interval,res_q.confidence_interval,res_Lambda.confidence_interval]


    def getTfromConcentration(c, x, M, n, D, q, Lambda):
        func = lambda t: M / np.sqrt(4 * np.pi * D * t) * np.exp(- (x - q * t / n) ** 2 / (4 * D * t) - Lambda * t) - c
        return fsolve(func, x0 = 60)

    if visualisation == True:
        
        results = spotpy.analyser.load_csv_results('HydrologyDREAM')
        param_all = results[['parn', 'parD', 'parq', 'parLambda']]
        param = param_all[-100:]
        param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))
        # Maybe take the mean?
        n, D, q, Lambda = param[15,0], param[15,1], param[15,2], param[15,3] #we should randomize this, but for the moment consider 15 random:)
        tcrit = np.floor(getTfromConcentration(2.5, 100, 200, n, D, q, Lambda)[0])
        tdet = tcrit - 10
        h = (tdet - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        dream_results_5(new_meas, plot = True)

    return(stds, conf_int, new_param)


results_close = chooseposition('close')
results_halfway = chooseposition('halfway')
results_far = chooseposition('far')
var_close = results_close[0]
var_halfway = results_halfway[0]
var_far = results_far[0]
conf_int_length_close = np.array([results_close[1][i][1] - results_close[1][i][0] for i in range(len(results_close[1]))])
conf_int_length_halfway = np.array([results_halfway[1][i][1] - results_halfway[1][i][0] for i in range(len(results_halfway[1]))])
conf_int_length_far = np.array([results_far[1][i][1] - results_far[1][i][0] for i in range(len(results_far[1]))])
variancetest_n = [levene(results_close[2][:,0], results_halfway[2][:,0]), 
                  levene(results_close[2][:,0], results_far[2][:,0]), levene(results_halfway[2][:,0], results_far[2][:,0])]
variancetest_D = [levene(results_close[2][:,1], results_halfway[2][:,1]), 
                  levene(results_close[2][:,1], results_far[2][:,1]), levene(results_halfway[2][:,1], results_far[2][:,1])]
variancetest_q = [levene(results_close[2][:,2], results_halfway[2][:,2]), 
                  levene(results_close[2][:,2], results_far[2][:,2]), levene(results_halfway[2][:,2], results_far[2][:,2])]
variancetest_Lambda = [levene(results_close[2][:,3], results_halfway[2][:,3]), 
                  levene(results_close[2][:,3], results_far[2][:,3]), levene(results_halfway[2][:,3], results_far[2][:,3])]

variancetest_n_3 = levene(results_close[2][:,0], results_halfway[2][:,0], results_far[2][:,0])
variancetest_D_3 = levene(results_close[2][:,1], results_halfway[2][:,1], results_far[2][:,1])
variancetest_q_3 = levene(results_close[2][:,2], results_halfway[2][:,2], results_far[2][:,2])
variancetest_Lambda_3 = levene(results_close[2][:,3], results_halfway[2][:,3], results_far[2][:,3])

print("Test whether the variances are significantly different:")
print(" ")
print("Pairwise: ")
print(" ")
print("n: ")
print(" ")
print(variancetest_n)
print(" ")
print("D: ")
print(" ")
print(variancetest_D)
print(" ")
print("q: ")
print(" ")
print(variancetest_q)
print(" ")
print("Lambda: ")
print(" ")
print(variancetest_Lambda)
print(" ")
print("Triples:")
print(" ")
print("n: ")
print(" ")
print(variancetest_n_3)
print(" ")
print("D: ")
print(" ")
print(variancetest_D_3)
print(" ")
print("q: ")
print(" ")
print(variancetest_q_3)
print(" ")
print("Lambda: ")
print(" ")
print(variancetest_Lambda_3)
print(" ")



labels = ['n', 'D', 'q', 'Lambda']
l = np.arange(len(labels))
width = 0.2
plt.bar(l - width, var_close, width, label = 'Close', yerr = conf_int_length_close)
plt.bar(l, var_halfway, width, label = 'Halfway', yerr = conf_int_length_halfway)
plt.bar(l + width, var_far, width, label = 'Far', yerr = conf_int_length_far)
plt.title('Standard deviations of parameters depending on the position of the new monitoring well relative to the drinking well')
plt.xticks(l, labels = labels)
plt.legend()
plt.show()



# stds_close = np.loadtxt("stds_close.txt", dtype='f', delimiter=' ') #pos here is about 88 m
# stds_halfway = np.loadtxt("stds_halfway.txt", dtype='f', delimiter=' ') #pos here is about 67 m
# stds_far = np.loadtxt("stds_far.txt", dtype='f', delimiter=' ') #pos here is about 57 m
# position_close = 88 # TO DO: write the positions in a file in assignment5_P and read them here
# position_halfway = 67
# position_far = 57

# meanstds_close = np.array([np.mean(stds_close[:,0]), np.mean(stds_close[:,1]), np.mean(stds_close[:,2]), np.mean(stds_close[:,3])])
# meanstds_halfway = np.array([np.mean(stds_halfway[:,0]), np.mean(stds_halfway[:,1]), np.mean(stds_halfway[:,2]), np.mean(stds_halfway[:,3])])
# meanstds_far = np.array([np.mean(stds_far[:,0]), np.mean(stds_far[:,1]), np.mean(stds_far[:,2]), np.mean(stds_far[:,3])])