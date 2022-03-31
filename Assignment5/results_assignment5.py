from tabnanny import check
import numpy as np
import spotpy
from scipy.optimize import fsolve
import sys
sys.path.insert(0, './Assignment1')
from functions import c
from dream_results_5 import dream_results_5
from assignment5_P import choose
import matplotlib.pyplot as plt
from scipy.stats import bootstrap, levene




# results = spotpy.analyser.load_csv_results('HydrologyDREAM')
# param_all = results[['parn', 'parD', 'parq', 'parLambda']]
# param = param_all[-1000:]
# param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))
# n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
# time2 = np.array([35, 38, 41, 44, 47, 50, 53])
# obs2 = np.zeros(7)
# obs2 = c(50, time2, 200, n, D, q, Lambda)
# print(obs2)



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
        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        tcrit = np.floor(getTfromConcentration(2.5, 100, 200, n, D, q, Lambda)[0])
        tdet = tcrit - 10
        h = (tdet - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        dream_results_5(new_meas, plot = True)

    return(stds, conf_int, new_param)


positions = ['88 - equidistant timepoints', '67 - equidistant timepoints', '57 - equidistant timepoints', '50 - equidistant timepoints', '62 - equidistant timepoints', '73 - equidistant timepoints', '57 - early timepoints', '57 - late timepoints', '67 - early timepoints', '67 - late timepoints','88 - early timepoints', '88 - late timepoints']
positions = ['50 - equidistant timepoints']
results = [chooseposition(positions[i]) for i in range(len(positions))]
var = [results[i][0] for i in range(len(results))]
conf_int_length = []
for j in range(len(results)):
    conf_int_length.append(np.array([results[j][1][i][1] - results[j][1][i][0] for i in range(4)]))
   
results_initial = spotpy.analyser.load_csv_results('HydrologyDREAM')
param_all = results_initial[['parn', 'parD', 'parq', 'parLambda']]
param = param_all[-1000:]
param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))
initial_stds = np.array([np.std(param[:,0]), np.std(param[:,1]), np.std(param[:,2]), np.std(param[:,3])])
print("Initial means:")
print(" ")
print(np.array([np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])]))
print("Initial variances")
print(" ")
print(initial_stds)

labels = ['n', 'D', 'q', 'Lambda']
l = np.arange(len(labels))
width = 0.06
k = len(var)
count = 0
plt.figure(1)
if k % 2 != 0:
    for i in range(-np.int(np.floor(k/2)), np.int(np.floor(k/2 + 1))):
        plt.bar(l + i*width, var[count], width, label = positions[count], yerr = conf_int_length[count])
        count += 1
    plt.bar(l, initial_stds, width/2, color = ["turquoise"], label = 'initial variance')
else:
    rng = [i/2 for i in range(-(k-1), k+1, 2)]
    for i in rng:
        plt.bar(l + i*width, var[count], width, label = positions[count], yerr = conf_int_length[count])
        count += 1
    plt.bar(l, initial_stds, width/2, color = ["turquoise"], label = 'initial variance')

plt.title('Standard deviations of parameters depending on the position of the new monitoring well')
plt.xticks(l, labels = labels)
plt.legend()

count = 0
plt.figure(2)
if k % 2 != 0:
    for i in range(-np.int(np.floor(k/2)), np.int(np.floor(k/2 + 1))):  
        plt.bar(l + i*width, var[count]/initial_stds, width, label = positions[count])
        count += 1
else:
    rng = [i/2 for i in range(-(k-1), k+1, 2)]
    for i in rng:
        plt.bar(l + i*width, var[count]/initial_stds, width, label = positions[count])
        count += 1

plt.axhline(y=1, color='r', linestyle='-')
plt.title('Percentage of standard deviations with respect to the initial standard deviations')
plt.xticks(l, labels = labels)
plt.legend()

plt.figure(3)
plt.hist(param[:,0], 30, label = "initial distribution")
plt.hist(results[0][2][:,0], 30, label = "distribution after the new measurements at 50 m")
plt.xlabel("n")
plt.legend(loc='upper right')

plt.figure(4)
plt.hist(param[:,1], 30, label = "initial distribution")
plt.hist(results[0][2][:,1], 30, label = "distribution after the new measurements at 50 m")
plt.xlabel("D")
plt.legend(loc='upper right')

plt.figure(5)
plt.hist(param[:,2], 30, label = "initial distribution")
plt.hist(results[0][2][:,2], 30, label = "distribution after the new measurements at 50 m")
plt.xlabel("q")
plt.legend(loc='upper right')

plt.figure(6)
plt.hist(param[:,3], 30, label = "initial distribution")
plt.hist(results[0][2][:,3], 30, label = "distribution after the new measurements at 50 m")
plt.xlabel("Lambda")
plt.legend(loc='upper right')


plt.show()


print(2)




# plt.bar(l - width, var_close, width, label = 'Close', yerr = conf_int_length_close)
# plt.bar(l, var_halfway, width, label = 'Halfway', yerr = conf_int_length_halfway)
# plt.bar(l + width, var_far, width, label = 'Far', yerr = conf_int_length_far)
# plt.title('Standard deviations of parameters depending on the position of the new monitoring well relative to the drinking well')
# plt.xticks(l, labels = labels)
# plt.legend()
# plt.show()

# results_close = chooseposition('close')
# results_halfway = chooseposition('halfway')
# results_far = chooseposition('far')
# results_50 = chooseposition('50')
# var_close = results_close[0]
# var_halfway = results_halfway[0]
# var_far = results_far[0]
# var_50 = results_50[0]
# conf_int_length_close = np.array([results_close[1][i][1] - results_close[1][i][0] for i in range(len(results_close[1]))])
# conf_int_length_halfway = np.array([results_halfway[1][i][1] - results_halfway[1][i][0] for i in range(len(results_halfway[1]))])
# conf_int_length_far = np.array([results_far[1][i][1] - results_far[1][i][0] for i in range(len(results_far[1]))])
# conf_int_length_50 = np.array([results_50[1][i][1] - results_50[1][i][0] for i in range(len(results_50[1]))])


# variancetest_n = [levene(results_close[2][:,0], results_halfway[2][:,0]), 
#                   levene(results_close[2][:,0], results_far[2][:,0]), levene(results_halfway[2][:,0], results_far[2][:,0])]
# variancetest_D = [levene(results_close[2][:,1], results_halfway[2][:,1]), 
#                   levene(results_close[2][:,1], results_far[2][:,1]), levene(results_halfway[2][:,1], results_far[2][:,1])]
# variancetest_q = [levene(results_close[2][:,2], results_halfway[2][:,2]), 
#                   levene(results_close[2][:,2], results_far[2][:,2]), levene(results_halfway[2][:,2], results_far[2][:,2])]
# variancetest_Lambda = [levene(results_close[2][:,3], results_halfway[2][:,3]), 
#                   levene(results_close[2][:,3], results_far[2][:,3]), levene(results_halfway[2][:,3], results_far[2][:,3])]

# variancetest_n_3 = levene(results_close[2][:,0], results_halfway[2][:,0], results_far[2][:,0])
# variancetest_D_3 = levene(results_close[2][:,1], results_halfway[2][:,1], results_far[2][:,1])
# variancetest_q_3 = levene(results_close[2][:,2], results_halfway[2][:,2], results_far[2][:,2])
# variancetest_Lambda_3 = levene(results_close[2][:,3], results_halfway[2][:,3], results_far[2][:,3])

# print("Test whether the variances are significantly different:")
# print(" ")
# print("Pairwise: ")
# print(" ")
# print("n: ")
# print(" ")
# print(variancetest_n)
# print(" ")
# print("D: ")
# print(" ")
# print(variancetest_D)
# print(" ")
# print("q: ")
# print(" ")
# print(variancetest_q)
# print(" ")
# print("Lambda: ")
# print(" ")
# print(variancetest_Lambda)
# print(" ")
# print("Triples:")
# print(" ")
# print("n: ")
# print(" ")
# print(variancetest_n_3)
# print(" ")
# print("D: ")
# print(" ")
# print(variancetest_D_3)
# print(" ")
# print("q: ")
# print(" ")
# print(variancetest_q_3)
# print(" ")
# print("Lambda: ")
# print(" ")
# print(variancetest_Lambda_3)
# print(" ")









#SCRAP PAPER

# stds_close = np.loadtxt("stds_close.txt", dtype='f', delimiter=' ') #pos here is about 88 m
# stds_halfway = np.loadtxt("stds_halfway.txt", dtype='f', delimiter=' ') #pos here is about 67 m
# stds_far = np.loadtxt("stds_far.txt", dtype='f', delimiter=' ') #pos here is about 57 m
# position_close = 88 # TO DO: write the positions in a file in assignment5_P and read them here
# position_halfway = 67
# position_far = 57

# meanstds_close = np.array([np.mean(stds_close[:,0]), np.mean(stds_close[:,1]), np.mean(stds_close[:,2]), np.mean(stds_close[:,3])])
# meanstds_halfway = np.array([np.mean(stds_halfway[:,0]), np.mean(stds_halfway[:,1]), np.mean(stds_halfway[:,2]), np.mean(stds_halfway[:,3])])
# meanstds_far = np.array([np.mean(stds_far[:,0]), np.mean(stds_far[:,1]), np.mean(stds_far[:,2]), np.mean(stds_far[:,3])])