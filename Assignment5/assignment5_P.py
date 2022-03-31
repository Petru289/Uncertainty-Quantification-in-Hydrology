import spotpy
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from dream_results_5 import dream_results_5
import sys
sys.path.insert(0, './Assignment1')
from functions import c


def choose(position):

    results = spotpy.analyser.load_csv_results('HydrologyDREAM')
    param_all = results[['parn', 'parD', 'parq', 'parLambda']]
    param = param_all[-1000:]
    param = np.column_stack((param['parn'], param['parD'], param['parq'], param['parLambda']))


    time = np.arange(1,101)
    cmax = np.zeros(param.shape[0])

    for i in range(len(cmax)):
        cmax[i] = np.max(c(100, time, 200, param[i,0], param[i,1], param[i,2], param[i,3]))

    #Pick only those values of cmax that are within one standard deviation from the mean: cmax \in (mean - std, mean + std)
    mean = np.mean(cmax)
    std =  np.std(cmax)
    ind = np.abs(cmax - mean) < std
    param = param[ind]


    def getTfromConcentration(c, x, M, n, D, q, Lambda):
        func = lambda t: M / np.sqrt(4 * np.pi * D * t) * np.exp(- (x - q * t / n) ** 2 / (4 * D * t) - Lambda * t) - c
        return fsolve(func, x0 = 60)

    def getXfromConcentration(c, t, M, n, D, q ,Lambda):
        A = M / (np.sqrt(4 * np.pi * D * t))
        x1 = q * t / n - np.sqrt(- 4 * D * t * (np.log(c / A) + Lambda * t))
        x2 = q * t / n + np.sqrt(- 4 * D * t * (np.log(c / A) + Lambda * t))
        return [x1, x2]

    limittimes = np.zeros(param.shape[0])
    positions = np.zeros(param.shape[0])

    #x where contaminant is detectable at day 36, i.e. earliest you can measure something -- gives farthest point (i.e. closest to the drinking well)
    #we can do measurements in the allowed time interval only before the center of mass passed the sampling point
    if position == '88 - equidistant timepoints':

        for i in range(param.shape[0]):
            tcrit = np.floor(getTfromConcentration(2.5, 100, 200, param[i,0], param[i,1], param[i,2], param[i,3])[0])
            tdet = tcrit - 10
            x_s = getXfromConcentration(10 ** (-5), 36, 200, param[i,0], param[i,1], param[i,2], param[i,3])
            limittimes[i] = tdet
            positions[i] = x_s[1]
        
        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = np.floor(np.mean(positions))
        pos = 88 #Sometimes pos varies with difference of 1, so, for stability of the results, we fix it here
        #Divide the interval to obtain 5 points in total
        h = (tdet - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)



    #x where contaminant is 0.5 at day 36. We want the maximum to be attained at the half time between 36 and tcrit - 10,
    #so that we can do measurements in that time interval before and after the center of mass passed the sampling point
    #we can figure out from previous simulations that this approximately happens when contaminant is about 0.5 at day 36
    if position == '67 - equidistant timepoints':

        for i in range(param.shape[0]):
            tcrit = np.floor(getTfromConcentration(2.5, 100, 200, param[i,0], param[i,1], param[i,2], param[i,3])[0])
            tdet = tcrit - 10
            x_s = getXfromConcentration(0.5, 36, 200, param[i,0], param[i,1], param[i,2], param[i,3])
            limittimes[i] = tdet
            positions[i] = x_s[1]

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = np.floor(np.mean(positions))
        pos = 67 # Sometimes pos varies with difference of 1, so, for stability of the results, we fix it here
        #Divide the interval to obtain 5 points in total
        h = (tdet - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)



#x where contaminant is 0.1 at day tcrit - 10. We want the maximum to be attained at approximately day 36,
#so that we can do measurements only after the center of mass passed the sampling point
#we can figure out from previous simulations that this approximately happens when contaminant is about 0.1 at day tcrit - 10
    if position == '57 - equidistant timepoints':

        for i in range(param.shape[0]):
            tcrit = np.floor(getTfromConcentration(2.5, 100, 200, param[i,0], param[i,1], param[i,2], param[i,3])[0])
            tdet = tcrit - 10
            x_s = getXfromConcentration(0.1, tdet, 200, param[i,0], param[i,1], param[i,2], param[i,3])

            limittimes[i] = tdet
            positions[i] = x_s[0]

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = np.floor(np.mean(positions))
        pos = 57 #Sometimes pos varies with difference of 1, so, for stability of the results, we fix it here
        #Divide the interval to obtain 5 points in total
        h = (tdet - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)

  

    if position == '50 - equidistant timepoints':

        for i in range(param.shape[0]):
            tcrit = np.floor(getTfromConcentration(2.5, 100, 200, param[i,0], param[i,1], param[i,2], param[i,3])[0])
            tdet = tcrit - 10
            x_s = getXfromConcentration(0.1, tdet, 200, param[i,0], param[i,1], param[i,2], param[i,3])

            limittimes[i] = tdet
            #positions[i] = x_s[0]

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = 50
        #Divide the interval to obtain 5 points in total
        h = (tdet - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)

    if position == '88 - late timepoints':

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = 88
        #Divide the interval to obtain 5 points in total
        h = (tdet - 46) / 4
        obstimes = np.array([46, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)

    if position == '57 - late timepoints':

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = 57
        #Divide the interval to obtain 5 points in total
        h = (tdet - 46) / 4
        obstimes = np.array([46, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)

    if position == '67 - late timepoints':

        for i in range(param.shape[0]):
            tcrit = np.floor(getTfromConcentration(2.5, 100, 200, param[i,0], param[i,1], param[i,2], param[i,3])[0])
            tdet = tcrit - 10
            x_s = getXfromConcentration(0.1, tdet, 200, param[i,0], param[i,1], param[i,2], param[i,3])

            limittimes[i] = tdet
            #positions[i] = x_s[0]

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = 67
        #Divide the interval to obtain 5 points in total
        h = (tdet - 46) / 4
        obstimes = np.array([46, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)

    if position == '67 - early timepoints':

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = 67
        #Divide the interval to obtain 5 points in total
        h = (tdet - 10 - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet - 10])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)
    
    if position == '57 - early timepoints':

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = 57
        #Divide the interval to obtain 5 points in total
        h = (tdet - 10 - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet - 10])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)

    if position == '88 - early timepoints':

        tdet = np.floor(np.mean(limittimes))
        tdet = 55 # We fix it here to make sure of the stability
        pos = 88
        #Divide the interval to obtain 5 points in total
        h = (tdet - 10 - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet - 10])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)
    
    
    if position == '62 - equidistant timepoints':

        tdet = 55 # We fix it here to make sure of the stability
        pos = 62
        #Divide the interval to obtain 5 points in total
        h = (tdet - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)

    if position == '73 - equidistant timepoints':

        tdet = 55 # We fix it here to make sure of the stability
        pos = 73
        #Divide the interval to obtain 5 points in total
        h = (tdet - 36) / 4
        obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])

        n, D, q, Lambda = np.mean(param[:,0]), np.mean(param[:,1]), np.mean(param[:,2]), np.mean(param[:,3])
        obs = c(pos, obstimes, 200, n, D, q, Lambda)
        new_meas = [pos, obstimes, obs]
        new_param = dream_results_5(new_meas)

    return(new_param, pos)
        
    
if __name__ == "main":
    choose('halfway')




# SCRAP PAPER



# for i in range(param.shape[0]):
#     concentration = c(100, time, 200, param[i,0], param[i,1], param[i,2], param[i,3])
#     plt.plot(time, concentration, alpha = 0.2, color = 'navy')

# plt.ylabel('concentration $[kg/m^3]$')
# plt.xlabel('time [days]')
# plt.show()

# mat = np.matrix(stds)
# with open('stds_far.txt','wb') as f:
#     for line in mat:
#         np.savetxt(f, line, fmt='%.2f')

# mat = np.matrix(means)
# with open('means_far.txt','wb') as f:
#     for line in mat:
#         np.savetxt(f, line, fmt='%.2f')

# with open('position_far.txt','wb') as f:
#     np.savetxt(f, pos, fmt='%.2f')




#     #x where contaminant is detectable at day 36, i.e. earliest you can measure something -- gives farthest point (i.e. closest to the drinking well)
#     #we can do measurements in the allowed time interval only before the center of mass passed the sampling point
#     if position == 'close':

#         for i in range(param.shape[0]):
#             tcrit = np.floor(getTfromConcentration(2.5, 100, 200, param[i,0], param[i,1], param[i,2], param[i,3])[0])
#             tdet = tcrit - 10
#             x_s = getXfromConcentration(10 ** (-5), 36, 200, param[i,0], param[i,1], param[i,2], param[i,3])
#             limittimes[i] = tdet
#             positions[i] = x_s[1]
        
#         tdet = np.floor(np.mean(limittimes))
#         pos = np.floor(np.mean(positions))
#         #Divide the interval to obtain 5 points in total
#         h = (tdet - 36) / 4
#         obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])
#         means = np.zeros(param.shape)
#         stds = np.zeros(param.shape)

#         for i in range(param.shape[0]): #param.shape[0]
#             n, D, q, Lambda = param[i,0], param[i,1], param[i,2], param[i,3]
#             obs = c(pos, obstimes, 200, n, D, q, Lambda)
#             new_meas = [pos, obstimes, obs]
#             new_param = dream_results_5(new_meas)
#             #means[i,:] = np.array([np.mean(new_param[:,j]) for j in range(means.shape[1])])
#             stds[i,:] = np.array([np.std(new_param[:,j]) for j in range(stds.shape[1])])
        
#         mat = np.matrix(stds)
#         with open('stds_close.txt','wb') as f:
#             for line in mat:
#                 np.savetxt(f, line, fmt='%.2f')
        
#         mat = np.matrix(means)
#         with open('means_close.txt','wb') as f:
#             for line in mat:
#                 np.savetxt(f, line, fmt='%.2f')

#         with open('position_close.txt','wb') as f:
#             np.savetxt(f, pos, fmt='%.2f')

#     #x where contaminant is 0.5 at day 36. We want the maximum to be attained at the half time between 36 and tcrit - 10,
#     #so that we can do measurements in that time interval before and after the center of mass passed the sampling point
#     #we can figure out from previous simulations that this approximately happens when contaminant is about 0.5 at day 36
#     if position == 'halfway':

#         for i in range(param.shape[0]):
#             tcrit = np.floor(getTfromConcentration(2.5, 100, 200, param[i,0], param[i,1], param[i,2], param[i,3])[0])
#             tdet = tcrit - 10
#             x_s = getXfromConcentration(0.5, 36, 200, param[i,0], param[i,1], param[i,2], param[i,3])
#             limittimes[i] = tdet
#             positions[i] = x_s[1]

#         tdet = np.floor(np.mean(limittimes))
#         pos = np.floor(np.mean(positions))
#         #Divide the interval to obtain 5 points in total
#         h = (tdet - 36) / 4
#         obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])
#         means = np.zeros(param.shape)
#         stds = np.zeros(param.shape)
        
#         for i in range(param.shape[0]): #param.shape[0]
#             n, D, q, Lambda = param[i,0], param[i,1], param[i,2], param[i,3]
#             obs = c(pos, obstimes, 200, n, D, q, Lambda)
#             new_meas = [pos, obstimes, obs]
#             new_param = dream_results_5(new_meas)
#             means[i,:] = np.array([np.mean(new_param[:,j]) for j in range(means.shape[1])])
#             stds[i,:] = np.array([np.std(new_param[:,j]) for j in range(stds.shape[1])])
        
#         mat = np.matrix(stds)
#         with open('stds_halfway.txt','wb') as f:
#             for line in mat:
#                 np.savetxt(f, line, fmt='%.2f')
        
#         mat = np.matrix(means)
#         with open('means_halfway.txt','wb') as f:
#             for line in mat:
#                 np.savetxt(f, line, fmt='%.2f')

#         with open('position_halfway.txt','wb') as f:
#             np.savetxt(f, pos, fmt='%.2f')

# #x where contaminant is 0.1 at day tcrit - 10. We want the maximum to be attained at approximately day 36,
# #so that we can do measurements only after the center of mass passed the sampling point
# #we can figure out from previous simulations that this approximately happens when contaminant is about 0.1 at day tcrit - 10
#     if position == 'far':

#         for i in range(param.shape[0]):
#             tcrit = np.floor(getTfromConcentration(2.5, 100, 200, param[i,0], param[i,1], param[i,2], param[i,3])[0])
#             tdet = tcrit - 10
#             x_s = getXfromConcentration(0.1, tdet, 200, param[i,0], param[i,1], param[i,2], param[i,3])

#             limittimes[i] = tdet
#             positions[i] = x_s[0]

#         tdet = np.floor(np.mean(limittimes))
#         pos = np.floor(np.mean(positions))
#         #Divide the interval to obtain 5 points in total
#         h = (tdet - 36) / 4
#         obstimes = np.array([36, np.floor(36 + h), np.floor(36 + 2*h), np.floor(36 + 3*h), tdet])
#         means = np.zeros(param.shape)
#         stds = np.zeros(param.shape)

#         for i in range(param.shape[0]): #param.shape[0]
#             n, D, q, Lambda = param[i,0], param[i,1], param[i,2], param[i,3]
#             obs = c(pos, obstimes, 200, n, D, q, Lambda)
#             new_meas = [pos, obstimes, obs]
#             new_param = dream_results_5(new_meas)
#             means[i,:] = np.array([np.mean(new_param[:,j]) for j in range(means.shape[1])])
#             stds[i,:] = np.array([np.std(new_param[:,j]) for j in range(stds.shape[1])])

#         mat = np.matrix(stds)
#         with open('stds_far.txt','wb') as f:
#             for line in mat:
#                 np.savetxt(f, line, fmt='%.2f')
        
#         mat = np.matrix(means)
#         with open('means_far.txt','wb') as f:
#             for line in mat:
#                 np.savetxt(f, line, fmt='%.2f')

#         with open('position_far.txt','wb') as f:
#             np.savetxt(f, pos, fmt='%.2f')






#These could go in the for loop if we want to print more info
# print("tcrit - 10:")
# print(tdet)
# print("x where contaminant is detectable (10^(-5)) at day 36:")
# print(x_s[1])
# conc1 = c(x_s[1], 36, 200, param[i,0], param[i,1], param[i,2], param[i,3]) #should give 10 ** (-5)
# conc2 = c(x_s[1], tdet, 200, param[i,0], param[i,1], param[i,2], param[i,3])
# print("Concentration on time t = 36:")
# print(conc1)
# print("Concentration on time tcrit - 10:")
# print(conc2)
# print("Max:")
# print(np.max(c(x_s[1], time, 200, param[i,0], param[i,1], param[i,2], param[i,3])))
# print("Day of maximum:")
# print(np.argmax(c(x_s[1], time, 200, param[i,0], param[i,1], param[i,2], param[i,3])))
# print(" ")


# print("tcrit - 10:")
# print(tdet)
# print("x where contaminant is 0.5 at day 36:")
# print(x_s[1])
# conc1 = c(x_s[1], 36, 200, param[i,0], param[i,1], param[i,2], param[i,3]) #should give 10 ** (-1)
# conc2 = c(x_s[1], tdet, 200, param[i,0], param[i,1], param[i,2], param[i,3])
# print("Concentration on time t = 36:")
# print(conc1)
# print("Concentration on time tcrit - 10:")
# print(conc2)
# print("Max:")
# print(np.max(c(x_s[1], time, 200, param[i,0], param[i,1], param[i,2], param[i,3])))
# print("Day of maximum:")
# print(np.argmax(c(x_s[1], time, 200, param[i,0], param[i,1], param[i,2], param[i,3])))
# print("Middle time between 36 and tcrit-10:")
# tmiddle = 36 + np.floor((tdet - 36) / 2)
# print(tmiddle)
# print(" ")


# print("tcrit - 10:")
# print(tdet)
# print("x where contaminant is 0.5 at day tcrit - 10:")
# print(x_s[0])
# conc1 = c(x_s[0], 36, 200, param[i,0], param[i,1], param[i,2], param[i,3]) 
# conc2 = c(x_s[0], tdet, 200, param[i,0], param[i,1], param[i,2], param[i,3]) #should give 0.5
# print("Concentration on time t = 36:")
# print(conc1)
# print("Concentration on time tcrit - 10:")
# print(conc2)
# print("Max:")
# print(np.max(c(x_s[0], time, 200, param[i,0], param[i,1], param[i,2], param[i,3])))
# print("Day of maximum:")
# print(np.argmax(c(x_s[0], time, 200, param[i,0], param[i,1], param[i,2], param[i,3])))
# print(" ")
