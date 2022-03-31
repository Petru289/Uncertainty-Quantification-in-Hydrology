import numpy as np
import chaospy as cp
from numpy.core.fromnumeric import size
import matplotlib.pyplot as plt
from functions import c
import parameters as par


x = par.x
t = par.t

fig, axes = plt.subplots(nrows = len(x), ncols = len(t))
fig.suptitle('Sobol indices', fontsize=20)
width = 0.3
axisLabelFontSize = 10
#fig2, axes2 = plt.subplots(nrows = len(x), ncols = len(t))
#fig2.tight_layout()

file = open('output2.txt', 'w')

#for N in [1000, 10000, 100000]: #sample size
N = 2**16

print("Sample size =", N)
print()
file.write("Sample size = ")
file.write(str(N))
file.write('\n')
file.write('\n')
file.write('\n')

for distIndex, distance in enumerate(x):
    for timeIndex, timestep in enumerate(t):
        # Get Matrix A 
        n = cp.Uniform(par.nRange[0], par.nRange[1])
        D = cp.Uniform(par.DRange[0], par.DRange[1])
        q = cp.Uniform(par.qRange[0], par.qRange[1])
        Lambda = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
        jointdist = cp.J(n, D, q, Lambda)
        A = np.transpose(jointdist.sample(size = N))
        #transpose, so that it's the same as in the algorithm description

        #print(A.shape)
        #print(np.max(A, axis = 0))
        #print(np.min(A, axis = 0))

        # Get matrix B
        n2 = cp.Uniform(par.nRange[0], par.nRange[1])
        D2 = cp.Uniform(par.DRange[0], par.DRange[1])
        q2 = cp.Uniform(par.qRange[0], par.qRange[1])
        Lambda2 = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
        jointdist2 = cp.J(n2, D2, q2, Lambda2)
        B = np.transpose(jointdist2.sample(size = N))

        
        #https://stackoverflow.com/questions/3059395/numpy-array-assignment-problem

        Cn = B.copy()
        CD = B.copy()
        Cq = B.copy()
        CLambda = B.copy()

        Cn[:,0] = A[:,0]
        CD[:,1] = A[:,1]
        Cq[:,2] = A[:,2]
        CLambda[:,3] = A[:,3]

        
        # Get concentrations
        yA = c(distance, timestep, par.M, A[:,0], A[:,1], A[:,2], A[:,3])
        yB = c(distance, timestep, par.M, B[:,0], B[:,1], B[:,2], B[:,3])
    

        yCn = c(distance, timestep, par.M, Cn[:,0], Cn[:,1], Cn[:,2], Cn[:,3])
        yCD = c(distance, timestep, par.M, CD[:,0], CD[:,1], CD[:,2], CD[:,3])
        yCq = c(distance, timestep, par.M, Cq[:,0], Cq[:,1], Cq[:,2], Cq[:,3])
        yCLambda = c(distance, timestep, par.M, CLambda[:,0], CLambda[:,1], CLambda[:,2], CLambda[:,3])

        

        #We need nicer notations


        print("x = ", distance, "time = ", timestep)
        file.write("x = ")
        file.write(str(distance))
        file.write(" m")
        file.write(" time = ")
        file.write(str(timestep))
        file.write(" days")
        file.write('\n')
        file.write('\n')

        # First order sobol
        file.write("First order Sobol indices")
        file.write('\n')
        print("First order Sobol indices")
        Sn =  (np.inner(yA, yCn) / N - np.mean(yA) * np.mean(yB)) / (np.mean(yA ** 2) - np.mean(yA) * np.mean(yB))
        file.write("n: ")
        file.write(str(Sn))
        file.write('\n')
        print("n: ", Sn)

        SD = (np.inner(yA, yCD) / N - np.mean(yA) * np.mean(yB)) / (np.mean(yA ** 2) - np.mean(yA) * np.mean(yB))
        file.write("D: ")
        file.write(str(SD))
        file.write('\n')
        print("D: ", SD)

        Sq = (np.inner(yA, yCq) / N - np.mean(yA) * np.mean(yB)) / (np.mean(yA ** 2) - np.mean(yA) * np.mean(yB))
        file.write("q: ")
        file.write(str(Sq))
        file.write('\n')
        print("q: ", Sq)

        SLambda = (np.inner(yA, yCLambda) / N - np.mean(yA) * np.mean(yB)) / (np.mean(yA ** 2) - np.mean(yA) * np.mean(yB))
        file.write("Lambda: ")
        file.write(str(SLambda))
        file.write('\n')
        file.write('\n')
        print("Lambda: ", SLambda)

        file.write("Total Sobol indices")
        file.write('\n')
        print("Total Sobol indices")
        STn = 1 - (np.inner(yB, yCn) / N - np.mean(yA) * np.mean(yB)) / (np.mean(yA ** 2) - np.mean(yA) * np.mean(yB))
        file.write("n: ")
        file.write(str(STn))
        file.write("\n")
        print("n: ", STn)

        STD = 1 - (np.inner(yB, yCD) / N - np.mean(yA) * np.mean(yB)) / (np.mean(yA ** 2) - np.mean(yA) * np.mean(yB))
        file.write("D: ")
        file.write(str(STD))
        file.write("\n")
        print("D: ", STD)

        STq = 1 - (np.inner(yB, yCq) / N - np.mean(yA) * np.mean(yB)) / (np.mean(yA ** 2) - np.mean(yA) * np.mean(yB))
        file.write("q: ")
        file.write(str(STq))
        file.write("\n")
        print("q: ", STq)

        STLambda = 1 - (np.inner(yB, yCLambda) / N - np.mean(yA) * np.mean(yB)) / (np.mean(yA ** 2) - np.mean(yA) * np.mean(yB))
        file.write("LAmbda: ")
        file.write(str(STLambda))
        file.write("\n")
        file.write("\n")
        print("Lambda: ", STLambda)

        print()
        
        #https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html

        labels = ['n', 'D', 'q', 'Lambda']
        sobol = [Sn, SD, Sq, SLambda]
        total_sobol = [STn, STD, STq, STLambda]
        l = np.arange(len(labels))
            

        if np.isnan(Sn) == True:
            axes[distIndex, timeIndex].bar(l - width / 2, [0, 0, 0, 0], width, label = 'First order')
            axes[distIndex, timeIndex].bar(l + width / 2, [0, 0, 0, 0], width, label = 'Total')
            axes[distIndex, timeIndex].set_title('x = {}, t = {}'.format(distance, timestep), fontsize = axisLabelFontSize)
            axes[distIndex, timeIndex].set_xticks(l)
            axes[distIndex, timeIndex].set_xticklabels(labels)
            axes[distIndex, timeIndex].set_ylim([0,1])
            
        else:
            axes[distIndex, timeIndex].bar(l - width / 2, sobol, width, label = 'First order')
            axes[distIndex, timeIndex].bar(l + width / 2, total_sobol, width, label = 'Total')
            axes[distIndex, timeIndex].set_title('x = {}, t = {}'.format(distance, timestep), fontsize = axisLabelFontSize)
            axes[distIndex, timeIndex].set_xticks(l)
            axes[distIndex, timeIndex].set_xticklabels(labels)
            axes[distIndex, timeIndex].set_ylim([0,1])
          

        fig.tight_layout()
axes[0, 2].legend(loc=(0.70,1.1))
plt.show()

file.close()

#Reference: https://en.wikipedia.org/wiki/Variance-based_sensitivity_analysis

# ToDo: clean up code
# ToDo: compare to salib library



