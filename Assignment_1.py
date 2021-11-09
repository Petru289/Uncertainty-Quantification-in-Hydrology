import numpy as np
import chaospy as cp
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def c(x, t, M, n, D, q, Lambda):
    c1 = M / np.sqrt(4 * np.pi * D * t)
    c2 = np.exp(-(x - q * t / n) ** 2 / (4 * D * t))
    c3 = np.exp(-Lambda * t)
    return c1 * c2 * c3


plt.rc('xtick', labelsize=8) 
plt.rc('ytick', labelsize=8)

x = [5, 50, 100]
t = [7, 50, 100]
M = 200
rownumber = 0
concentrations = np.zeros(shape = (9, 100))
fig, axes = plt.subplots(nrows = 3, ncols = 3)
plt.subplots_adjust(left=0.1, bottom = 0.1, right = 0.9, top = 0.9, wspace = 0.7, hspace = 0.8)
#fig.tight_layout()
pp = PdfPages('Plots.pdf')
for i in range(0, 3):
    for j in range(0, 3):
        n = cp.Uniform(0.3, 0.5)
        D = cp.Uniform(0.1, 0.7)
        q = cp.Uniform(0.15, 0.5)
        Lambda = cp.Uniform(0, 0.03)
        jointdist = cp.J(n, D, q, Lambda)
        jointsample = jointdist.sample(size = 100)
        concentrations[rownumber, :] = c(x[i], t[j], M, jointsample[0,:], jointsample[1,:], jointsample[2,:], jointsample[3,:])
        mean = np.round(cp.mean(concentrations[rownumber, :]), 2)
        std =  np.round(np.std(concentrations[rownumber, :]), 2)
        print(cp.mean(concentrations[rownumber, :]), np.std(concentrations[rownumber, :]))
        axes[i, j].hist(concentrations[rownumber,:], 30)
        axes[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes[i, j].set_xlabel('concentration (mean = {}, sd = {})'.format(mean, std), fontsize = 8)
        rownumber += 1
fig.savefig(pp, format = "pdf")
plt.show()
#plt.savefig('plots.pdf')




#Sensitivity Analysis


#n varies

def dcdn(x, t, M, n, D, q, Lambda):
    return c(x, t, M, n, D, q, Lambda) * q * (- x + q * t / n)/(2 * D * n ** 2)

D = (0.1 + 0.7)/2
q = (0.15 + 0.5)/2
Lambda = (0 + 0.03)/2
rownumber = 0
concentrations = np.zeros(shape = (9, 100))
dconcentrations = np.zeros(shape = (9, 100))
fig1, axes1 = plt.subplots(nrows = 3, ncols = 3)
fig1.tight_layout()
fig2, axes2 = plt.subplots(nrows = 3, ncols = 3)
fig2.tight_layout()
fig3, axes3 = plt.subplots(nrows = 3, ncols = 3)
fig3.tight_layout()

for i in range(0, 3):
    for j in range(0, 3):
        n = cp.Uniform(0.3, 0.5)
        n_sample = n.sample(size = 100)
        concentrations[rownumber, :] = c(x[i], t[j], M, n_sample, D, q, Lambda)
        dconcentrations[rownumber, :] = dcdn(x[i], t[j], M, n_sample, D, q, Lambda)
        mean = np.round(cp.mean(concentrations[rownumber, :]), 2)
        std =  np.round(np.std(concentrations[rownumber, :]), 2)
        print(cp.mean(concentrations[rownumber, :]), np.std(concentrations[rownumber, :]))
        axes1[i, j].hist(concentrations[rownumber,:], 30)
        axes1[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes1[i, j].set_xlabel('concentration, mean = {}, sd = {}'.format(mean, std), fontsize = 8)

        axes2[i, j].scatter(n_sample, concentrations[rownumber,:], s = 5, c = "firebrick")
        axes2[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes2[i, j].set_ylabel('concentration', fontsize = 8)
        axes2[i, j].set_xlabel('porosity (n)', fontsize = 8)

        axes3[i, j].scatter(n_sample, dconcentrations[rownumber,:], s = 5, c = "darkorange")
        axes3[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes3[i, j].set_ylabel('dc', fontsize = 8)
        axes3[i, j].set_xlabel('porosity (n)', fontsize = 8)

        rownumber += 1

fig1.savefig(pp, format = "pdf")
fig2.savefig(pp, format = "pdf")
fig3.savefig(pp, format = "pdf")
plt.show()




#D varies

def dcdD(x, t, M, n, D, q, Lambda): #to be computed
    return 

n = (0.3 + 0.5)/2
q = (0.15 + 0.5)/2
Lambda = (0 + 0.03)/2
rownumber = 0
concentrations = np.zeros(shape = (9, 100))
dconcentrations = np.zeros(shape = (9, 100))
fig4, axes4 = plt.subplots(nrows = 3, ncols = 3)
fig4.tight_layout()
fig5, axes5 = plt.subplots(nrows = 3, ncols = 3)
fig5.tight_layout()
fig6, axes6 = plt.subplots(nrows = 3, ncols = 3)
fig6.tight_layout()

for i in range(0, 3):
    for j in range(0, 3):
        D = cp.Uniform(0.1, 0.7)
        D_sample = D.sample(size = 100)
        concentrations[rownumber, :] = c(x[i], t[j], M, n, D_sample, q, Lambda)
        dconcentrations[rownumber, :] = dcdD(x[i], t[j], M, n, D_sample, q, Lambda)
        mean = np.round(cp.mean(concentrations[rownumber, :]), 2)
        std =  np.round(np.std(concentrations[rownumber, :]), 2)
        print(cp.mean(concentrations[rownumber, :]), np.std(concentrations[rownumber, :]))
        axes4[i, j].hist(concentrations[rownumber,:], 30)
        axes4[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes4[i, j].set_xlabel('concentration, mean = {}, sd = {}'.format(mean, std), fontsize = 8)

        axes5[i, j].scatter(D_sample, concentrations[rownumber,:], s = 5, c = "firebrick")
        axes5[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes5[i, j].set_ylabel('concentration', fontsize = 8)
        axes5[i, j].set_xlabel('dispersion (D)', fontsize = 8)

        axes6[i, j].scatter(D_sample, dconcentrations[rownumber,:], s = 5, c = "darkorange")
        axes6[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes6[i, j].set_ylabel('dc', fontsize = 8)
        axes6[i, j].set_xlabel('dispersion (D)', fontsize = 8)

        rownumber += 1

fig4.savefig(pp, format = "pdf")
fig5.savefig(pp, format = "pdf")
fig6.savefig(pp, format = "pdf")
plt.show()






#q varies

def dcdq(x, t, M, n, D, q, Lambda): #to be computed
    return 

n = (0.3 + 0.5) / 2
D = (0.1 + 0.7) / 2
Lambda = (0 + 0.03) / 2
rownumber = 0
concentrations = np.zeros(shape = (9, 100))
dconcentrations = np.zeros(shape = (9, 100))
fig7, axes7 = plt.subplots(nrows = 3, ncols = 3)
fig7.tight_layout()
fig8, axes8 = plt.subplots(nrows = 3, ncols = 3)
fig8.tight_layout()
fig9, axes9 = plt.subplots(nrows = 3, ncols = 3)
fig9.tight_layout()

for i in range(0, 3):
    for j in range(0, 3):
        q = cp.Uniform(0.15, 0.5)
        q_sample = q.sample(size = 100)
        concentrations[rownumber, :] = c(x[i], t[j], M, n, D, q_sample, Lambda)
        dconcentrations[rownumber, :] = dcdq(x[i], t[j], M, n, D, q_sample, Lambda)
        mean = np.round(cp.mean(concentrations[rownumber, :]), 2)
        std =  np.round(np.std(concentrations[rownumber, :]), 2)
        print(cp.mean(concentrations[rownumber, :]), np.std(concentrations[rownumber, :]))
        axes7[i, j].hist(concentrations[rownumber,:], 30)
        axes7[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes7[i, j].set_xlabel('concentration, mean = {}, sd = {}'.format(mean, std), fontsize = 8)

        axes8[i, j].scatter(q_sample, concentrations[rownumber,:], s = 5, c = "firebrick")
        axes8[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes8[i, j].set_ylabel('concentration', fontsize = 8)
        axes8[i, j].set_xlabel('discharge (q)', fontsize = 8)

        axes9[i, j].scatter(q_sample, dconcentrations[rownumber,:], s = 5, c = "darkorange")
        axes9[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes9[i, j].set_ylabel('dc', fontsize = 8)
        axes9[i, j].set_xlabel('discharge (q)', fontsize = 8)

        rownumber += 1

fig7.savefig(pp, format = "pdf")
fig8.savefig(pp, format = "pdf")
fig9.savefig(pp, format = "pdf")
plt.show()




#Lambda varies

def dcdLambda(x, t, M, n, D, q, Lambda): #to be computed
    return 

q = (0.15 + 0.5) / 2
n = (0.3 + 0.5) / 2
D = (0.1 + 0.7) / 2
rownumber = 0
concentrations = np.zeros(shape = (9, 100))
dconcentrations = np.zeros(shape = (9, 100))
fig10, axes10 = plt.subplots(nrows = 3, ncols = 3)
fig10.tight_layout()
fig11, axes11 = plt.subplots(nrows = 3, ncols = 3)
fig11.tight_layout()
fig12, axes12 = plt.subplots(nrows = 3, ncols = 3)
fig12.tight_layout()

for i in range(0, 3):
    for j in range(0, 3):
        Lambda = cp.Uniform(0, 0.03)
        Lambda_sample = Lambda.sample(size = 100)
        concentrations[rownumber, :] = c(x[i], t[j], M, n, D, q, Lambda_sample)
        dconcentrations[rownumber, :] = dcdLambda(x[i], t[j], M, n, D, q, Lambda_sample)
        mean = np.round(cp.mean(concentrations[rownumber, :]), 2)
        std =  np.round(np.std(concentrations[rownumber, :]), 2)
        print(cp.mean(concentrations[rownumber, :]), np.std(concentrations[rownumber, :]))
        axes10[i, j].hist(concentrations[rownumber,:], 30)
        axes10[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes10[i, j].set_xlabel('concentration, mean = {}, sd = {}'.format(mean, std), fontsize = 8)

        axes11[i, j].scatter(Lambda_sample, concentrations[rownumber,:], s = 5, c = "firebrick")
        axes11[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes11[i, j].set_ylabel('concentration', fontsize = 8)
        axes11[i, j].set_xlabel('decay (Lambda)', fontsize = 8)

        axes12[i, j].scatter(Lambda_sample, dconcentrations[rownumber,:], s = 5, c = "darkorange")
        axes12[i, j].set_title('x = {}, t = {}'.format(x[i], t[j]), fontsize = 8)
        axes12[i, j].set_ylabel('dc', fontsize = 8)
        axes12[i, j].set_xlabel('decay (Lambda)', fontsize = 8)

        rownumber += 1

fig10.savefig(pp, format = "pdf")
fig11.savefig(pp, format = "pdf")
fig12.savefig(pp, format = "pdf")
plt.show()



pp.close()



