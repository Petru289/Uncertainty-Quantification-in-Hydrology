# Global sensitivity analysis using the SAlib library.
# Because of the required sample size this can take a while.

import enum
from SALib.sample import saltelli
from SALib.analyze import sobol
import sys
sys.path.insert(0, './Assignment1')
import numpy as np
import parameters as par
from functions import c
import matplotlib.pylab as plt
import seaborn as sns

# Define the model input
problem = {
    'num_vars': 4,
    'names': ['n', 'D', 'q', 'Lambda'],
    'bounds': [[par.nRange[0], par.nRange[1]],
                [par.DRange[0], par.DRange[1]],
                [par.qRange[0], par.qRange[1]],
                [par.LambdaRange[0], par.LambdaRange[1]]]
}

# Because we use a saltelli sampler the sample size will be N (2D + 2), with D being the number of input parameters
N = 2**16

# Prepare the figure
fig, axes = plt.subplots(nrows = len(par.x), ncols = len(par.t), figsize=(8.0*1.4, 5.0*1.4))
fig.suptitle('Second order sobol indices', fontsize=20)
fig.subplots_adjust(left=0.05, bottom = 0.1, right = 0.95, top = 0.9, wspace = 0.1, hspace = 0.6)
width = 0.3
axisLabelFontSize = 10


for xIndex, x in enumerate(par.x):
    for tIndex, t in enumerate(par.t):

        # Sample
        param_values = saltelli.sample(problem, N)
        Y = np.zeros([param_values.shape[0]])

        # Evaluate the model
        for i, X in enumerate(param_values):
            Y[i] = c(x, t, par.M ,X[0],X[1],X[2],X[3])

        # Analyze the model
        Si = sobol.analyze(problem, Y, print_to_console=False)

        axes[xIndex, tIndex].set_title('x = {}, t = {}'.format(x, t), fontsize = 12)
        sns.heatmap(np.around(Si['S2'], 4), ax=axes[xIndex,tIndex], linewidth=0.5, cmap="YlGnBu", xticklabels=['n','D','q','$\lambda$'], yticklabels=['n','D','q','$\lambda$'], annot=True, linewidths=2, linecolor='black', vmin=0, vmax=0.5)

# This migth not work on windows because it uses a different figure manager
mng = plt.get_current_fig_manager()
mng.full_screen_toggle()
plt.savefig('salib.png', dpi=300)


#plt.show()
