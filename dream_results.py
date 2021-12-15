import numpy as np
import spotpy
from dream_setup import spotpy_setup
from dream_run import dream_run

spot_setup = spotpy_setup()
results = dream_run()
spotpy.analyser.plot_parametertrace_algorithms(results, 'dream', spot_setup)