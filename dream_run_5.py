import numpy as np
import spotpy
from dream_setup_5 import spotpy_setup




def dream_run_5(new_meas):

    results = []
    timeout = 10 #Given in Seconds

    #see spotpy._algorithm.py for comment on this
    random_state = 2
    np.random.seed(random_state)

    parallel = "seq" #sequential
    dbformat = "csv"
    spot_setup = spotpy_setup(new_meas)

    #See spotpy.dream.py
    rep = 5000
    nChains = 5 #default is 5
    convergence_limit = 1.2 #default is 1.2

    nCr = 3 #default is 3
    runs_after_convergence = 20000 #default is 100

    sampler = spotpy.algorithms.dream(spot_setup, parallel=parallel, dbname='HydrologyDREAM_5', dbformat=dbformat, db_precision = np.float16, save_sim = True, sim_timeout=timeout)
    r_hat = sampler.sample(repetitions = rep,  nChains = nChains, convergence_limit = convergence_limit, runs_after_convergence = runs_after_convergence)
    

    return r_hat
    # dbname: str
    #* Name of the database where parameter, objectivefunction value and simulation results will be saved.

    #   dbformat: str
    #   * ram: fast suited for short sampling time. no file will be created and results are saved in an array.
    #   * csv: A csv file will be created, which you can import afterwards.

    #   parallel: str
    #   * seq: Sequentiel sampling (default): Normal iterations on one core of your cpu.
    #   * mpi: Message Passing Interface: Parallel computing on cluster pcs (recommended for unix os).

    #   save_sim: boolean #default is True
    #   * True:  Simulation results will be saved
    #   * False: Simulation results will not be saved

if __name__ == "main":
    dream_run_5()