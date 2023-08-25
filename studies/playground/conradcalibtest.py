import pathosim as inf
import behaviour as bh
import pandas as pd
import sciris as sc
import os
import multiprocessing

ukdata = pd.read_csv(r"C:\ukearlydata - Sheet1.csv")

pars = sc.objdict(
    pop_size = 1000,
    pop_type = 'behaviour_module',
    start_day = '2020-01-30',
    n_days = 50,
    verbose = 0,
)

MyNewPathogen = inf.SARS_COV_2(100)   
sim = inf.Sim(pars=pars, pathogens = [MyNewPathogen], verbose = 1, datafile = ukdata)

if __name__ == '__main__':
    # Parameters to calibrate -- format is best, low, high
    calib_pars = dict(
        pop_size      = [pars.pop_size, 1000, 10000]
)
    # Run the calibration
    if os.path.exists('covasim_calibration.db'):
        os.remove('covasim_calibration.db')
                  
    calib = inf.calibrate(calib_pars=calib_pars, total_trials=50, keep_db = True)

