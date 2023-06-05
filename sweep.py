import infection as inf 
from infection import covid_testing as ct
import behaviour as bh
from behaviour import BehaviourModel
import copy
import numpy as np
import os


def parameter_sweep(nominal_dict, sweep_dict, file_loc, n_runs=3, parallel=False, build_pop=False):
    '''
    Function that conducts parameter sweeps and runs many simulations.

    Args: 
        nominal_dict   (dict) : Contains a specified set of nominal parameters, which should be defaulted to when other parameters are being swept
        sweep_dict     (dict) : Contains the parameter values to sweep over. Values should be in the format: [initial_value, final_value, num_steps]
        two_sweep      (list) : List of 2-element sets. Elements are parameter names, and each set indicates a duo of parameters to sweep in conjunction.
        n_runs         (int)  : Size of multisim for each parameter
        file_loc       (str)  : Location to save all the msims
    '''
    # Sweep over simulation parameters first
    for param_name in sweep_dict['sim_sweep'].keys(): 
        test_pars = nominal_dict['test_pars'].copy()
        sim_pars = nominal_dict['sim_pars'].copy()

        sweep_vals = sweep_dict['sim_sweep'][param_name]  # Get the values to sweep
        for val in sweep_vals:  # Sweep over the values
            default = sim_pars['pars'][param_name]
            sim_pars['pars'][param_name] = val

            testobjs = build_tests(test_pars)
            sim = build_sim(sim_pars, testobjs, build_pop)
            save_name = f'simsweep_{param_name}_{str(val)}_{str(n_runs)}'
            run_multisim(sim, n_runs, file_loc, save_name, parallel)

            sim_pars['pars'][param_name] = default  # Reset


    # Then sweep over testing parameters for each test object
    for testname in sweep_dict['test_sweep'].keys(): 
        for param_name in sweep_dict['test_sweep'][testname]:
            test_pars = nominal_dict['test_pars'].copy()  # Ensures that parameters not being swept are returned to default
            sim_pars = nominal_dict['sim_pars'].copy()
            
            sweep_vals = sweep_dict['test_sweep'][testname][param_name]
            for val in sweep_vals: 
                default = test_pars[testname][param_name]
                test_pars[testname][param_name] = val

                testobjs = build_tests(test_pars)
                sim = build_sim(sim_pars, testobjs, build_pop)
                save_name = f'{testname}sweep_{param_name}_{str(val)}_{str(n_runs)}'
                run_multisim(sim, n_runs, file_loc, save_name, parallel)

                test_pars[testname][param_name] = default  # Reset
            
    return 

def build_tests(test_pars):
    options = ['PCR_disc', 'RAT_disc', 'RAT_surv']
    testobjs = []
    for testname, pars in test_pars.items(): 
        if not(testname in options): 
            raise RuntimeError(f'Test parameter dictionary specifies a non-existent test object. Options are {options}')
        
        # These are the only test objects available at the moment
        if testname == 'RAT_disc': 
            testobjs.append(ct.RAT_disc(**pars))  # Unpack dictionary parameters as keyword arguments. Keys need to match argument name. 
        elif testname == 'PCR_disc': 
            testobjs.append(ct.PCR_disc(**pars))
        elif testname == 'RAT_surv':
            testobjs.append(ct.RAT_surv(**pars))

    return testobjs

def build_sim(sim_pars, testobjs, build_pop):
    '''
    Build a simulation object. Works faster if the population dictionary is directly supplied through sim_pars. Otherwise, a new population will need to be initialized each time. However, this may be inevitable
    if you want to sweep population parameters (such as population size).
    Args: 
        sim_pars   (dict)  : Dictionary of all simulation parameters, except for testingn objects 
        testobjs   (list)  : List of defined test objects
    '''
    sim_pars['pars']['testing'] = testobjs
    if build_pop:  # This is currently slow (probably)
        pop = build_pop(sim_pars['people'])  # Population parameters supplied here
        sim_pars['people'] = pop.popdict
        sim_pars['pars']['pop'] = pop
    return inf.Sim(**sim_pars)

def build_pop(pop_pars):
    return BehaviourModel(pars=pop_pars)

def run_multisim(sim, n_runs, file_loc, save_name, parallel): 
    '''
    Given a single sim, run it multiple times
    '''
    msim = inf.MultiSim(sim, parallel=parallel)
    msim.run(n_runs)
    save_path = os.path.join(file_loc, save_name) + '.msim'
    msim.save(save_path, keep_people=True) 
    return



if __name__ == '__main__':
    # Set population and simulation defaults
    pop_type = 'behaviour_module'
    pop_size = 100000
    n_days = 200
    country_location='usa'
    state_location='New_York'
    county = 'Kings_County'

    # Generate synthetic population 
    popPars=dict(
        n=pop_size,
        with_school_types=True,
        school_mixing_type='age_and_class_clustered',
        country_location=country_location,
        state_location=state_location,
        location=county
    )
    pop = BehaviourModel(pars=popPars)

    # Construct the nominal parameter dictionary
    # Primary 
    pars = dict(pop_size=100000,  # Investigate: sims don't work with different population sizes
            pop_type='behaviour_module',
            pop=pop,
            n_days=200,
            bkg_ILI=0.04,
            symp_prev_multiplier=1)

    RATpars = dict(sensitivity=0.70,
            specificity=0.99,
            capacity=5000)
    PCRpars = dict(sensitivity=0.85, 
            specificity=1.0, 
            capacity=400, 
            RAT_ind=1)  # Index in which RAT_disc object will be placed in test_pars dictinoary

    # Secondary 
    test_pars = dict(PCR_disc = PCRpars,
                    RAT_disc = RATpars)

    sim_pars = dict(pars=pars, people=pop.popdict)

    # Tertiary
    nominal_dict = dict(sim_pars=sim_pars,
                        test_pars=test_pars)

    # Construct the parameter sweeping dictionary 
    # Primary 
    PCRsweep = dict(sensitivity=np.arange(0.7, 0.91, 0.04), 
                specificity=np.arange(0.8, 1.01, 0.04), 
                capacity=np.arange(250, 550, 50))

    RATsweep = dict(sensitivity=np.arange(0.6, 0.81, 0.04),
                specificity=np.arange(0.8, 1.01, 0.04),
                capacity=np.arange(4100, 5900, 300))

    # Secondary 
    test_sweep = dict(PCR_disc = PCRsweep,
                    RAT_disc = RATsweep)

    sim_sweep = dict(bkg_ILI=np.arange(0.03, 0.06, 0.005),
            symp_prev_multiplier=np.arange(0.6, 1.2, 0.1))

    # Tertiary
    sweep_dict = dict(sim_sweep=sim_sweep,
                    test_sweep=test_sweep)

    loc = '/Users/Ritchie/Desktop/SeroTracker/Code/Stage 7/parameter_sweep_msims_v2'
    parameter_sweep(nominal_dict, sweep_dict, loc)
    

