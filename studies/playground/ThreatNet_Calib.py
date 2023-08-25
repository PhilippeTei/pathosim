import optuna
import numpy as np
import behaviour as bh
import pathosim as inf 
import pandas as pd
import os

# Load your empirical data
empirical_df = pd.read_csv("C:\Simdata\ThreatNet_Test_Data.csv")
empirical_data = empirical_df['infections_column'].values

def compute_error(beta, seed, pop_size):
    pars = dict(
        pop_size=pop_size,
        pop_type='random',
        n_days=100,
        start_day='2020-01-30',
        enable_surveillance=False,
        enable_syndromic_testing=False,
    )

    results_list = []
    for i in range(10):
        Coronavirus = inf.SARS_COV_2(seed)
        Coronavirus.configure_transmission(beta=beta)
        sim = inf.Sim(pars=pars, pathogens=[Coronavirus], verbose=0, rand_seed=np.random.randint(0, 2**32 - 1, dtype=np.int64))
        sim.run()
        results_list.append(sim.results[0]['cum_infections'].values)

    results_array = np.array(results_list)
    mean_results = np.mean(results_array, axis=0)
    
    # Compute the mean squared error (MSE) between the empirical and simulated data
    error = np.mean((mean_results - empirical_data) ** 2)
    return error

# Define the objective function
def objective(trial):
   # Sample a value for each parameter
    beta = trial.suggest_float('beta', 0.008, 0.012)
    seed = trial.suggest_int('seed', 1, 10)
    pop_size = trial.suggest_int('pop_size', 4500, 5500)
    
    return compute_error(beta, seed, pop_size)

#periodically save the study object


# Create a study object and specify the direction is 'minimize'.
#optuna.storage.get_storage("sqlite:///ThreatNet_Calib.db").delete_study(study_name="ThreatNet_Calib")
study = optuna.create_study(storage="sqlite:///ThreatNet_Calib.db", study_name="ThreatNet_Calib", direction='minimize', load_if_exists=True)

# Optimize the study, the objective function is passed in as the first argument.
study.optimize(objective, n_trials=100)

# Results
print('Number of finished trials: ', len(study.trials))
print('Best trial:')
trial = study.best_trial

print(f"Value of beta: {trial.params['beta']}")
print(f"Value of seed: {trial.params['seed']}")
print(f"Value of pop_size: {trial.params['pop_size']}")
print('Minimum Mean Squared Error: ', trial.value)