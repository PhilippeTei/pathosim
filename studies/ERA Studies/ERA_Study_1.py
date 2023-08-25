import numpy as np
import behaviour as bh
import pathosim as inf 
import matplotlib.pyplot as plt
import pandas as pd
import os
import sqlite3

#create dataframe
results_df = pd.DataFrame(columns=['Simulation_Setup', 'Mean_Detection', 'Mean_Costs', 'Mean_Runs'])

#define output directory
output_directory = "c:\\Outputs\\Simulations\\Results\\Replication_Main\\Study_1"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Base parameters for all simulations
base_pars = dict(
    pop_size = 5284,
    pop_type = 'random',
    n_days = 100,
    start_day = '2020-01-30'
)

# Random surveillance setup with placeholders
rand_pars = base_pars.copy()
rand_pars['enable_surveillance'] = True
rand_pars['surveillance_test_percent'] = None  # Placeholder
rand_pars['enable_random_testing'] = True
rand_pars['random_test_frequency'] = None  # Placeholder

# Define Iteration Lists:   
random_test_sample_values = [0.0002, 0.002, 0.02]
random_test_frequency_values = [1, 7, 14, 28]

def generate_label(pars):
    random_test_sample_size = pars.get('surveillance_test_percent', None)
    random_test_frequency = pars.get('random_test_frequency', None)

    # Generate a label for the simulation setup based on the parameters
    label = f"Random_Test_Sample_Size_{random_test_sample_size}_Random_Test_Frequency_{random_test_frequency}"

    return label

def run_simulation(pars, loops=50):

    # Initialize a list to hold the results of each simulation
    results_list = []
    detection_list = []
    run_list = []
    cost_list = []
    infections_at_detection_list = []

    # Run 100 simulations
    for i in range(loops):
        Coronavirus = inf.SARS_COV_2(4)
        Coronavirus.configure_transmission(beta=0.011457989455949544)
        sim = inf.Sim(pars=pars, pathogens=[Coronavirus], verbose=1, rand_seed=np.random.randint(0, 2**32 - 1, dtype=np.int64))
        sim.run()
        results_list.append(sim.results[0]['cum_infections'].values)
        detection_time = sim.confirmed_detection_time
        run_list.append(sim.total_runs)
        cost_list.append(sim.total_costs)
        if detection_time is not None:
            detection_list.append(detection_time)
            infections_at_detection = sim.results[0]['cum_infections'].values[int(detection_time) - 1]
            infections_at_detection_list.append(infections_at_detection)
            print(f'Simulation {i+1} complete with detection on day {detection_time}')
        else:
            print(f'Simulation {i+1} complete with no detection.')

    # Convert the list to a numpy array for easier calculations
    results_array = np.array(results_list)
    detection_array = np.array(detection_list)
    run_array = np.array(run_list)
    cost_array = np.array(cost_list)
    infections_at_detection_array = np.array(infections_at_detection_list)

    #initialise mean detection, mean runs, and mean costs
    mean_detection = np.nan
    mean_runs = np.nan
    mean_costs = np.nan
    mean_infections_at_detection = np.nan

    # only calculate mean if there are detections
    if len(detection_array) > 0:
        mean_detection = np.nanmean(detection_array)
        mean_runs = np.nanmean(run_array)
        mean_costs = np.nanmean(cost_array)
        mean_infections_at_detection = np.nanmean(infections_at_detection_array)

    # Function to get confidence intervals from samples
    def get_confidence_intervals(samples, upper=95):
        intervals = np.zeros((3, samples.shape[1]))
        intervals[0, :] = np.percentile(samples, (100-upper)/2, axis=0)
        intervals[1, :] = np.percentile(samples, 50, axis=0)
        intervals[2, :] = np.percentile(samples, upper + (100-upper)/2, axis=0)
        return intervals

    # Get the confidence intervals
    intervals = get_confidence_intervals(results_array)

    #generate label
    label = generate_label(pars)

    # Extracting the lower bound, median, and upper bound
    def plot_confidence_intervals(intervals, days):
        # Extracting the lower bound, median, and upper bound
        lower_bound = intervals[0, :]
        median = intervals[1, :]
        upper_bound = intervals[2, :]
        
        plt.figure(figsize=(10, 6))
        
        # Plot the median
        plt.plot(days, median, label='Median', color='blue')
        
        # Fill the area between the lower and upper confidence intervals
        plt.fill_between(days, lower_bound, upper_bound, color='blue', alpha=0.2, label='95% CI')
        
        # Plot the mean detection day if there was a detection
        if len(detection_array) > 0:
            plt.axvline(mean_detection, color='red', linestyle='--', label=f'Mean Detection Day: {int(round(mean_detection))}')
            cost_text = f"Mean Cost of Program: £{mean_costs:,.2f}"
            sims_text = f"Number of Simulations: {len(results_array)}"
            # Choose a suitable position for the annotations at the top left
            plt.annotate(cost_text, xy=(0.025, 0.80), xycoords='axes fraction', fontsize=10, color='green')
            plt.annotate(sims_text, xy=(0.025, 0.75), xycoords='axes fraction', fontsize=10, color='green')

        plt.xlabel('Day')
        plt.ylabel('Number of Infections')
        plt.title('Number of Infections with 95% Confidence Intervals')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plot_filename = os.path.join(output_directory, f"ThreatNet_{label}.png")
        plt.savefig(plot_filename, dpi=300)
        plt.close()

        return plot_filename

    #calculate mean results
    mean_results = np.mean(results_array, axis=0)

    # After generating the intervals, plot them:
    days = np.arange(0, pars['n_days'] + 1)
    intervals = get_confidence_intervals(results_array)
    plot_confidence_intervals(intervals, days)

    #Save to CSV
    df = pd.DataFrame({
        'Day': sim.t,
        'Mean_Infections': mean_results
    })

    filename = os.path.join(output_directory, f"Mean_Results_{label}.csv")
    df.to_csv(filename, index=False)

    #print mean detection day, cost of program, and number of simulations if there was a detection
    if len(detection_array) > 0:
        print(f'Mean Detection Day: {mean_detection:.2f}')
        print(f'Mean Number of Sequencing Runs: {mean_runs:.2f}')
        print(f'Mean Cost of Program: £{mean_costs:,.2f}')
        print(f'Number of Simulations: {len(results_array)}')
        print(f'Mean Infections at Detection: {mean_infections_at_detection:.2f}')

    return {
        'Simulation_Setup': generate_label(pars),
        'Mean_Detection': mean_detection,
        'Mean_Costs': mean_costs,
        'Mean_Runs': mean_runs,
        'Mean_Infections_at_Detection': mean_infections_at_detection
    }

# Run the simulations
for random_test_sample_size in random_test_sample_values:
    for random_test_frequency in random_test_frequency_values:
        rand_pars['surveillance_test_percent'] = random_test_sample_size
        rand_pars['random_test_frequency'] = random_test_frequency
        result = run_simulation(rand_pars, 50)
        results_df = pd.concat([results_df, pd.DataFrame([result], columns=results_df.columns)], ignore_index=True)

#export to csv
results_csv_path = os.path.join(output_directory, "cumulative_results.csv")
results_df.to_csv(results_csv_path, index=False)

#TOTAL: 12 simulations; 100 runs each; 1200 runs total