import numpy as np
import behaviour as bh
import pathosim as inf 
import matplotlib.pyplot as plt
import pandas as pd
import os

#import data
calib_data = pd.read_csv(f'C:\Outputs\Calibration Data\ThreatNet_Test_Data - Sheet2.csv')
calib_days = calib_data['day'].values
calib_infections = calib_data['infections_column'].values

pars = dict(
    pop_size = 5284,
    pop_type = 'random',
    n_days = 100,
    start_day = '2020-01-30',
    enable_surveillance = False,
    enable_syndromic_testing = False,
    syndromic_pooling_enabled = False,
    syndromic_pool_size = 5
)

# Initialize a list to hold the results of each simulation
results_list = []
detection_list = []
run_list = []
cost_list = []

# Run 10 simulations
for i in range(100):
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
        print(f'Simulation {i+1} complete with detection on day {detection_time}')
    else:
        print(f'Simulation {i+1} complete with no detection.')

# Convert the list to a numpy array for easier calculations
results_array = np.array(results_list)
detection_array = np.array(detection_list)
run_array = np.array(run_list)
cost_array = np.array(cost_list)

# only calculate mean if there are detections
if len(detection_array) > 0:
    mean_detection = np.nanmean(detection_array)
    mean_runs = np.nanmean(run_array)
    mean_costs = np.nanmean(cost_array)

# Function to get confidence intervals from samples
def get_confidence_intervals(samples, upper=95):
    intervals = np.zeros((3, samples.shape[1]))
    intervals[0, :] = np.percentile(samples, (100-upper)/2, axis=0)
    intervals[1, :] = np.percentile(samples, 50, axis=0)
    intervals[2, :] = np.percentile(samples, upper + (100-upper)/2, axis=0)
    return intervals

# Get the confidence intervals
intervals = get_confidence_intervals(results_array)

def plot_confidence_intervals(intervals, days):
    # Extracting the lower bound, median, and upper bound
    lower_bound = intervals[0, :]
    median = intervals[1, :]
    upper_bound = intervals[2, :]
    
    plt.figure(figsize=(10, 6))
    
    # plot the reference line
    if calib_days is not None and calib_infections is not None:
        plt.plot(calib_days, calib_infections, 'o-', color='green', label='Threat Net Infection Progression Estimate Data')

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
    plot_filename = 'ThreatNet_Replicated {}.png'
    plt.savefig(plot_filename, dpi=300)

#calculate mean results
mean_results = np.mean(results_array, axis=0)

#Save to CSV
df = pd.DataFrame({
    'Day': sim.t,
    'Mean_Infections': mean_results
})

filename = 'mean_results.csv'
if os.path.exists(filename):
    os.remove(filename)

df.to_csv('mean_results.csv', index=False)

#print mean detection day, cost of program, and number of simulations if there was a detection
if len(detection_array) > 0:
    print(f'Mean Detection Day: {mean_detection:.2f}')
    print(f'Mean Number of Sequencing Runs: {mean_runs:.2f}')
    print(f'Mean Cost of Program: £{mean_costs:,.2f}')
    print(f'Number of Simulations: {len(results_array)}')


# Plot the results
days = np.arange(0, pars['n_days'] + 1)
intervals = get_confidence_intervals(results_array, upper=95)
plot_confidence_intervals(intervals, days)