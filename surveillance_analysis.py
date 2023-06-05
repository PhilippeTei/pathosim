import numpy as np 
import copy 
from scipy.stats import pearsonr 
from sklearn.metrics import mean_squared_error as MSE

# Compute outcome measures for a multisim input
def get_timeliness(msim, i):
    '''
    Get timeliness of a testobj, indexed by i (sim.pars['testing'][i])
    '''
    return dot_crosscorr(compute_msim_incidence(msim), compute_msim_pos(msim)[i])[0]

def get_sensitivity(msim, i):
    total_true_inc = np.sum(compute_msim_incidence(msim))
    total_captured = compute_msim_pos(msim)[i]
    return round(np.sum(total_captured) / total_true_inc, 3)

def get_correlation(msim, i, method='dot'): 
    true_inc = compute_msim_incidence(msim)
    captured = compute_msim_pos(msim)[i]
    return round(adjusted_pearsonr(true_inc, captured, method=method), 3)

def get_underascertainment(msim, i): 
    '''
    Inverse operation of get_sensitivity
    '''
    total_true_inc = np.sum(compute_msim_incidence(msim))
    total_captured = compute_msim_pos(msim)[i]
    return round(total_true_inc / np.sum(total_captured), 3)

# Compute outcome measures for a single sim input
def get_timeliness_single(sim): 
    return dot_crosscorr(sim.results['new_infections'].values, sim.pars['testing'][0].pos_history)[0]

def get_sensitivity_single(sim):
    return round(np.sum(sim.pars['testing'][0].pos_history) / np.sum(sim.results['new_infections'].values), 3)

def get_correlation_single(sim): 
    return round(adjusted_pearsonr(sim.results['new_infections'].values, sim.pars['testing'][0].pos_history, method='dot'), 3)


# Timeliness analysis
def dot_crosscorr(x, y, n_days=200, max_delay=20):
    '''
    Function to compute time lagged cross correlation  using dot products. Slide one signal and then 
    compute dot product of two signals, maximum of this repeated computation corresponds to the 
    estimated lag. 
    
    Args: 
        x         (np.array): Ground truth array 
        y         (np.array): Surveillance array 
        max_delay (int)     : Expected to be an upper limit on maximum possible signal
        
    Returns: 
        index   (int): Estimated lag of the surveillance signal
        dot_lst (arr): List of all computed dot products for visualization purposes 
    
    '''
    dot_lst = [] 
    delay_lst = np.arange(-max_delay, max_delay+1)

    start = max_delay
    end = n_days - max_delay

    for d in delay_lst:
        dot_lst.append(np.dot(x[start:end], y[start+d:end+d]))  # Think of it as a sliding window
        
    return np.argmax(dot_lst) - max_delay, dot_lst  # Delay corresponding to largest dot product, max_delay corrects the index


def pr_crosscorr(x, y, n_days=200, max_delay=20):
    '''
    Function to compute time lagged cross correlation using Pearson correlation. Slide one signal 
    and then compute dot product of two signals, maximum of this repeated computation corresponds 
    to the estimated lag. 
    
    Args: 
        x       (np.array): Ground truth array 
        y       (np.array): Surveillance array 
        start   (int): Start index
        end     (int): End index 
        delay   (int): Max delay 
        
    Returns: 
        index   (int): Estimated lag of the surveillance signal
        pr_lst (arr): List of all computed Pearson correlations for visualization purposes 
    
    '''
    pr_lst = []  # Hold the Pearson correlations
    delay_lst = np.arange(-max_delay, max_delay+1)

    start = max_delay
    end = n_days - max_delay
    
    for d in delay_lst:
        pr_lst.append(pearsonr(x[start:end], y[start+d:end+d])[0]) 
        
    return np.argmax(pr_lst) -  max_delay, pr_lst  # Delay which corresponds to largest Pearson correlation 


# Completeness analysis
def adjusted_pearsonr(x, y, n_days=200, max_delay=20, method='average'): 
    '''
    max_delay       (int) : Must be the same as the max_delay used for timeliness
    '''

    dot_delay = dot_crosscorr(x, y, n_days, max_delay)[0]
    pr_delay = pr_crosscorr(x, y, n_days, max_delay)[0]
    avg_delay = int((dot_delay + pr_delay) / 2)

    start = max_delay 
    end = n_days - max_delay

    if method == 'average': 
        return pearsonr(x[start:end], y[start+avg_delay:end+avg_delay])[0]
    elif method == 'dot': 
        return pearsonr(x[start:end], y[start+dot_delay:end+dot_delay])[0]
    elif method == 'pearson': 
        return pearsonr(x[start:end], y[start+pr_delay:end+pr_delay])[0]
    else: 
        raise RuntimeError('Method must be either average, dot, or pearson')


def adjusted_RMSE(x, y, n_days=200, max_delay=20, method='average'): 
    '''
    max_delay       (int) : Must be the same as the max_delay used for timeliness
    '''

    dot_delay = dot_crosscorr(x, y, n_days, max_delay)[0]
    pr_delay = pr_crosscorr(x, y, n_days, max_delay)[0]
    avg_delay = int((dot_delay + pr_delay) / 2)

    start = max_delay 
    end = n_days - max_delay

    if method == 'average': 
        return MSE(x[start:end], y[start+avg_delay:end+avg_delay], squared=False)
    elif method == 'dot': 
        return MSE(x[start:end], y[start+dot_delay:end+dot_delay], squared=False)
    elif method == 'pearson': 
        return MSE(x[start:end], y[start+pr_delay:end+pr_delay], squared=False)
    else: 
        raise RuntimeError('Method must be either average, dot, or pearson')


# Functions for generating tornado sensitivity plot


# Helper functions
def compute_msim_incidence(msim): 
    '''
    Return average incidence of infection of an msim.
    '''
    incidences = np.array([sim.results['new_infections'].values for sim in msim.sims])
    return np.mean(incidences, axis=0)

def compute_msim_pos(msim):
    '''
    Return a list containing the average number of positive cases for each testobj over time, given an msim. 
    '''
    pos_avg_lst = []
    for i in range(len(msim.sims[0].pars['testing'])): 
        pos_avg_lst.append(np.mean(np.array([sim.pars['testing'][i].pos_history for sim in msim.sims]), axis=0))
    return pos_avg_lst

def compute_msim_pipeline(msim):
    '''
    Return a list of dummy testobjs, each holding an average of pipeline (eligibility, seek, allocate) signals. 
    '''
    dummy_lst = []
    for i in range(len(msim.sims[0].pars['testing'])):  # Iterate through each type of testing object 
        testobj_lst = [sim.pars['testing'][i] for sim in msim.sims]  # Collect the particular testing object from each simulation 

        dummy = copy.deepcopy(msim.sims[0].pars['testing'][i])  # Create a dummy testobj to hold the averaged data

        cg = dict()
        sg = dict()
        ag = dict()

        for key in dummy.crit_groups.keys(): 
            cg[key] = np.mean(np.array([testobj.crit_groups[key] for testobj in testobj_lst]), axis=0)
        for key in dummy.seek_groups.keys(): 
            sg[key] = np.mean(np.array([testobj.seek_groups[key] for testobj in testobj_lst]), axis=0)
        for key in dummy.alloc_groups.keys(): 
            ag[key] = np.mean(np.array([testobj.alloc_groups[key] for testobj in testobj_lst]), axis=0)

        # Main signals
        dummy.crit_groups = cg
        dummy.seek_groups = sg
        dummy.alloc_groups = ag

        # Other signals
        dummy.tests_consumed = np.mean(np.array([testobj.tests_consumed for testobj in testobj_lst]), axis=0)
        
        dummy_lst.append(dummy)

    return dummy_lst


if __name__ == '__main__': 

    # Test cases 
    ground_truth =    np.array([-2, -1, 0, 1, 2, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7])
    test_arr =         np.array([0,  1, 2, 3, 4, 5, 6, 7, 6,  5,  4,  3,  2,  1,  0, -1])

    print("Timeliness")
    print("Inner product cross correlation: ", dot_crosscorr(ground_truth, test_arr, n_days=len(ground_truth), max_delay=4)[0])
    print("Pearson cross correlation: ", pr_crosscorr(ground_truth, test_arr, n_days=len(ground_truth), max_delay=4)[0])

    print()
    print("Completeness")
    print("Time-adjusted Pearson: ", adjusted_pearsonr(ground_truth, test_arr, n_days=len(ground_truth), max_delay=4))
    print("Time-adjusted RMSE: ", adjusted_RMSE(ground_truth, test_arr, n_days=len(ground_truth), max_delay=4))

