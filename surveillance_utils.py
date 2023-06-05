# Utility functions imported into surveillance_systems.py 

# Some of these imports are unnecessary 
import numpy as np
import pandas as pd

from . import utils as cvu
from .background_ILI import * 

# def retest_eligibility(surv):  # Eventually it may be wise to implement this generally
#     return

###############################
### Preprocessing functions ###
###############################

def process_UK_PCR_capacity(fname='/Users/Ritchie/Desktop/SeroTracker/Data/UK_PCR.csv', pop=67.22e6, window=10, capacity_per=1e5):
    '''
    Gets daily test capacity per 100k people, using United Kingdom data. Smooths the resulting test capacity curve. This function works specifically for United Kingdom 
    data obtained here: https://coronavirus.data.gov.uk/details/testing

    Args: 
        fname           (str): csv file containing unprocessed test capacity data
        pop             (int): UK population 
        window          (int): Smoothing window
        capacity_per()  (int): Desired rate of tests (e.g. tests per 100k)

    Returns: 
        smooth_cap  (np.array): Smoothed test capacity
        dates       (np.array): Dates of each test capacity
        
    '''
    data = pd.read_csv(fname)
    dates = data['date'].iloc[::-1].values  # Corresponds to capacity by index 
    capacity = data['plannedPCRCapacityByPublishDate'].iloc[::-1].values
    capacity *= (capacity_per / pop)  # Number of tests per 100k people 
    smooth_cap = np.convolve(capacity, np.ones(window), 'valid') / window  # Smooth the capacity curve
    return smooth_cap, dates


# def process_England_antigen_capacity():
# TODO: Implement this function, data located here, "Rapid lateral flow testes reported": https://coronavirus.data.gov.uk/details/testing?areaType=nation&areaName=England



###############################
###### Testing functions ######
###############################

def PCR_test(sim, i, se, sp): 
    '''
    Administer PCR test to individual i in the population. Assume people have COVID-19 only if they are exposed (i.e. any time between infection to recovery)

    Args: 
        sim     (sim): Simulation object
        i       (int): Index of individual 
        se      (float): Test sensitivity 
        sp      (float): Test specificity

    Returns: 
        result  (bool): Truth value of test
    '''
    if sim.people.exposed[i]: 
        return np.random.uniform(0,1) < se
    else:
        return np.random.uniform(0,1) > sp 


def antigen_test(sim, i, se, sp): 
    '''
    Administer antigen test to individual i in the population. Assume people have COVID-19 only if they are exposed (i.e. any time between infection to recovery)

    Args: 
        sim     (sim): Simulation object
        i       (int): Index of individual 
        se      (float): Test sensitivity 
        sp      (float): Test specificity

    Returns: 
        result  (bool): Truth value of test
    '''
    if sim.people.exposed[i]: 
        return np.random.uniform(0,1) < se
    else:
        return np.random.uniform(0,1) > sp 


def calculate_retest_delay(label, sim, i, pos_delay, neg_delay): 
    '''
    Calculate how long until individual i in the population can retest again. Only use this if the individual just tested on the current day. 

    Args: 
        surveillance    (surveillance): Surveillance object
        sim             (sim): Simulation object
        i               (int): Index of individual 
        pos_delay       (int): Number of days before they can retest, if they tested positive
        neg_delay       (int): Number of days before they can retest, if they tested negative
    Returns: 
        retest_delay    (int): Number of days before they can retest
    '''
    test_result = sim.people[label][i]
    if test_result == True: 
        return pos_delay
    else: 
        return neg_delay



###############################
###### Contact functions ######
###############################

def get_contact_counts(sim): 
    '''
    Get number of contacts for each individual initialized in the population, by iterating through the edge list. Note that an individual i may appear in both 
    p1 and p2 lists. Furthermore, the contacts of individual i in p1 are unique from the contacts of individual i in p2. Thus we cannot get all conatcts of individual i 
    just by iterating through one of p1 or p2. 
    
    Args: 
        sim             (sim): Initialized simulation object

    Returns: 
        contact_count   (np.array): Array of net contact counts for each age 
    '''
    contact_count_dict = dict()
    for key in sim.people.contacts.keys(): 
        contact_count_dict[key] = np.full(sim.pars['pop_size'], 0)

        for p1 in sim.people.contacts[key]['p1']: 
            contact_count_dict[key][p1] += 1  # Count number of contacts

        for p2 in sim.people.contacts[key]['p2']: 
            contact_count_dict[key][p2] += 1

    contact_count = np.full(sim.pars['pop_size'], 0)
    for key in contact_count_dict.keys(): 
        print(key)
        contact_count += contact_count_dict[key]
        
    return contact_count



###############################
#### Allocation functions #####
###############################

#TODO: Move to surveillance_allocation.py

def distribute_wearables(sim, rel_age_distr): 
    '''
    Distribute wearables among population, using fixed age distribution. 

    Args: 
        w_surv  (surveillance): Wearable surveillance object
        sim     (sim): Simulation object

    Returns: 
        w_inds  (list): List of 
    '''
    ages = sim.people.age
    # Get boolean indices of age groups 
    A_bool = np.bitwise_and(ages >= 18, ages <= 34)
    B_bool = np.bitwise_and(ages >= 35, ages <= 49)
    C_bool = np.bitwise_and(ages >= 50, ages <= 64)
    D_bool = np.bitwise_and(ages >= 65, ages <= 74)
    E_bool = ages >= 75

    # Get the actual indices of age groups 
    A_ind = cvu.true(A_bool)
    B_ind = cvu.true(B_bool)
    C_ind = cvu.true(C_bool)
    D_ind = cvu.true(D_bool)
    E_ind = cvu.true(E_bool)

    # Get indices of who gets a wearable, distributed by age
    A_w = np.random.choice(A_ind, int(len(A_ind) * rel_age_distr['A, 18-34']), replace=False)  
    B_w = np.random.choice(B_ind, int(len(B_ind) * rel_age_distr['B, 35-49']), replace=False)
    C_w = np.random.choice(C_ind, int(len(C_ind) * rel_age_distr['C, 50-64']), replace=False)
    D_w = np.random.choice(D_ind, int(len(D_ind) * rel_age_distr['D, 65-74']), replace=False)
    E_w = np.random.choice(E_ind, int(len(E_ind) * rel_age_distr['E, 75+']), replace=False)

    w_inds = np.concatenate((A_w, B_w, C_w, D_w, E_w), axis=None)
    num_w = len(w_inds)

    return w_inds, num_w



###############################
##### Analysis functions ######
###############################

# General analysis functions that do not directly relate to test allocation or mechanism

def log_symptomatic(surv, sim, meets_criteria): 
    '''
    Get all people who are symptomatic across time. 
    '''
    if sim.t == 0: 
        setattr(surv, 'symptomatic_log', [])

    surv.symptomatic_log.append(len(meets_criteria['sx']))


#TODO: Move everything below to surveillance_allocation.py

def relative_w_age_distr(w_surveillance, sim): 
    '''
    Get relative age distribution of an initialized wearable surveillance system. Does not work for other surveillance systems. 
    
    Args: 
        surveillance (surveillance): Wearable surveillance object
        sim          (sim): Initialized simulation object
    '''
    # Get ages
    ages = []
    for ind in w_surveillance.ind: 
        age = sim.people.age[ind]
        ages.append(age)
    w_ages = np.array(ages)

    all_ages = sim.people.age

    # Compute relative wearable age distribution
    rel_Z = sum(w_ages < 18) / sum(all_ages < 18)
    rel_A = sum(np.bitwise_and(w_ages >= 18, w_ages <= 34)) / sum(np.bitwise_and(all_ages >= 18, all_ages <= 34))  # (n people in group with wearables) / (m people in group)
    rel_B = sum(np.bitwise_and(w_ages >= 35, w_ages <= 49)) / sum(np.bitwise_and(all_ages >= 35, all_ages <= 49))
    rel_C = sum(np.bitwise_and(w_ages >= 50, w_ages <= 64)) / sum(np.bitwise_and(all_ages >= 50, all_ages <= 64))
    rel_D = sum(np.bitwise_and(w_ages >= 65, w_ages <= 74)) / sum(np.bitwise_and(all_ages >= 65, all_ages <= 74))
    rel_E = sum(w_ages >= 75) / sum(all_ages >=75)

    x = np.array(['Z', 'A', 'B', 'C', 'D', 'E'])
    y = np.array([rel_Z, rel_A, rel_B, rel_C, rel_D, rel_E])
    
    return x, y


def relative_infect_by_age(surveillance, sim): 
    '''
    Get proportions of each age group exposed at each time step. Needs to be run repeatedly during the simulation. 
    Args: 
        surveillance (surveillance): Surveillance object 
        sim          (sim): Simulation object
    '''
    keys = ['rel_Z', 'rel_A', 'rel_B', 'rel_C', 'rel_D', 'rel_E']
    if sim.t == 0: 
        setattr(surveillance, 'rel_infect', dict()) # Initialize tracker dictionary
        for k in keys: 
            surveillance.rel_infect[k] = []
    elif sim.t > 0: 
        ages = sim.people.age

        # Get all individuals exposed at current time step
        exposed = sim.people.exposed  # Boolean array

        # Compute absolute proportions of people exposed
        abs_Z = sum(np.bitwise_and(exposed, ages < 18))
        abs_A = sum(np.bitwise_and(exposed, np.bitwise_and(ages >= 18, ages <= 34)))  # Exposed and age belongs in range [18, 34]
        abs_B = sum(np.bitwise_and(exposed, np.bitwise_and(ages >= 35, ages <= 49)))
        abs_C = sum(np.bitwise_and(exposed, np.bitwise_and(ages >= 50, ages <= 64)))
        abs_D = sum(np.bitwise_and(exposed, np.bitwise_and(ages >= 65, ages <= 74)))
        abs_E = sum(np.bitwise_and(exposed, ages >= 75))

        # Compute relative proportions of people exposed
        rel_Z = abs_Z / sum(ages < 18)
        rel_A = abs_A / sum(np.bitwise_and(ages >= 18, ages <= 34))
        rel_B = abs_B / sum(np.bitwise_and(ages >= 35, ages <= 49))
        rel_C = abs_C / sum(np.bitwise_and(ages >= 50, ages <= 64))
        rel_D = abs_D / sum(np.bitwise_and(ages >= 65, ages <= 74))
        rel_E = abs_E / sum(ages >= 75)

        values = [rel_Z, rel_A, rel_B, rel_C, rel_D, rel_E]

        for k, v in zip(keys, values): 
            surveillance.rel_infect[k].append(v)

    return


def get_users_by_age(surveillance, sim): 
    '''
    Get the number of wearable users by age. 
    '''
    user_ages = sim.people.age[surveillance.ind]

    Z = sum(user_ages < 18)
    A = sum(np.bitwise_and(user_ages >= 18, user_ages <= 34))
    B = sum(np.bitwise_and(user_ages >= 35, user_ages <= 49))
    C = sum(np.bitwise_and(user_ages >= 50, user_ages <= 64))
    D = sum(np.bitwise_and(user_ages >= 65, user_ages <= 74))
    E = sum(user_ages >= 75)

    age_dict = dict()
    keys = ['Z', 'A', 'B', 'C', 'D', 'E']   
    values = [Z, A, B, C, D, E]

    for k, v in zip(keys, values): 
        age_dict[k] = v

    return age_dict


def abs_wearable_signal_by_age(surveillance, ages, sim): 
    '''
    Get the wearable surveillance signal stratified by age. (Consider implementing a similar function for the other surveillance systems)
    Args: 
        surveillance (surveillance): Surveillance object 
        sim          (sim): Simulation object
    '''
    keys = ['Z', 'A', 'B', 'C', 'D', 'E']
    if sim.t == 0:
        setattr(surveillance, 'abs_w_signal', dict()) # Initialize tracker dictionary
        for k in keys: 
            surveillance.abs_w_signal[k] = []

    elif sim.t > 0: 
        warned_ages = ages  # All individuals who were warned, represented by their age
        Z = sum(warned_ages < 18)
        A = sum(np.bitwise_and(warned_ages >= 18, warned_ages <= 34))
        B = sum(np.bitwise_and(warned_ages >= 35, warned_ages <= 49))
        C = sum(np.bitwise_and(warned_ages >= 50, warned_ages <= 64))
        D = sum(np.bitwise_and(warned_ages >= 65, warned_ages <= 74))
        E = sum(warned_ages >= 75)

        values = [Z, A, B, C, D, E]

        for k, v in zip(keys, values): 
            surveillance.abs_w_signal[k].append(v)
            
    return


def rel_wearable_signal_by_age(surveillance, ages, sim):
    '''
    Returns relative daily distribution of wearable warnings, stratified by age. 

    '''
    abs_wearable_signal_by_age(surveillance, ages, sim)

    keys = ['Z', 'A', 'B', 'C', 'D', 'E']
    if sim.t == 0: 
        setattr(surveillance, 'rel_w_signal', dict())
        setattr(surveillance, 'users_by_age', get_users_by_age(surveillance, sim))
        for k in keys: 
            surveillance.rel_w_signal[k] = []

    elif sim.t > 0: 
        for k in keys: 
            if surveillance.users_by_age[k] == 0: 
                surveillance.rel_w_signal[k].append(0)
            else: 
                surveillance.rel_w_signal[k].append(surveillance.abs_w_signal[k][-1] / surveillance.users_by_age[k])
    return











