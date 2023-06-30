# Functions for defining test allocation
import numpy as np
import copy

from itertools import combinations
from functools import reduce
from pathosim import symptoms as symptoms
from pathosim.background_ILI import * 

#### Main ####
# random_allocation
# proportional_allocation
# cross_allocation 
# long_allocation 


#### Utility ####
# (For proportional_allocation)
# compute_members 
# compute_proportions
# compute_criteria_powset

# (For cross/long_allocation)
# get_targets 
# sample_targets
# calculate_test_days
# get_age_group_members

# Main call 
def allocate_tests(sim, testobj, seekers_dict=None, seekprobs=None, mode='proportional'): 
    if mode == 'proportional': 
        if isinstance(testobj.capacity, np.ndarray) or isinstance(testobj.capacity, list):
            capacity = int(testobj.capacity[sim.t] * sim.pars['pop_size'])
        else: 
            capacity = int(testobj.capacity * sim.pars['pop_size'])
        return proportional_allocation(capacity, seekers_dict, seekprobs, testobj.name, sim.t)

    elif mode == 'basic': 
        if isinstance(testobj.capacity, np.ndarray) or isinstance(testobj.capacity, list):
            capacity = int(testobj.capacity[sim.t] * sim.pars['pop_size'])
        else: 
            capacity = int(testobj.capacity * sim.pars['pop_size'])
        return basic_allocation(capacity, seekers_dict)

    elif mode == 'cross': 
        if not(testobj.system == 'RAT_surv'): 
            raise RuntimeWarning('Cross-sectional allocation mode selected, but testing system is not RAT active surveillance')
        return cross_allocation(testobj, sim)

    elif mode == 'long':
        if not(testobj.system == 'RAT_surv'): 
            raise RuntimeWarning('Longitudinal allocation mode selected, but testing system is not RAT active surveillance')
        return long_allocation(testobj, sim)
        
    else: 
        raise RuntimeError("Invalid allocation policy supplied. Options are: 'proportional', 'basic', 'cross', 'long'")


def basic_allocation(test_capacity, seekers_dict): 
    '''
    Each day, get people who are seeking a PCR test. Randomly determine who gets a test. There is no temporal dimension, people do not stay in any sort of queue if 
    they do not get a test today.  

    Args: 
        test_capacity    (int)       : Number of tests available to allocate
        seekers_dict     (dict)      : Dictionary of people from each criterion who is seeking a test

    Returns: 
        testers          (np.array)  : Array of tester UIDs, where max(length) = test capacity
    '''
    # Unpack test seekers
    test_seekers = set()  # Need uniqueness since some test seekers may satisfy multiple criteria, and we do not want to count them twice
    for value in seekers_dict.values():
        test_seekers.update(value)

    # Determine who actually gets a test
    if len(test_seekers) <= test_capacity:  # All test-seekers can get tested (very unlikely)
        testers = np.fromiter(test_seekers, dtype=int)
    else: 
        testers = np.random.choice(list(test_seekers), test_capacity, replace=False)  # Insufficient capacity, randomly choose a lucky few test-seekers

    return testers

#TODO: Should we guarantee tests for individuals who screen? 
def proportional_allocation(test_capacity, seekers_dict, seekprobs, testobj_name, day): 
    '''
    Compute the power set of criteria. For each power set element (a set of criteria), get the number of people who precisely met these criteria 
    at the eligibility assessment, and divide by the total number of people in the population. This proportion (at most) of tests will be allocated to
    those belonging to the power set element. Note that the power set contains the empty set. This corresponds to people who did not satisfy any criteria, 
    which is untrue of any test seeker. 

    Recall if S = len(X) then power_set(X) has length 2^S

    Args: 
        sim    (simulation) : Simulation object
        test_capacity   (int) : Number of tests to distribute on this day 
        seekers_dict    (dict) : Keys are criteria, values are test-seekers who meet the criteria
        criteria        (list) : 
    '''
    valid_criteria = set([criterion for criterion in seekers_dict.keys() if not(seekprobs[criterion] is None)])  # Do not account for cont_conf and sx sets, since they have subsets that are already accounted for
    
    criteria_powset = compute_criteria_powset(valid_criteria)

    # Get test seekers and the powerset elements they correspond to
    seeker_crit_subset = compute_members(valid_criteria, criteria_powset, seekers_dict)

    # Compute proportions and allocate
    test_proportions = compute_proportions(seeker_crit_subset)

    testers = set()
    for sname in test_proportions.keys(): 
        subset_seekers = np.array(list(seeker_crit_subset[sname]))  # Test seekers who meet the criteria described by subset
        subset_n_tests = int(test_proportions[sname] * test_capacity)

        if subset_n_tests < len(subset_seekers):
            subset_testers = set(np.random.choice(subset_seekers, subset_n_tests, replace=False).tolist())
        else: 
            subset_testers = subset_seekers

        testers.update(subset_testers)
    # Redistribute extras by randomly sampling subsets according to proportions
    extra_tests = test_capacity - len(testers)
    while extra_tests > 0: 

        total_not_tested = sum([len(seeker_crit_subset[sname].difference(testers)) for sname in test_proportions.keys()])
        if total_not_tested == 0: 
            # print(f'{extra_tests} unused {testobj_name} tests on day {day} because not enough people sought a test.')
            break

        probs = [test_proportions[sname] for sname in test_proportions.keys()]
        extra_sname = np.random.choice(list(test_proportions.keys()), 1, p=probs)[0]
        not_tested = seeker_crit_subset[extra_sname].difference(testers)

        if len(not_tested) > 0: 
            testers.add(np.random.choice(list(not_tested)))
            extra_tests -= 1

    return np.fromiter(testers, dtype=int)


def cross_allocation(testobj, sim): 
    '''
    At some frequency, randomly choose a subset of people from the targeted populations.

    '''
    cur = sim.t

    if cur == 0: 
        setattr(testobj, 'all_testers', None)
        setattr(testobj, 'test_days', None)

    if cur % testobj.frequency == 0:  # E.g. if f=7, then we resample after 7 days have passed (0, 1, 2, 3, 4, 5, 6)
        testobj.all_testers = sample_targets(testobj, sim)  # Resample according to defined frequency
        testobj.test_days = calculate_test_days(cur, testobj.frequency, testobj.all_testers)
    
    return testobj.all_testers['all'][testobj.test_days['all'] == cur]  #TODO: Figure out how to test specific/combinations of population layers. Could return a dictionary of those who are testing today. 


def long_allocation(testobj, sim): 
    '''
    On the first day, randomly choose a subset of people from the targeted population(s). At some frequency, re-test each individual. 
    Example: 
        f = 7 days. On day 0, individuals A, B, C are tested. On the 7th day after their initial test (day 7), they get tested again. 
    '''
    cur = sim.t

    if cur == 0: 
        setattr(testobj, 'all_testers', None)
        setattr(testobj, 'test_days', None)

        testobj.all_testers = sample_targets(testobj, sim)  # These individuals will be followed throughout the outbreak for testing
        testobj.test_days = calculate_test_days(cur, testobj.frequency, testobj.all_testers)  # Initialize the first day of testing for each selected individual

    return testobj.all_testers['all'][(cur - testobj.test_days['all']) % testobj.frequency == 0]


def compute_members(valid_criteria, criteria_powset, criteria_dict): 
    '''
    Given a powerset of criteria, figure out which people precisely meet each subset of criteria (no more or less). 

    Args: 
        valid_criteria      (set)  : Set of valid criteria that was used to compute the power set
        criteria_powset     (dict) : Keys are subset names (e.g. 'subset1'), values are the actual criteria subsets
        criteria_dict       (dict) : Keys are criteria, values are sets of people meeting those criteria

    Returns: 
        meets_crit_subsets  (dict) : Keys are subset names, values are people who precisely meet the criteria described in the subset
    '''
    # Get the number of people who precisely meet the criteria specified by each power set element.
    meets_crit_subset = dict()
    for subset_key in criteria_powset.keys(): 
        crit_subset_copy = copy.deepcopy(criteria_powset[subset_key])  # Criteria we want seekers to meet. A copied object avoids bugs with referencing when we do pop(). 
        complement = valid_criteria.difference(crit_subset_copy)  # Criteria we don't want seekers to meet. crit_subset_copy should be a subset of valid_criteria.

        members = criteria_dict[crit_subset_copy.pop()]  # Choose any criterion to start with, pop() can only happen after line above 

        # Get people who meet at least the prescribed subset criteria 
        if len(crit_subset_copy) > 0:  # May be zero if it only contained one element initially, in that case then we already have all members 
            for criterion in crit_subset_copy: 
                members = members.intersection(criteria_dict[criterion])  # Remove anyone in 

        # Remove anyone who meets other criteria
        if len(complement) > 0:  # May be zero when the crit_subset contains all criteria, in that case then there is nothing to remove
            for criterion in complement: 
                members = members - members.intersection(criteria_dict[criterion])

        meets_crit_subset[subset_key] = members  # Meets all and only the subset criteria

    return meets_crit_subset

    
def compute_proportions(crit_subset_dict): 
    '''
    Determine which proportion of tests should be delegated to each subset element of the criteria powerset

    Args: 
        eligible_crit_subset  (dict) : Keys are powerset subset names, values are corresponding test-eligible people who meet the criteria subset
    '''
    proportions = dict()
    n_eligible = 0
    for subset_ppl in crit_subset_dict.values(): 
        n_eligible += len(subset_ppl)

    for sname in crit_subset_dict.keys(): 
        if n_eligible != 0: 
            proportions[sname] = len(crit_subset_dict[sname]) / n_eligible
        else: 
            proportions[sname] = 0  # If no one is eligible, then no one should be given a test

    return proportions
    

def compute_criteria_powset(valid_criteria): 
    ''' 
    Given a set of criteria, compute and return the powerset.
    
    Args: 
        valid_criteria  (set)   : Set of valid criteria to compute the power set from. 

    Returns: 
        criteria_powset (dict)  : Keys identify subsets, values are all subsets of valid_criteria
    '''
    criteria_powset = dict()
    for i in range(1, len(valid_criteria) + 1):  # We skip the empty set, since this corresponds to someone satisfying no criteria
        subsets_size_i = list(combinations(valid_criteria, i))
        for j in range(len(subsets_size_i)): 
            criteria_powset['subset' + str(i) + str(j)] = set(subsets_size_i[j])

    return criteria_powset


def get_targets(testobj, sim): 
    '''
    Get target populations. Eventually this may expand into its own class (akin to TestEligible) if we want more detailed targetting (e.g. classroom level)

    NOTE: For now, testing is only designed around targetting one population layer. Not multiple, even though the infrastructure is partly there. 
    '''
    target_layers = testobj.population 

    target_uid = dict()
    for l in target_layers: 
        if l == "all":
            target_uid[l] = sim.people.uid
        elif l == "household": 
            target_uid[l] = sim.people.contacts['h'].members
        elif l == "school": 
            target_uid[l] = sim.people.contacts['s'].members
        elif l == "workplace": 
            target_uid[l] = sim.people.contacts['w'].members
        else: 
            raise RuntimeWarning("Invalid population layer supplied. Options are: 'all', 'household', 'school', 'workplace'")

    return target_uid


def sample_targets(testobj, sim):
    '''
    Randomly select some proportion of each targeted population layer. 

    NOTE: To allow multiple selected population layers, we need to define what happens when a person is selected in multiple layers
    '''
    target_uid = get_targets(testobj, sim)
    p = testobj.participation

    sampled_uid = dict()
    if not(testobj.age_report is None):  # Sample according to custom defined distributions NOTE: This feature is currently only developed for the case that we are ONLY targetting 'all'
        age_group_members = get_age_group_members(testobj, sim, target_uid)
        age_distribution = testobj.all_age_reports[testobj.age_report]['distribution']

        # if sim.t == 0: 
        #     setattr(testobj, 'age_distribution_history', dict())
        #     for interval in age_group_members.keys(): 
        #         testobj.age_distribution_history[interval] = []

        for l in target_uid.keys():  # 'all' is the only option
            sampled_uid[l] = set()
            for age_group, interval in zip(age_distribution, age_group_members.keys()):
                n_tests = int(age_group[2] * p[l] * len(target_uid[l]))  # Age prevalence * participation * population size
                sampled_uid[l].update(set(np.random.choice(age_group_members[interval], n_tests, replace=False).tolist()))
                # testobj.age_distribution_history[interval].append(n_tests)
            
            sampled_uid[l] = np.array(list(sampled_uid[l]))

    else:  # Pure random sampling
        for l in target_uid.keys(): 
            sampled_uid[l] = np.random.choice(target_uid[l], int(len(target_uid[l]) * p[l]), replace=False)

    return sampled_uid


def calculate_test_days(t, f, sampled_uid):
    '''
    Args: 
        t             (int)      : Current time in the simulation
        f             (int)      : Frequency of testing
        sampled_uid   (dict)     : Dictionary containing arrays of sampled individuals in different population layers
        mode          (str)      : Whether to use cross sectional or longitudinal strategy
    '''
    test_days = dict()

    for l in sampled_uid.keys(): 
        test_days[l] = t + np.random.choice(f, len(sampled_uid[l]))  # Maximum choice is f - 1
    return test_days


def get_age_group_members(testobj, sim, target_uid):
    '''
    Args: 
        testobj         (surveillance) : RAT surveillance object 
        target_uid    (dict)         : Dictionary of targetted individuals for different population layers
    '''
    age_distribution = testobj.all_age_reports[testobj.age_report]['distribution']

    # Find everyone who belongs to each age group
    age_group_members = dict()
    for age_group in age_distribution:
        if age_group[1] != -1:  # -1 symbolizes "greater than"
            interval = str(int(age_group[0])) + '-' + str(int(age_group[1]))
            age_group_members[interval] = np.logical_and(sim.people.age >= age_group[0], sim.people.age <= age_group[1]).nonzero()[0]
        else: 
            interval = str(int(age_group[0])) + '+'
            age_group_members[interval] = (sim.people.age >= age_group[0]).nonzero()[0]

    # NOTE: This feature is currently not being used. However, the architecture may be useful in the future. 
    # Find everyone who ALSO is a target - this is done in case the target layers do not contain everyone in the population
    if not('all' in target_uid.keys()):
        combined_targets = reduce(np.union1d, (targets for targets in target_uid.values()))  # Compute union of all layers in target_uid
        for interval in age_group_members.keys(): 
            age_group_members[interval] = np.intersect1d(age_group_members[interval], combined_targets)
    
    return age_group_members
