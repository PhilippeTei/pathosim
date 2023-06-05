import numpy as np
from .background_ILI import * 

###############################
##### Diagnostic testing ######
###############################

def basic_PCR_test(surv, sim, uids, se, sp, analysis=True): 
    '''
    Administer PCR tests to a list of users. Assume people have COVID-19 only if they are exposed (i.e. any time between infection to recovery). 
    Conduct analysis if desired. Record whether they were symptomatic at the time of testing.

    Args: 
        sim                (sim): Simulation object
        uids          (np.array): Agent indices to test
        se               (float): Test sensitivity 
        sp               (float): Test specificity

    Returns: 
        test_results  (np.array): Boolean array of test results
    '''
    if sim.t == 0: 
        setattr(surv, 'neg_was_sx', np.full(sim.pars['pop_size'], 0))  # 0: Never tested before, 1: Tested and was sx, 2: Tested and was asx
        setattr(surv, 'pos_was_asx', np.full(sim.pars['pop_size'], 0))

    test_results = np.full(len(uids), False)
    test_results[sim.people.exposed[uids]] = np.random.uniform(0, 1, sum(sim.people.exposed[uids])) < se
    test_results[~sim.people.exposed[uids]] = np.random.uniform(0, 1, sum(~sim.people.exposed[uids])) > sp

    # for i in range(len(uids)): 
    #     if sim.people.exposed[uids[i]]: 
    #         test_results[i] = np.random.uniform(0, 1) < se
    #     else: 
    #         test_results[i] = np.random.uniform(0, 1) > sp

    # Identify who is sx/asx at the time of testing and record
    neg_asx = uids[np.logical_and(~test_results, np.logical_and(~sim.people.symptomatic[uids], ~sim.people.symptomatic_ILI[uids]))]
    neg_sx = uids[np.logical_and(~test_results, np.logical_or(sim.people.symptomatic[uids], sim.people.symptomatic_ILI[uids]))]
    pos_sx = uids[np.logical_and(test_results, np.logical_or(sim.people.symptomatic[uids], sim.people.symptomatic_ILI[uids]))]
    pos_asx = uids[np.logical_and(test_results, np.logical_and(~sim.people.symptomatic[uids], ~sim.people.symptomatic_ILI[uids]))]

    surv.neg_was_sx[neg_asx] = 2
    surv.neg_was_sx[neg_sx] = 1
    surv.pos_was_asx[pos_sx] = 1
    surv.pos_was_asx[pos_asx] = 2


    if analysis: 
        tests_consumed(surv, sim, len(uids))

    return test_results



def basic_RAT_test(surv, sim, uids, se, sp, analysis=True): 
    '''
    Administer RAT tests to a list of users. Assume people have COVID-19 only if they are exposed (i.e. any time between infection to recovery). 
    Conduct analysis if desired. Record whether they were symptomatic or not at the time of testing.

    Args: 
        sim                (sim): Simulation object
        uids          (np.array): Agent indices to test
        se               (float): Test sensitivity 
        sp               (float): Test specificity

    Returns: 
        test_results  (np.array): Boolean array of test results
    '''
    if sim.t == 0: 
        setattr(surv, 'neg_was_sx', np.full(sim.pars['pop_size'], 0))  # 0: Never tested before, 1: Tested and was sx, 2: Tested and was asx
        setattr(surv, 'pos_was_asx', np.full(sim.pars['pop_size'], 0))

    test_results = np.full(len(uids), False)
    test_results[sim.people.exposed[uids]] = np.random.uniform(0, 1, sum(sim.people.exposed[uids])) < se
    test_results[~sim.people.exposed[uids]] = np.random.uniform(0, 1, sum(~sim.people.exposed[uids])) > sp

    # for i in range(len(uids)): 
    #     if sim.people.exposed[uids[i]]: 
    #         test_results[i] = np.random.uniform(0, 1) < se
    #     else: 
    #         test_results[i] = np.random.uniform(0, 1) > sp
    
    # Identify who is sx/asx at the time of testing and record
    neg_asx = uids[np.logical_and(~test_results, np.logical_and(~sim.people.symptomatic[uids], ~sim.people.symptomatic_ILI[uids]))]
    neg_sx = uids[np.logical_and(~test_results, np.logical_or(sim.people.symptomatic[uids], sim.people.symptomatic_ILI[uids]))]
    pos_sx = uids[np.logical_and(test_results, np.logical_or(sim.people.symptomatic[uids], sim.people.symptomatic_ILI[uids]))]
    pos_asx = uids[np.logical_and(test_results, np.logical_and(~sim.people.symptomatic[uids], ~sim.people.symptomatic_ILI[uids]))]

    surv.neg_was_sx[neg_asx] = 2
    surv.neg_was_sx[neg_sx] = 1
    surv.pos_was_asx[pos_sx] = 1
    surv.pos_was_asx[pos_asx] = 2

    # Test cases
    # print(len(neg_sx), len(neg_asx), sum(~test_results))
    # print(len(neg_sx) + len(neg_asx) == sum(test_results))
    # print(len(cur_sx) + len(cur_asx) == len(uids))  # True every time
    # print(sum(sim.people.exposed[uids]) + sum(~sim.people.exposed[uids]) == len(uids))


    if analysis: 
        tests_consumed(surv, sim, len(uids))

    return test_results



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
###### Testing analysis #######
###############################


def tests_consumed(surv, sim, tests): 
    if sim.t == 0: 
        setattr(surv, 'tests_consumed', [])

    surv.tests_consumed.append(tests)