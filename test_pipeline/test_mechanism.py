import numpy as np 

def LOD_test(testobj, sim, uids): 
    '''
    Administer tests to a list of users. Assume people have COVID-19 only if they are exposed (i.e. any time between infection to recovery). 
    Conduct analysis if desired. Record whether they were symptomatic at the time of testing.

    Args: 
        sim                (sim): Simulation object
        uids          (np.array): Agent indices to test
        se               (float): Test sensitivity 
        sp               (float): Test specificity

    Returns: 
        test_results  (np.array): Boolean array of test results
    '''
    LOD = testobj.LOD
    test_results = sim.people.viral_load[uids] > LOD
    # Identify who is sx/asx at the time of testing and record
    neg_asx = uids[np.logical_and(~test_results, np.logical_and(~sim.people.symptomatic[uids], ~sim.people.symptomatic_ILI[uids]))]
    neg_sx = uids[np.logical_and(~test_results, np.logical_or(sim.people.symptomatic[uids], sim.people.symptomatic_ILI[uids]))]
    pos_sx = uids[np.logical_and(test_results, np.logical_or(sim.people.symptomatic[uids], sim.people.symptomatic_ILI[uids]))]
    pos_asx = uids[np.logical_and(test_results, np.logical_and(~sim.people.symptomatic[uids], ~sim.people.symptomatic_ILI[uids]))]

    testobj.neg_was_sx[neg_asx] = 2
    testobj.neg_was_sx[neg_sx] = 1
    testobj.pos_was_asx[pos_sx] = 1
    testobj.pos_was_asx[pos_asx] = 2

    return test_results

def basic_test(testobj, sim, uids): 
    '''
    Administer tests to a list of users. Assume people have COVID-19 only if they are exposed (i.e. any time between infection to recovery). 
    Conduct analysis if desired. Record whether they were symptomatic at the time of testing.

    Args: 
        sim                (sim): Simulation object
        uids          (np.array): Agent indices to test
        se               (float): Test sensitivity 
        sp               (float): Test specificity

    Returns: 
        test_results  (np.array): Boolean array of test results
    '''
    se = testobj.sensitivity
    sp = testobj.specificity

    test_results = np.full(len(uids), False)
    test_results[sim.people.exposed[uids]] = np.random.uniform(0, 1, sum(sim.people.exposed[uids])) < se
    test_results[~sim.people.exposed[uids]] = np.random.uniform(0, 1, sum(~sim.people.exposed[uids])) > sp

    # Identify who is sx/asx at the time of testing and record
    neg_asx = uids[np.logical_and(~test_results, np.logical_and(~sim.people.symptomatic[uids], ~sim.people.symptomatic_ILI[uids]))]
    neg_sx = uids[np.logical_and(~test_results, np.logical_or(sim.people.symptomatic[uids], sim.people.symptomatic_ILI[uids]))]
    pos_sx = uids[np.logical_and(test_results, np.logical_or(sim.people.symptomatic[uids], sim.people.symptomatic_ILI[uids]))]
    pos_asx = uids[np.logical_and(test_results, np.logical_and(~sim.people.symptomatic[uids], ~sim.people.symptomatic_ILI[uids]))]

    testobj.neg_was_sx[neg_asx] = 2
    testobj.neg_was_sx[neg_sx] = 1
    testobj.pos_was_asx[pos_sx] = 1
    testobj.pos_was_asx[pos_asx] = 2

    return test_results

# Not being used
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
