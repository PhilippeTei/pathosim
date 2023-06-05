# Functions for assessing test eligibility
import numpy as np
import numba as nb
import copy
import json
from infection import utils as cvu


def check_criteria(base_testobj, testobj, sim, **crit_kwargs): 
        '''
        Given the criterion names, call respective criterion functions to identify who meets each criterion. Apply restrictions at the end.
        '''
        all_eligible = dict()  # Will be a dictionary of sets
        RAT = testobj.RAT
        criteria = base_testobj.criteria

        for name in criteria: 
            # Diagnostic
            if name == 'cont_conf': 
                all_eligible['cont_conf'] = contact_conf(sim)
            if name == 'cont_conf_quar':
                all_eligible['cont_conf_quar'] = contact_conf_quar(sim)
            if name == 'cont_vuln': 
                all_eligible['cont_vuln'] = contact_vuln(sim)

            if name == 'cont_sx': 
                all_eligible['cont_sx'] = contact_sx(sim, mode='all')
            if name == 'cont_cs_sx': 
                all_eligible['cont_cs_sx'] = contact_sx(sim, mode='specific')
            if name == 'cont_ncs_sx':
                all_eligible['cont_ncs_sx'] = contact_sx(sim, mode='nonspecific')
            
            if name == 'sx': 
                all_eligible['sx'] = sx(sim, mode='all')
            if name == 'cs_sx': 
                all_eligible['cs_sx'] = sx(sim, mode='specific')
            if name == 'ncs_sx': 
                all_eligible['ncs_sx'] = sx(sim, mode='nonspecific')
            if name == 'sx_quar':
                all_eligible['sx_quar'] = sx_quar(sim)
            if name == 'sev_crit':
                all_eligible['sev_crit'] = sev_crit(sim)

            if name == 'other':
                if 'other' in crit_kwargs.keys(): 
                    all_eligible['other'] = other(sim, **crit_kwargs['other'])
                else:
                    all_eligible['other'] = other(sim)

            # NOTE: Assumed that the first three if statements in this block will only be triggered by PCR. Fourth statement will only be triggered by RAT. 
            # May adjust this assumption later. For example, confirmatory testing given a positive RAT is useful during low prevalence. 
            if name == 'pos_RAT': 
                all_eligible['pos_RAT'] = RAT_testers(sim, RAT, mode='pos_RAT')
            if name == 'pos_RAT_asx': 
                all_eligible['pos_RAT_asx'] = RAT_testers(sim, RAT, mode='pos_RAT_asx')
            if name == 'neg_RAT_sx': 
                all_eligible['neg_RAT_sx'] = RAT_testers(sim, RAT, mode='neg_RAT_sx')
            if name == 'neg_RAT': 
                all_eligible['neg_RAT'] = RAT_testers(sim, testobj, mode='neg_RAT')

            # Screening
            if name == 'work': 
                if 'work' in crit_kwargs.keys():
                    all_eligible['work'] = work(sim, base_testobj, **crit_kwargs['work'])
                else:
                    all_eligible['work'] = work(sim, base_testobj)

            # Smartwatch
            if name == 'swa_RAT':
                all_eligible['swa_RAT'] = sw_alerted_RAT(sim)
            if name == 'swa_PCR':
                all_eligible['swa_PCR'] = sw_alerted_PCR(sim)

        # Apply retest restrictions to those seeking diagnostic testing and finally return
        if base_testobj.system == 'Diagnostic': 
            return apply_restrictions(testobj, all_eligible) 
        else: 
            return all_eligible


def contact_conf(sim):
    '''
    Select people who have contacted confirmed case today
    '''
    contacted_inds = cvu.true(sim.people.date_known_contact == sim.t)
    contacts = set(contacted_inds.tolist())
    return contacts

def contact_conf_quar(sim, d=7):
    '''
    Select people who have contacted confirmed case in last 7 days and they are not quarantined)
    '''
    prev_contacted = np.logical_and(sim.t - sim.people.date_known_contact <= d, sim.t - sim.people.date_known_contact >= 0)  # logical_and not strictly necessary
    not_quarantined = sim.people.quarantined == False
    return set(cvu.true(np.logical_and(prev_contacted, not_quarantined)).tolist())

def contact_conf_old(sim, testobj, d=3): 
        '''
        Old version of contact_conf.
        Get people who were in contact with a confirmed case in the past d days. Select a subset who seeks a test. 
            * Note that this criterion inherently requires people to have tested positive previously. If this is the only eligibility criterion 
              used to determine test seekers and testers, then no one will end up testing. 

        Args: 
            sim  (simulation)  : Simulation object 
            testobj            : Test regime object
            d    (int)         : Used to get people who were confirmed positive in past d days

        Returns: 
            seekers (set)       : People who satisfy the criterion and will seek a test
        '''
        cur = sim.t

        # Get people who have tested positive previously (and are not dead) with the current Test regime system
        pos_inds = cvu.true(np.logical_and(testobj.date_positive >= 0, ~sim.people.dead))
        conf_uid = pos_inds[np.logical_and(cur - testobj.date_positive[pos_inds] <= d, cur - testobj.date_positive[pos_inds] >= 0)]

        # Get people who have tested positive previously (and are not dead) with the complement surveillance system
        if testobj.hybrid and (testobj.system == 'RAT_disc' or testobj.system == 'RAT_surv'): 
            PCR_pos_inds = cvu.true(np.logical_and(testobj.PCR.date_positive >= 0, ~sim.people.dead))
            PCR_conf_uid = PCR_pos_inds[np.logical_and(cur - testobj.PCR.date_positive[PCR_pos_inds] <= d, cur - testobj.PCR.date_positive[PCR_pos_inds] >= 0)]
            conf_uid = np.union1d(conf_uid, PCR_conf_uid)
        elif testobj.hybrid and (testobj.system == 'PCR_disc' or testobj.system == 'PCR_surv'):
            RAT_pos_inds = cvu.true(np.logical_and(testobj.RAT.date_positive >= 0, ~sim.people.dead))
            RAT_conf_uid = RAT_pos_inds[np.logical_and(cur - testobj.RAT.date_positive[RAT_pos_inds] <= d, cur - testobj.RAT.date_positive[RAT_pos_inds] >= 0)]
            conf_uid = np.union1d(conf_uid, RAT_conf_uid)
        elif not(testobj.hybrid) and not(testobj.RAT is None):
            raise RuntimeError('RAT complement was provided, but PCR testobj is not hybrid.')
        elif not(testobj.hybrid) and not(testobj.PCR is None):
            raise RuntimeError('PCR complement was provided, but RAT testobj is not hybrid.')

        # For each layer, get contacts of these people
        contacts = set()  # One person may have contacted several confirmed positives, which may increase probability of test-seeking. We are currently not accounting for this.
        for layer in sim.people.contacts.keys(): 
            contacts.update(sim.people.contacts[layer].find_contacts(conf_uid, as_array=False))  # Adds empty list if agent i does not exist in the layer, or if layer is empty. Sped up by Numba.

        return contacts

def contact_vuln(sim):  
    '''
    Get people who are in contact with vulnerable people in the household. Vulnerability defined in terms of odds ratio of susceptibility. 
    '''
    # Get people who are vulnerable (and not dead) 
    vuln_uid = cvu.true(np.logical_and(sim.people.rel_sus > 1.0, ~sim.people.dead))  # See parameters.py for assignment of susceptibility odds ratios by age bracket

    # Get their household contacts. It may be more realistic to only get some subset of the household contacts.
    return sim.people.contacts['h'].find_contacts(vuln_uid, as_array=False)  # Sped up by Numba

def contact_sx(sim, d=3, mode='all'): 
    '''
    Get people who were in contact with a symptomatic person, in the last d days. 
    These people have contacted someone with either COVID specific or nonspecific symptoms.
    This function only works because contact networks are contacted every single day. 
    
    Args: 
        sim            (sim) : Simulation object 
        d              (int) : Contacts of symptomatic people, in the last d days, satisfy this criterion
        mode           (str) : Choose return option

    Returns: 
        contacts       (set) : Set of contacts. 
                                    If mode == 'all', returns set of people who have contacted at least one symptomatic person. 
                                    If mode == 'specific', returns set of people who have contacted at least one person with specific COVID symptoms.
                                    If mode == 'nonspecific', returns set of people who have contacted at least one person with nonspecific COVID symptoms. 
    '''
    return_options = ['all','specific', 'nonspecific']

    if not(mode in return_options): 
        raise RuntimeError('Return choice is invalid, the options are: "all", "specific", "nonspecific"')

    all_sx = cvu.true(np.logical_or(sim.people.symptomatic, sim.people.symptomatic_ILI))  # Anyone who is currently symptomatic with COVID or ILI

    if mode == 'all': 
        all_contacts = set() 
        for layer in sim.people.contacts.keys(): 
            all_contacts.update(sim.people.contacts[layer].find_contacts(all_sx, as_array=False))  # as_array=True will sort indices, which is unnecessary
        return all_contacts

    elif mode == 'specific':
        return stratify_contacts_sx(sim, all_sx, 's')

    elif mode == 'nonspecific': 
        return stratify_contacts_sx(sim, all_sx, 'ns')

def sx(sim, mode='all'): 
    '''
    Get people who are currently displaying symptoms. These people either have COVID specific or nonspecific symptoms.

    Args: 
        sim          (sim) : Simulation object 
        mode         (str) : Choose return option 

    Returns: 
        symptomatic  (set) : Set of people who are symptomatic with COVID or ILI. 
                                If mode == 'all', returns set of people who have any symptoms
                                If mode == 'specific', returns set of people who have COVID specific symptoms
                                If mode == 'nonspecific', returns set of people who have non COVID specific symptoms
    '''
    return_options = ['all', 'specific', 'nonspecific']

    if not(mode in return_options): 
        raise RuntimeError('Return choice is invalid, the options are: "all", "specific, "nonspecific"')

    all_sx = cvu.true(np.logical_or(sim.people.symptomatic, sim.people.symptomatic_ILI))

    if mode == 'all': 
        return set(all_sx.tolist())

    elif mode == 'specific': 
        return set(stratify_sx(sim.people.symptom_types, all_sx, 's').tolist())

    elif mode == 'nonspecific': 
        return set(stratify_sx(sim.people.symptom_types, all_sx, 'ns').tolist())

def sx_quar(sim, d=6):
    '''
    Get people who are currently quarantined and the date of symptom onset was more than d days ago.
    Asymptomatic individuals will never invoke this criterion (date_symptomatic and date_symptomatic_ILI are both np.nan).
    '''
    sx_covid = np.logical_and(sim.t - sim.people.date_symptomatic >= d, sim.people.symptomatic)
    sx_ILI = np.logical_and(sim.t - sim.people.date_symptomatic_ILI >= d, sim.people.symptomatic_ILI)
    sx = np.logical_or(sx_covid, sx_ILI)
    sx_quar = cvu.true(np.logical_and(sim.people.quarantined == True, sx))
    return set(sx_quar.tolist())

def sev_crit(sim): 
    '''
    Get people who are currently severe or critical
    '''
    sev_crit = np.logical_or(sim.people.severe, sim.people.critical)
    return set(sev_crit.tolist())

def RAT_testers(sim, RAT, mode='pos_RAT', d=3): 
    '''
    Interface between passive PCR and active RAT Test regime in hybrid model. Get people who either: 
        1) Received positive RAT in past d days 
        2) Received positive RAT in past d days and is asx OR received negative RAT in past d days and is sx
    '''
    if RAT is None: 
        raise RuntimeError('RAT_testers called but no RAT system was provided.')
    cur = sim.t
    date_positive = RAT.date_positive
    date_negative = RAT.date_negative

    # Get people who had a positive/negative RAT previously (and are not dead) 
    pos_inds = cvu.true(np.logical_and(date_positive >= 0, ~sim.people.dead))
    neg_inds = cvu.true(np.logical_and(date_negative >= 0, ~sim.people.dead))

    # Get people who were confirmed positive/negative by a RAT test in the last d days
    pos_d_inds = pos_inds[np.logical_and(cur - date_positive[pos_inds] <= d, cur - date_positive[pos_inds] >= 0)]  # 2nd part of "and" statement is needed since reporting date may be in the future
    neg_d_inds = neg_inds[np.logical_and(cur - date_negative[neg_inds] <= d, cur - date_negative[neg_inds] >= 0)]

    # print(len(neg_d_inds) == len(neg_d_inds[RAT.neg_was_sx[neg_d_inds] == 1]) + len(neg_d_inds[RAT.neg_was_sx[neg_d_inds] == 2]))  # Sanity check

    if mode == 'pos_RAT': 
        # Return set of individuals who got a positive RAT test in past d days
        return set(pos_d_inds.tolist())

    elif mode == 'pos_RAT_asx':  
        # Return set of individuals who got a positive RAT test in past d days and was asymptomatic at time of testing
        return set(pos_d_inds[RAT.pos_was_asx[pos_d_inds] == 2].tolist())  # to
        
    elif mode == 'neg_RAT_sx': 
        # Return set of individuals who got a negative RAT test in past d days and was symptomatic at time of testing
        return set(neg_d_inds[RAT.neg_was_sx[neg_d_inds] == 1].tolist())

    elif mode == 'neg_RAT': 
        return set(neg_d_inds.tolist())

def other(sim, p=0.05):
    return set(np.random.choice(sim.pars['pop_size'], int(p * sim.pars['pop_size']), replace=False).tolist())

def work(sim, base_testobj, p=0.10):
    '''
    Get people who have to do mandatory workplace testing. Randomly select some number of workplaces. Every week, these people will seek a test. 
    To consider: Do we want them to have priority access to testing resources? 

    Args: 
        p    (float) : Proportion of workplaces to undergo mandatory testing
    '''
    if sim.t == 0:
        # Identify which workplaces are chosen to do mandatory testing, and save the workers
        test_wpid = np.random.choice(len(base_testobj.workplaces), int(base_testobj.n_workplaces * p), replace=False)

        worker_uids = np.array([], dtype=int)
        for wpid in test_wpid:
            worker_uids = np.union1d(worker_uids, base_testobj.workplaces[wpid]['member_uids']).astype(int)
        
        worker_test_dates = np.random.choice(7, len(worker_uids))

        setattr(base_testobj, 'wuid', worker_uids)  # Set of all workers
        setattr(base_testobj, 'wuid_td', worker_test_dates)  # Days of the week they are assigned to test

    # Remove anyone who died
    base_testobj.wuid_td = base_testobj.wuid_td[~sim.people.dead[base_testobj.wuid]]  # Order matters
    base_testobj.wuid = base_testobj.wuid[~sim.people.dead[base_testobj.wuid]]

    # Same workplaces, and thus workers, test each time
    return set(base_testobj.wuid[base_testobj.wuid_td == sim.t % 7].tolist())

def sw_alerted_RAT(sim):
    '''
    Get people who received a smartwatch alert. Only call this for the RAT test object.
    The functions are similar, but they are managed differently so that we can specify differing test seeking probabilities more explicitly.
    '''
    return set(cvu.true(sim.people.alerted).tolist())

def sw_alerted_PCR(sim):
    '''
    Get people who received a smartwatch alert. Only call this for the PCR test object. 
    '''
    return set(cvu.true(sim.people.alerted).tolist())


#### Helper functions ####

# This implementation is not flexible to different symptom configurations 
@nb.njit((nb.boolean[:], nb.boolean[:], nb.boolean[:], nb.int64[:]))
def get_cs_sx(s0, s1, s2, sx_uid):
    '''
    Get symptomatic individuals who have at least one of the three COVID-specific symptoms.
    Associativity of 'or': s0 V (s1 V s2) is equivalent to (s0 V s1) V s2.
    '''
    return sx_uid[np.logical_or(s0, np.logical_or(s1, s2))]


@nb.njit((nb.boolean[:], nb.boolean[:], nb.boolean[:], nb.int64[:]))
def get_ncs_sx(s0, s1, s2, sx_uid): 
    return sx_uid[~np.logical_or(s0, np.logical_or(s1, s2))]

def stratify_sx(symptom_types, sx_uid, mode): 
    '''
    Helper function. Stratify those who have symptoms by whether their symptoms are specific or nonspecific. 

    Args: 
        symptom_types      (np.array)  : 2D array where rows are symptoms and columns are uids. Elements are bools (whether individual has a symptom)
        sx_uid             (np.array)  : UID of symptomatic people
        mode               (str)       : Whether to get specific or nonspecific
    '''
    # Define specific and nonspecific symptoms. Specific symptoms have higher associated probability of test seeking. 

    if mode == 's': 
        return get_cs_sx(symptom_types[0, sx_uid], symptom_types[1, sx_uid], symptom_types[2, sx_uid], sx_uid)
        # cs_uid = sx_uid[np.logical_and(symptom_types[0, sx_uid], np.logical_or(symptom_types[1, sx_uid], symptom_types[2, sx_uid]))]  # Fever AND Cough or Sore Throat
        # return set(cs_uid.tolist())

    elif mode == 'ns': 
        return get_ncs_sx(symptom_types[0, sx_uid], symptom_types[1, sx_uid], symptom_types[2, sx_uid], sx_uid)
        # ncs_uid = sx_uid[~np.logical_and(symptom_types[0, sx_uid], np.logical_or(symptom_types[1, sx_uid], symptom_types[2, sx_uid]))]  # not(Fever AND Cough or Sore Throat)
        # return set(ncs_uid.tolist())


def stratify_contacts_sx(sim, sx_uid, mode): 
    '''
    Helper function. Get contacts of people who are symptomatic with either COVID-specific or non COVID-specific symptoms.
    Same function exists in symptoms.py

    Args: 
        sim          (sim) : Simulation object
        sx_uid       (set) : Set of people with symptoms (either due to COVID or ILI)

    Returns: 
        cs_contacts  (set) : Set of people who have contacted someone with COVID-specific symptoms
        ncs_contacts (set) : Set of people who have contacted someone with non-COVID-specific symptoms
    '''
    if mode == 's': 
        cs_uid = stratify_sx(sim.people.symptom_types, np.fromiter(sx_uid, int), 's')
        cs_contacts = set()
        for layer in sim.people.contacts.keys(): 
            cs_contacts.update(sim.people.contacts[layer].find_contacts(cs_uid, as_array=False))
        return cs_contacts

    elif mode == 'ns': 
        ncs_uid = stratify_sx(sim.people.symptom_types, np.fromiter(sx_uid, int), 'ns')
        ncs_contacts = set()
        for layer in sim.people.contacts.keys(): 
            ncs_contacts.update(sim.people.contacts[layer].find_contacts(ncs_uid, as_array=False))
        return ncs_contacts


def apply_restrictions(testobj, eligible_dict): 
    '''
    Remove individuals from the eligibility dictionary according to a set of hard restrictions. 
    Used to restrict retest frequency.
    '''
    restricted = set((testobj.retest_tracker > 0).nonzero()[0].tolist())  # Negative values are fine
   
    # Iterate through each diagnostic criteria and remove anyone who is restricted. Screening individuals are not restricted
    for criterion in testobj.di_criteria:
        eligible_dict[criterion] = eligible_dict[criterion] - eligible_dict[criterion].intersection(restricted)

    return eligible_dict

    

# # TODO: Speed this up with Numba, try to see if an array implementation is helpful
# def stratify_sx(symptom_types, sx_uid, mode): 
#     '''
#     Helper function. Stratify those who have symptoms by whether their symptoms are specific or nonspecific. 

#     Args: 
#         symptom_types      (np.array)  : Array of sets, containing the symptoms each person has
#         sx_uid             (np.array)  : UID of symptomatic people
#         mode               (str)       : Whether to get specific or nonspecific
#     '''
#     # Define specific and nonspecific symptoms. Specific symptoms have higher associated probability of test seeking. 

#     if mode == 's': 
#         cs_uid = set() 
#         for uid in sx_uid: 
#             if ('fever' in symptom_types[uid]) and (('cough' in symptom_types[uid]) or ('sore_throat' in symptom_types[uid])):
#                 cs_uid.add(uid)
#         return cs_uid

#     elif mode == 'ns': 
#         ncs_uid = set() 
#         for uid in sx_uid: 
#             if not(('fever' in symptom_types[uid]) and (('cough' in symptom_types[uid]) or ('sore_throat' in symptom_types[uid]))):
#                 ncs_uid.add(uid)
#         return ncs_uid

# def stratify_contacts_sx(sim, sx_uid, mode): 
#     '''
#     Helper function. Get contacts of people who are symptomatic with either COVID-specific or non COVID-specific symptoms.
#     Same function exists in symptoms.py

#     Args: 
#         sim          (sim) : Simulation object
#         sx_uid       (set) : Set of people with symptoms (either due to COVID or ILI)

#     Returns: 
#         cs_contacts  (set) : Set of people who have contacted someone with COVID-specific symptoms
#         ncs_contacts (set) : Set of people who have contacted someone with non-COVID-specific symptoms
#     '''
#     if mode == 's': 
#         cs_uid = stratify_sx(sim.people.symptom_types, sx_uid, 's')
#         cs_contacts = set()
#         for layer in sim.people.contacts.keys(): 
#             cs_contacts.update(sim.people.contacts[layer].find_contacts(np.fromiter(cs_uid, dtype=int), as_array=False))
#         return cs_contacts

#     elif mode == 'ns': 
#         ncs_uid = stratify_sx(sim.people.symptom_types, sx_uid, 'ns')
#         ncs_contacts = set()
#         for layer in sim.people.contacts.keys(): 
#             ncs_contacts.update(sim.people.contacts[layer].find_contacts(np.fromiter(ncs_uid, dtype=int), as_array=False))
#         return ncs_contacts


