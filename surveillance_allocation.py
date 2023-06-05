import numpy as np
import copy
import json

from itertools import combinations
from functools import reduce
from . import utils as cvu
from . import symptoms as symptoms
from infection.background_ILI import * 

# Abbreviations: 
# sx: symptomatic 
# cs: COVID-specific 
# ncs: non-COVID-specific


###############################
##### Check test criteria #####
###############################
#TODO: We may want to run test seekers through other filters, such as swab delay (delay from symptomaticity to swab)
#TODO: We should move wearable distribution to this file

class TestEligible:
    '''
    Used to check who can get a PCR test in the population. 
        1) Define a set of criteria
        2) Check who meets the criteria -> all_eligible stores a dictionary of people who satisfy each criterion
    
    Args: 
        criteria      (list)    : Defines the criteria to check
        RAT           (RAT) : Optional argument used to interface PCR and RAT testing
    '''

    def __init__(self, surv, sim, criteria=None, RAT=None): 
        # valid_criteria = {
        #     'personal': ['cont_conf', 'cont_sx', 'cont_vuln', 'sx_COVID', 'sx_ILI' ],
        #     'mandatory': ['travel', 'work']}

        if criteria is None: 
            self.criteria = [
                'cont_conf', 
                'cont_vuln',
                'cont_sx', 
                'cont_cs_sx',
                'cont_ncs_sx',
                'sx',
                'cs_sx',
                'ncs_sx',
                'pos_RAT',
                'pos_RAT_asx',
                'neg_RAT_sx',
                'travel',
                'work']
            
        else: 
            self.criteria = criteria

        self.RAT = RAT
        self.all_eligible = self.check_criteria(surv, sim)


    def check_criteria(self, surv, sim): 
        '''
        Given the criterion names, call respective criterion functions to identify who meets each criterion. 

        Args: 
            sim     (simulation)    : Initialized simulation oboject
            surv    (surveillance)  : The PCR surveillance system being used 
        '''
        all_eligible = dict()  # Will be a dictionary of sets

        for name in self.criteria:  # Not sure if there's a cleaner way to do this.
            if name == 'cont_conf': 
                all_eligible['cont_conf'] = self.contact_conf(sim, surv)
            if name == 'cont_vuln': 
                all_eligible['cont_vuln'] = self.contact_vuln(sim)

            if name == 'cont_sx': 
                all_eligible['cont_sx'] = self.contact_sx(sim, mode='all')
            if name == 'cont_cs_sx': 
                all_eligible['cont_cs_sx'] = self.contact_sx(sim, mode='specific')
            if name == 'cont_ncs_sx':
                all_eligible['cont_ncs_sx'] = self.contact_sx(sim, mode='nonspecific')
            
            if name == 'sx': 
                all_eligible['sx'] = self.sx(sim, mode='all')
            if name == 'cs_sx': 
                all_eligible['cs_sx'] = self.sx(sim, mode='specific')
            if name == 'ncs_sx': 
                all_eligible['ncs_sx'] = self.sx(sim, mode='nonspecific')

            if name == 'pos_RAT': 
                all_eligible['pos_RAT'] = self.RAT_testers(sim, mode='pos_RAT')
            if name == 'pos_RAT_asx': 
                all_eligible['pos_RAT_asx'] = self.RAT_testers(sim, mode='pos_RAT_asx')
            if name == 'neg_RAT_sx': 
                all_eligible['neg_RAT_sx'] = self.RAT_testers(sim, mode='neg_RAT_sx')

            if name == 'travel': 
                all_eligible['travel'] = self.travel(sim)
            if name == 'work': 
                all_eligible['work'] = self.work(sim, surv)
            
        return all_eligible


    def contact_conf(self, sim, surv, d=3): 
        '''
        Get people who were in contact with a confirmed case in the past d days. Select a subset who seeks a test. 
            * Note that this criterion inherently requires people to have tested positive previously. If this is the only eligibility criterion 
              used to determine test seekers and testers, then no one will end up testing. 

        Args: 
            sim  (simulation)  : Simulation object 
            surv (surveillance): Surveillance object
            d    (int)         : Used to get people who were confirmed positive in past d days

        Returns: 
            seekers (set)       : People who satisfy the criterion and will seek a test
        '''
        cur = sim.t

        # Get people who have tested positive previously (and are not dead) with the current surveillance system
        pos_inds = cvu.true(np.logical_and(surv.date_positive >= 0, ~sim.people.dead))
        conf_uid = pos_inds[np.logical_and(cur - surv.date_positive[pos_inds] <= d, cur - surv.date_positive[pos_inds] >= 0)]

        # Get people who have tested positive previously (and are not dead) with passive PCR surveillance
        if surv.hybrid and (surv.system == 'PassiveRAT' or surv.system == 'ActiveRAT'): 
            PCR_pos_inds = cvu.true(np.logical_and(surv.PCR.date_positive >= 0, ~sim.people.dead))
            PCR_conf_uid = PCR_pos_inds[np.logical_and(cur - surv.PCR.date_positive[PCR_pos_inds] <= d, cur - surv.PCR.date_positive[PCR_pos_inds] >= 0)]
            conf_uid = np.union1d(conf_uid, PCR_conf_uid)

        # For each layer, get contacts of these people - i.e., people who satisfy the criterion
        contacts = set()  # One person may have contacted several confirmed positives, which may increase probability of test-seeking. We are currently not accounting for this.
        for layer in sim.people.contacts.keys(): 
            contacts.update(sim.people.contacts[layer].find_contacts(conf_uid))  # Adds empty list if agent i does not exist in the layer, or if layer is empty. Sped up by Numba.

        return contacts


    def contact_vuln(self, sim):  
        '''
        Get people who are in contact with vulnerable people in the household. Vulnerability defined in terms of odds ratio of susceptibility. 
        '''
        # Get people who are vulnerable (and not dead) 
        vuln_uid = cvu.true(np.logical_and(sim.people.rel_sus > 1.0, ~sim.people.dead))  # See parameters.py for assignment of susceptibility odds ratios by age bracket

        # Get their household contacts. It may be more realistic to only get some subset of the household contacts.
        contacts = set()
        contacts.update(sim.people.contacts['h'].find_contacts(vuln_uid))  # Sped up by Numba
        
        return contacts


    def contact_sx(self, sim, d=3, mode='all'): 
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

        #TODO: Look for an elegant implementation of "contact with person in the last d days"
        # sx_COVID = set(cvu.true(sim.people.symptomatic))  # Anyone who is currently symptomatic (and not dead)
        # # sx_COVID.update(cvu.true(np.bitwise_and(sim.people.date_recovered + d < sim.t, ~np.isnan(sim.people.date_symptomatic))))  # Anyone who recovered, but was symptomatic in last d days

        # sx_ILI = set(cvu.true(sim.people.symptomatic_ILI))  # Do the same with ILI 
        # # sx_COVID.update(cvu.true(np.bitwise_and(sim.people.date_recovered_ILI + d < sim.t, ~np.isnan(sim.people.date_symptomatic_ILI))))

        # all_sx = sx_COVID.union(sx_ILI)

        all_sx = set(cvu.true(np.logical_or(sim.people.symptomatic, sim.people.symptomatic_ILI)))  # Anyone who is currently symptomatic with COVID or ILI

        if mode == 'all': 
            all_contacts = set() 
            for layer in sim.people.contacts.keys(): 
                all_contacts.update(sim.people.contacts[layer].find_contacts(np.array(list(all_sx))))  # Sped up by Numba
            return all_contacts

        elif mode == 'specific':
            return self.stratify_contacts_sx(sim, all_sx)[0]

        elif mode == 'nonspecific': 
            return self.stratify_contacts_sx(sim, all_sx)[1]


    def sx(self, sim, mode='all'): 
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

        sx_COVID = set(cvu.true(sim.people.symptomatic))
        sx_ILI = set(cvu.true(sim.people.symptomatic_ILI))
        all_sx = sx_COVID.union(sx_ILI)
        cs_sx, ncs_sx = self.stratify_sx(sim, all_sx)

        if mode == 'all': 
            return all_sx

        elif mode == 'specific': 
            return cs_sx

        elif mode == 'nonspecific': 
            return ncs_sx


    def RAT_testers(self, sim, mode='pos_RAT', d=3): 
        '''
        Interface between passive PCR and active RAT surveillance in hybrid model. Get people who either: 
            1) Received positive RAT in past d days 
            2) Received positive RAT in past d days and is asx OR received negative RAT in past d days and is sx
        '''
        cur = sim.t
        date_positive = self.RAT.date_positive
        date_negative = self.RAT.date_negative

        # Get people who had a positive/negative RAT previously (and are not dead) 
        pos_inds = cvu.true(np.logical_and(date_positive >= 0, ~sim.people.dead))
        neg_inds = cvu.true(np.logical_and(date_negative >= 0, ~sim.people.dead))

        # Get people who were confirmed positive/negative by a RAT test in the last d days
        pos_d_inds = pos_inds[np.logical_and(cur - date_positive[pos_inds] <= d, cur - date_positive[pos_inds] >= 0)]  # 2nd part of "and" statement is needed since reporting date may be in the future
        neg_d_inds = neg_inds[np.logical_and(cur - date_negative[neg_inds] <= d, cur - date_negative[neg_inds] >= 0)]

        if mode == 'pos_RAT': 
            # Return set of individuals who got a positive RAT test in past d days
            return set(pos_d_inds)

        elif mode == 'pos_RAT_asx':  
            # Return set of individuals who got a positive RAT test in past d days and was asymptomatic at time of testing. Two cases here:
                # Case 1: They became symptomatic within these d days
                # Case 2: They remained asymptomatic within these d days
                # Intuitively, a confirmatory test would be warranted in both cases
            return set(pos_d_inds[self.RAT.pos_was_asx[pos_d_inds] == 2])
            
            
        elif mode == 'neg_RAT_sx': 
            # Return set of individuals who got a negative RAT test in past d days and was symptomatic at time of testing
            if sim.t == 0: 
                setattr(self.RAT, 'neg_was_sx_hist', [])
            self.RAT.neg_was_sx_hist.append(sum(self.RAT.neg_was_sx == 1))
            return set(neg_d_inds[self.RAT.neg_was_sx[neg_d_inds] == 1])
            


    
    # TODO: Figure out how to implement this when the population is static. Merge with "other".
    def travel(self, sim): 
        '''
        Get people who have to do mandatory travel testing. Currently not modelled since the population remains static. 
        '''
        return set()

    
    #TODO: Stratify workplaces by industry, though we should ask if this is too detailed
    def work(self, sim, surv, p=0.10):
        '''
        Get people who have to do mandatory workplace testing. Randomly select some number of workplaces. Every week, these people will seek a test. 
        To consider: Do we want them to have priority access to testing resources? 

        Args: 
            p    (float) : Proportion of workplaces to undergo mandatory testing
        '''
        if sim.t == 0:
            # Identify which workplaces are chosen to do mandatory testing, and save the workers
            test_wpid = np.random.choice(len(surv.pop.workplaces), int(surv.pop.n_workplaces * p), replace=False)

            worker_uids = np.array([])
            for wpid in test_wpid:
                worker_uids = np.union1d(worker_uids, surv.pop.workplaces[wpid]['member_uids']).astype(int)
            
            worker_test_dates = np.random.choice(7, len(worker_uids))

            setattr(surv, 'wuid', worker_uids)  # Set of all workers
            setattr(surv, 'wuid_td', worker_test_dates)  # Days of the week they are assigned to test


        # Remove anyone who died
        surv.wuid_td = surv.wuid_td[~sim.people.dead[surv.wuid]]  # Order matters
        surv.wuid = surv.wuid[~sim.people.dead[surv.wuid]]

        # Same workplaces, and thus workers, test each time
        return set(list(surv.wuid[surv.wuid_td == sim.t % 7]))
        

    def stratify_sx(self, sim, sx_uid): 
        '''
        Helper function. Stratify those who have symptoms by whether their symptoms are specific or nonspecific. 

        Args: 
            sim      (sim)      : Simulation object
            sx_uid  (np.array) : UID of symptomatic people
        '''
        # Define specific and nonspecific symptoms. Specific symptoms have higher associated probability of test seeking. 
        cs_sx = set(['fever', 'cough', 'sore_throat', 'short_breath', 'loss_taste_smell', 'muscle_joint_pain', 'chills_shaking'])
        ncs_sx = set(['tight_chest', 'fatigue', 'nausea', 'runny nose'])

        cs_uid = set() 
        ncs_uid = set() 

        for uid in sx_uid: 
            symps = sim.people.symptom_types[uid]
            has_cs = ('fever' in symps) and (('cough' in symps) or ('sore_throat' in symps))

            if has_cs:
                cs_uid.add(uid)
            elif not(has_cs): 
                ncs_uid.add(uid)

        return cs_uid, ncs_uid


    def stratify_contacts_sx(self, sim, sx_uid): 
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
        cs_uid, ncs_uid = self.stratify_sx(sim, sx_uid)

        cs_contacts = set()
        ncs_contacts = set()

        for layer in sim.people.contacts.keys(): 
            cs_contacts.update(sim.people.contacts[layer].find_contacts(np.array(list(cs_uid))))
            ncs_contacts.update(sim.people.contacts[layer].find_contacts(np.array(list(ncs_uid))))

        return cs_contacts, ncs_contacts

    

# sx_COVID and sx_ILI could go under an analysis function
    # def sx_COVID(self, sim):
    #     '''
    #     Get people who are symptomatic and have COVID-19 (this is just people who are symptomatic in Covasim). We only look at people who are mildly symptomatic. 
    #     This is probably not a realistic criterion, since people won't know if they have COVID or ILI
    #     '''
    #     # sx_c_uid = set(cvu.true(np.bitwise_and(sim.people.symptomatic, ~sim.people.severe)))  # Critical people are also considered severe, so this takes care of critical people
        
    #     sx_c_uid = set(cvu.true(sim.people.symptomatic))  # ALL symptoms, whether mild or severe

    #     return sx_c_uid


    # def sx_ILI(self, sim): 
    #     '''
    #     Get people who have non-COVID ILI symptoms
    #     This is probably not a realistic criterion, since people won't know if they have COVID or ILI
    #     '''
    #     sx_nc_uid = set(cvu.true(sim.people.symptomatic_ILI))
    #     return sx_nc_uid


##############################
######## Get seekers #########
##############################
class TestSeek: 
    '''
    Out of eligible test seekers, determine who seeks a test. 
        1) Define test seeking probabilities associated with each criterion
        2) Probabilistically determine who seeks a test

    Args: 
        seek_probs      (dict) : Defines the test seeking probabilities
    '''
    def __init__(self, eligible_dict, seek_probs=None): 

        if seek_probs is None: 
            self.seek_probs = {
                'cont_conf': 0.8, 
                'cont_vuln': 0.2,
                'cont_sx': None, 
                'cont_cs_sx': 0.3,
                'cont_ncs_sx': 0.1,
                'sx': None, 
                'cs_sx': 0.8,
                'ncs_sx': 0.4,
                'pos_RAT': 0.5,
                'pos_RAT_asx': 0.8,
                'neg_RAT_sx': 0.8,
                'travel': 1.0,
                'work': 1.0
            }
        else: 
            self.seek_probs = seek_probs

        self.all_seekers = self.get_all_seekers(eligible_dict)

    
    def get_seekers(self, candidates, criterion): 
        ''' 
        Takes as input a set of potential test-seekers. Determines if they actually will seek a test, probabilistically. 

        Args: 
            candidates  (set): Set of people who satisfy a test-seeking criterion
            criterion   (str): Name of the criterion, used to key self.criteria (dict)

        Returns: 
            seekers     (set): Set of people who will seek a test
        '''
        seekers = set()
        prob = self.seek_probs[criterion]

        if not(prob is None):  # cs_sx and ncs_sx are a partition of sx, and we already calculate who seeks a test among cs_sx/ncs_sx. Hence we don't calculate seekers for sx, as there would be double-counting. 
            for c in candidates: 
                if np.random.binomial(1, prob): 
                    seekers.add(c) 

        return seekers

    
    def get_all_seekers(self, eligible_dict):
        '''
        Look at who satisfies the test criteria, and determine which of these people seek a test

        Args: 
            sim              (sim)           : Simulation object 
            surv             (surveillance)  : Surveillance system
            eligible_dict    (dict)          : Dictionary containing sets of people who meet each criterion

        Returns: 
            test_seekers  (dict)    : Dictionary of test seekers for each criterion
        '''
        test_seekers = dict()

        for criterion in eligible_dict.keys():
                test_seekers[criterion] = self.get_seekers(eligible_dict[criterion], criterion)

        return test_seekers



###############################
####### Allocate tests ########
###############################

class TestAllocate: 
    '''
    Among test seekers, this class helps identify who receives a test. For now, the implementation is basic, and randomly selects among test seekers. 
    '''
    def __init__(self, test_capacity, eligible_dict, seekers_dict, criteria, seek_probs, allocation_type='proportional'): 
        # Placeholder for now, until we come up with a more complex algorithm
        if allocation_type == 'basic': 
            self.all_testers = self.basic_PCR_allocation(test_capacity, seekers_dict)
        elif allocation_type == 'proportional': 
            self.all_testers = self.proportional_PCR_allocation(test_capacity, eligible_dict, seekers_dict, criteria, seek_probs)
        return


    def basic_PCR_allocation(self, test_capacity, seekers_dict): 
        '''
        Each day, get people who are seeking a PCR test. Randomly determine who gets a test. There is no temporal dimension, people do not stay in any sort of queue if 
        they do not get a test today.  
        Conduct analysis if desired.

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
            testers = np.array(list(test_seekers))
        else: 
            testers = np.random.choice(list(test_seekers), test_capacity, replace=False)  # Insufficient capacity, randomly choose a lucky few test-seekers

        return testers


    def proportional_PCR_allocation(self, test_capacity, eligible_dict, seekers_dict, criteria, seek_probs): 
        '''
        Compute the power set of criteria. For each power set element (a set of criteria), get the number of people who precisely met these criteria 
        at the eligibility assessment, and divide by the total number of people in the population. This proportion (at most) of tests will be allocated to
        those belonging to the power set element. Note that the power set contains the empty set. This corresponds to people who did not satisfy any criteria, 
        which is untrue of any test seeker. 

        Recall if S = len(X) then power_set(X) has length 2^S

        Args: 
            sim    (simulation) : Simulation object
            test_capacity   (int) : Number of tests to distribute on this day 
            eligibile_dict  (dict) : Keys are criteria, values are people who meet the criteria
            seekers_dict    (dict) : Keys are criteria, values are test-seekers who meet the criteria
            criteria        (list) : 

        '''
        valid_criteria = set([criterion for criterion in criteria if not(seek_probs[criterion] is None)])  # Do not account for cont_conf and sx sets, since they have subsets that are already accounted for
        
        criteria_powset = self.compute_criteria_powset(valid_criteria)

        # # Get eligibile people and the powerset elements they correspond to
        # eligible_crit_subset = self.compute_members(valid_criteria, criteria_powset, eligible_dict)  # This was used for allocation proportional to eligibility

        # Get test seekers and the powerset elements they correspond to
        seeker_crit_subset = self.compute_members(valid_criteria, criteria_powset, seekers_dict)

        # Compute proportions and allocate
        # test_proportions = self.compute_proportions(eligible_crit_subset)
        test_proportions = self.compute_proportions(seeker_crit_subset)

        testers = set()
        for sname in test_proportions.keys(): 
            subset_seekers = np.array(list(seeker_crit_subset[sname]))  # Test seekers who meet the criteria described by subset
            subset_n_tests = int(test_proportions[sname] * test_capacity)
            # subset_n_tests = math.ceil(test_proportions[sname] * test_capacity)

            if subset_n_tests <= len(subset_seekers):
                subset_testers = set(np.random.choice(subset_seekers, subset_n_tests, replace=False))
            else: 
                print(sname, "does not have enough seekers!")
                print("They were allocated", subset_n_tests, "tests")
                print("But they only have", len(subset_seekers), "seekers")
                print(sname, "criteria consists of", criteria_powset[sname])
                print()
                subset_testers = subset_seekers

            testers.update(subset_testers)

        # Redistribute extras by randomly sampling subsets according to proportions
        extra_tests = test_capacity - len(testers)
        while extra_tests > 0: 
            probs = [test_proportions[sname] for sname in test_proportions.keys()]
            extra_sname = np.random.choice(list(test_proportions.keys()), 1, p=probs)[0]
            not_already_tested = seeker_crit_subset[extra_sname].difference(testers)

            if len(not_already_tested) > 0: 
                testers.add(np.random.choice(list(not_already_tested)))
                extra_tests -= 1
            
        return np.array(list(testers))

    
    def compute_members(self, valid_criteria, criteria_powset, criteria_dict): 
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
            complement = valid_criteria.difference(crit_subset_copy)  # Criteria we don't want seekers to meet
    
            members = criteria_dict[crit_subset_copy.pop()]  # Choose any criterion to start with, pop() can only happen after line above 

            # Get people who meet at least the prescribed subset criteria 
            if len(crit_subset_copy) > 0:  # May be zero if it only contained one element initially, in that case then we already have all members 
                for criterion in crit_subset_copy: 
                    members = members.intersection(criteria_dict[criterion])

            # Remove anyone who meets other criteria
            if len(complement) > 0:  # May be zero when the crit_subset contains all criteria, in that case then there is nothing to remove
                for criterion in complement: 
                    members = members - members.intersection(criteria_dict[criterion])

            meets_crit_subset[subset_key] = members  # Meets all and only the subset criteria

        return meets_crit_subset

    
    def compute_proportions(self, crit_subset_dict): 
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
            proportions[sname] = len(crit_subset_dict[sname]) / n_eligible

        return proportions
        
    
    def compute_criteria_powset(self, valid_criteria): 
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


class RAT_allocation:  # Merge this into TestAllocate
    '''
    We currently do not define test eligibility and test seeking modules for RAT testing. We directly allocate tests according to an active surveillance policy. 

    Args: 
        asurv       (RAT) : Initialized RAT surveillance class
    '''
    def __init__(self, asurv, sim):
        if asurv.policy == 'cross_sectional': 
            self.cur_testers = self.cross_allocation(asurv, sim)
        elif asurv.policy == 'longitudinal': 
            self.cur_testers = self.long_allocation(asurv, sim)
        else: 
            raise RuntimeError('Invalid policy supplied, options are: "cross_allocation" or "longitudinal"')

        fname = '/Users/Ritchie/Desktop/cv-test/covasim/covasim/data/age_distribution_reports.JSON'
        with open(fname, 'r') as f: 
            self.age_distr = json.load(f)
            
        return 


    def cross_allocation(self, asurv, sim): 
        '''
        At some frequency, randomly choose a subset of people from the targeted populations.

        '''
        cur = sim.t

        if cur == 0: 
            setattr(asurv, 'all_testers', None)
            setattr(asurv, 'test_days', None)

        if cur % asurv.frequency == 0:  # E.g. if f=7, then we resample after 7 days have passed (0, 1, 2, 3, 4, 5, 6)
            asurv.all_testers = self.sample_targets(asurv, sim)  # Resample according to defined frequency
            asurv.test_days = self.calculate_test_days(cur, asurv.frequency, asurv.all_testers)
        
        return asurv.all_testers['all'][asurv.test_days['all'] == cur]  #TODO: Figure out how to test specific/combinations of population layers. Could return a dictionary of those who are testing today. 


    def long_allocation(self, asurv, sim): 
        '''
        On the first day, randomly choose a subset of people from the targeted population(s). At some frequency, re-test each individual. 
        Example: 
            f = 7 days. On day 0, individuals A, B, C are tested. On the 7th day after their initial test (day 7), they get tested again. 
        '''
        cur = sim.t

        if cur == 0: 
            setattr(asurv, 'all_testers', None)
            setattr(asurv, 'test_days', None)

            asurv.all_testers = self.sample_targets(asurv, sim)  # These individuals will be followed throughout the outbreak for testing
            asurv.test_days = self.calculate_test_days(cur, asurv.frequency, asurv.all_testers)  # Initialize the first day of testing for each selected individual

        return asurv.all_testers['all'][(cur - asurv.test_days['all']) % asurv.frequency == 0]


    def get_targets(self, asurv, sim): 
        '''
        Get target populations. Eventually this may expand into its own class (akin to TestEligible) if we want more detailed targetting (e.g. classroom level)

        NOTE: For now, testing is only designed around targetting one population layer. Not multiple, even though the infrastructure is partly there. 
        '''
        target_layers = asurv.population 

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


    def sample_targets(self, asurv, sim):
        '''
        Randomly select some proportion of each targeted population layer. 

        NOTE: To allow multiple selected population layers, we need to define what happens when a person is selected in multiple layers
        '''
        target_uid = self.get_targets(asurv, sim)
        p = asurv.participation

        sampled_uid = dict()
        if not(asurv.age_report is None):  # Sample according to custom defined distributions NOTE: This feature is currently only developed for the case that we are ONLY targetting 'all'
            age_group_members = self.get_age_group_members(asurv, sim, target_uid)
            age_distribution = asurv.all_age_reports[asurv.age_report]['distribution']

            # if sim.t == 0: 
            #     setattr(asurv, 'age_distribution_history', dict())
            #     for interval in age_group_members.keys(): 
            #         asurv.age_distribution_history[interval] = []

            for l in target_uid.keys():  # 'all' is the only option
                sampled_uid[l] = set()
                for age_group, interval in zip(age_distribution, age_group_members.keys()):
                    n_tests = int(age_group[2] * p[l] * len(target_uid[l]))  # Age prevalence * participation * population size
                    sampled_uid[l].update(set(list(np.random.choice(age_group_members[interval], n_tests, replace=False))))
                    # asurv.age_distribution_history[interval].append(n_tests)
                
                sampled_uid[l] = np.array(list(sampled_uid[l]))

        else:  # Pure random sampling
            for l in target_uid.keys(): 
                sampled_uid[l] = np.random.choice(target_uid[l], int(len(target_uid[l]) * p[l]), replace=False)

        return sampled_uid

    
    def calculate_test_days(self, t, f, sampled_uid):
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


    def get_age_group_members(self, asurv, sim, target_uid):
        '''
        Args: 
            asurv         (surveillance) : RAT surveillance object 
            target_uid    (dict)         : Dictionary of targetted individuals for different population layers
        '''
        age_distribution = asurv.all_age_reports[asurv.age_report]['distribution']

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
                    
                    



###############################
##### Allocation analysis #####
###############################

# General analysis

def TestEligibleGroups(sim, surv, criteria, all_eligible): 
    '''
    Get the number of people who meet each criteria each day. 
    
    Args: 
        sim        (sim)          : Simulation object 
        surv       (surveillance) : Surveillance system
        criteria   (list)         : Eligibility criteria for testing
        all_eligible  (dict)    : Values are sets containing people who meet each criterion 
    '''
    if sim.t == 0: 
        setattr(surv, 'crit_groups', dict())
        for criterion in criteria: 
            surv.crit_groups[criterion] = []

    for criterion in surv.crit_groups.keys():
        surv.crit_groups[criterion].append(len(all_eligible[criterion])) 


def TestSeekGroups(sim, surv, criteria, all_seekers): 
    '''
    For each testing criterion, get the number of test seekers belonging to the group. 

    Args: 
        sim             (sim) : Simualtion object 
        surv   (surveillance) : Surveillance system 
        criteria       (list) : Test criteria 
        all_seekers    (dict) : Values are sets containing people who seek a test for each criterion
    '''
    if sim.t == 0: 
        setattr(surv, 'seek_groups', dict())
        for criterion in criteria: 
            surv.seek_groups[criterion] = []

    # Get all test seekers at the current state of the simulation
    for criterion in surv.seek_groups.keys(): 
        surv.seek_groups[criterion].append(len(all_seekers[criterion]))


def TestAllocateGroups(sim, surv, testers, criteria, all_eligible): 
    '''
    Among those that received a PCR test, identify which criteria each person satisfied - this is their "reason" for testing. 

    Args: 
        sim             (sim)           : Simulation object 
        surv            (surveillance)  : Surveillance object 
        testers         (np.array)      : UID of each person who got tested 
        criteria        (list)          : Test criteria
        all_eligible  (dict)          : Values are sets containing people who meet each criterion
    '''
    if sim.t == 0: 
        setattr(surv, 'alloc_groups', dict())
        for criterion in criteria: 
            surv.alloc_groups[criterion] = []
    
    testers = set(testers)

    for criterion in surv.alloc_groups.keys():
        surv.alloc_groups[criterion].append(len(testers.intersection(all_eligible[criterion])))

    return





