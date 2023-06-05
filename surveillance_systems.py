import numpy as np
import sciris as sc
import json

import infection.surveillance_allocation as sa
import infection.surveillance_testing as st

from infection.surveillance_utils import *
from scipy import stats

# NBU = Not Being Used


class Surveillance: 
    '''
    Base class for implementing all surveillance objects. 

    Args: 
        initial_unit   (any): Initial unit to fill surveillance arrays with
    '''

    def __init__(self, initial_unit):
        self.system = self.__class__.__name__ # Name of the type of data to be collected by surveillance strategy 
        self.initialized = False # Whether or not it has been initialized
        self.finalized = False # Whether or not it has been finalized 
        self.surv_history = []  # Contain surveillance history, NBU
        self.initial_unit = initial_unit  # Base type to fill surveillance array, can be bool (e.g. False), int (e.g. 0), objects, etc. 
        return


    def __call__(self, *args, **kwargs):
        # Makes Surveillance(sim) equivalent to Surveillance.apply(sim)
        if not self.initialized:  # pragma: no cover
            errormsg = f'Surveillance (label={self.system}, {type(self)}) has not been initialized'
            raise RuntimeError(errormsg)
        return self.apply(*args, **kwargs)


    def initialize(self, sim): 
        '''
        Initialize the Surveillance object
        '''
        sim.people.unlock()  # Make sure people is unlocked to allow adding keys and editing values 
        sim.people[self.system] = np.full(sim.people.pars['pop_size'], self.initial_unit) # Initialize ground truth agent data
        self.initialized = True
        return


    def finalize(self, sim):
        '''
        Finalize the Surveillance object 
        '''
        sim.people.lock()  # Prevent further changes
        if self.finalized: # pragma: no cover
            raise RuntimeError('Surveillance already finalized')  # Raise an error because finalizing multiple times has a high probability of producing incorrect results e.g. applying rescale factors twice
        self.finalized = True
        return


    def apply(self, sim): 
        '''
        Function to be implemented by user.
        '''
        raise NotImplementedError


class PassivePCR(Surveillance):
    '''
    This PCR class is intended to recapitulate PCR allocation more realistically than the symptom_PCR class. 
    Tests funnel into several categories: 
        1) Personal testing 
            - If someone was in contact with a confirmed case
            - If someone 
        2) Mandatory testing 
        3) Unused
    
    '''
    def __init__(self,
                initial_unit,
                pop,
                sensitivity=0.85, 
                specificity=1.00, 
                criteria=None,
                seek_probs=None,
                capacity=None,
                pos_retest_delay=90,  # NBU
                neg_retest_delay=14,  # NBU
                latency_shape=2,  
                latency_scale=1,
                static_latency=None,
                RAT_ind=None,
                hybrid=False,
                analysis=True):
        '''
        Args: 
            pop                  (sp.Pop)    : SynthPops Pop object, contains detailed workplace/school information
            sensitivity          (float)     : PCR test sensitivity
            specificity          (float)     : PCR test specificity 
            criteria             (list)      : List of criteria that qualify an individual to seek a test
            seek_probs           (dict)      : Dictionary of criteria associated with test-seeking probabilities
            swab_delay           (int)       : Number of days agent gets tested after becoming symptomatic. Default gamma distributed. 
            result_delay         (int)       : Number of days before test results are reported.
            capacity             (int)       : Test capacity per day. Defaults to UK daily test capacity if not provided. 
            pos_retest_delay     (int)       : Number of days before a person who tests positive can retest. Independent of whether they really are positive. 
            neg_retest_delay     (int)       : Number of days before a person who tests negative can retest. Independent of whether they really are negative. 
            capacity             (any)       : Can be int or list
            RAT_ind              (int)       : Specify which RAT surveillance system to interface with if supplied 
            analysis             (bool)      : Whether to call analysis functions or not
        '''
        super().__init__(initial_unit)
        self.sensitivity = sensitivity
        self.specificity = specificity 
        self.criteria = criteria
        self.seek_probs = seek_probs
        self.pos_retest_delay = pos_retest_delay  # TODO: Research this parameter
        self.neg_retest_delay = neg_retest_delay  # TODO: Research this parameter
        self.latency_shape = latency_shape
        self.latency_scale = latency_scale
        self.static_latency = static_latency
        self.capacity = capacity
        self.pop = pop
        self.RAT_ind = RAT_ind
        self.analysis = analysis
        return


    def initialize(self, sim):
        '''
        Args: 
            sim         (sim): Simulation object. Do not initialize simulation with start/end day, use n_day to avoid potential error)
        '''
        super().initialize(sim)  # Unlock people object and set initialized as true 

        # Initialize count trackers
        # self.surv_history = np.full(sim.pars['n_days'] + 1, self.initial_unit)  # NBU, turn on if using old surveillance systems
        self.pos_history = np.full(sim.pars['n_days'] + 1, self.initial_unit)  # The simulation always goes on for one extra day, hence +1
        self.neg_history = np.full(sim.pars['n_days'] + 1, self.initial_unit)

        # Initialize date trackers
        self.date_positive = np.full(sim.pars['pop_size'], -1)
        self.date_negative = np.full(sim.pars['pop_size'], -1)

        # Initialize overflow trackers 
        self.pos_overflow = 0
        self.neg_overflow = 0

        # Initialize other trackers
        self.retest_tracker = np.full(sim.pars['pop_size'], 0)  # Track whether an individual is allowed to retest 

        # Initialize interface with RAT system
        if not(self.RAT_ind is None): 
            # Two-way interface with active RAT surveillance system 
            self.RAT = sim.pars['surveillance'][self.RAT_ind] 
            self.RAT.PCR = self

            # Turn on lights
            self.RAT.hybrid = True
            self.hybrid = True

         # WARNING: If testing capacity is too high, eventually everyone will have to wait for retest, and code will break
        if self.capacity is None: 
            self.capacity = process_UK_PCR_capacity()[0] * (sim.pars['pop_size'] / 1e5) # UK capacity is calculated per 100k people

        return


    def finalize(self, sim):
        '''
        Do any finalizing work here (e.g. resetting to save memory)
        '''
        super().finalize(sim)
        return


    def apply(self, sim, analysis=True): 
        cur = sim.t

        # Assign capacity
        if sc.isnumber(self.capacity):
            tests_left = self.capacity
        else: 
            tests_left = int(self.capacity[cur])  # Number of individuals to test (assuming the capacity estimate is correct, and all tests are used)

        # Determine who is eligible to receive a test 
        if not(self.RAT_ind is None):
            eligible_module = sa.TestEligible(self, sim, criteria=self.criteria, RAT=self.RAT)
        else: 
            eligible_module = sa.TestEligible(self, sim, criteria=self.criteria)
        eligible_dict = eligible_module.all_eligible

        # Determine who seeks a test
        seeker_module = sa.TestSeek(eligible_dict, seek_probs=self.seek_probs)
        seekers_dict = seeker_module.all_seekers
        
        # Determine who is allocated a test
        tester_module = sa.TestAllocate(tests_left, eligible_dict, seekers_dict, eligible_module.criteria, seeker_module.seek_probs, allocation_type='proportional')
        test_uid = tester_module.all_testers  # Numpy array 

        # Administer the tests
        test_results = st.basic_PCR_test(self, sim, test_uid, self.sensitivity, self.specificity)

        # Do some analyses
        if self.analysis: 
            # General analysis 
            # log_symptomatic(self, sim, eligible_dict)

            # Allocation analysis
            sa.TestEligibleGroups(sim, self, eligible_module.criteria, eligible_dict)
            sa.TestSeekGroups(sim, self, eligible_module.criteria, seekers_dict)
            sa.TestAllocateGroups(sim, self, test_uid, eligible_module.criteria, eligible_dict)

        if len(test_results) == 0:  # Nothing need to be updated if no tests were conducted
            return
        
        # Record the positive results with a result latency term (swab delay is not included in this term, as we are not delaying test administration)
        # Generate test latencies
        if not(self.static_latency is None): 
            test_latencies = np.full(len(test_results), self.static_latency)
            unique_latencies = [self.static_latency]
        else: 
            test_latencies = np.trunc(stats.gamma.rvs(a=self.latency_shape, scale=self.latency_scale, size=len(test_results)))  # Associate each test result with a latency
            unique_latencies = np.unique(test_latencies).astype(np.intc)

        for l in unique_latencies:  # This should be a short loop
            if cur + l <= sim.pars['n_days']: 

                # Track case counts
                self.pos_history[cur + l] += sum(test_results[test_latencies == l])
                self.neg_history[cur + l] += sum(~test_results[test_latencies == l])  # Evaluates ~test_results before boolean indexing

                # Track dates
                self.date_positive[test_uid[np.logical_and(test_results, test_latencies == l)]] = cur + l  # Pick people who returned positive result and was assigned latency l 
                self.date_negative[test_uid[np.logical_and(~test_results, test_latencies == l)]] = cur + l  # Pick people who returned negative result and was assigned latency l 

            else:  # Pool tests that are recorded after the hypothetical outbreak ends
                self.pos_overflow += sum(test_results[test_latencies == l])
                self.neg_overflow += sum(~test_results[test_latencies == l])

        return 

# TODO: Add a new eligibility criteria that allows people to retest if they got a positive RAT test -- circular 
class PassiveRAT(Surveillance):

    def __init__(self,
                initial_unit,
                pop,
                sensitivity=0.70, 
                specificity=0.99, 
                criteria=None,
                seek_probs=None,
                capacity=0.50,
                pos_retest_delay=90,  # NBU
                neg_retest_delay=14,  # NBU
                latency_shape=2,  
                latency_scale=1,
                static_latency=None,
                hybrid=False,
                analysis=True):
        '''
        Args: 
            pop                  (sp.Pop)    : SynthPops Pop object, contains detailed workplace/school information
            sensitivity          (float)     : PCR test sensitivity
            specificity          (float)     : PCR test specificity 
            criteria             (list)      : List of criteria that qualify an individual to seek a test
            seek_probs           (dict)      : Dictionary of criteria associated with test-seeking probabilities
            swab_delay           (int)       : Number of days agent gets tested after becoming symptomatic. Default gamma distributed. 
            result_delay         (int)       : Number of days before test results are reported.
            capacity             (int)       : Test capacity per day. Defaults to UK daily test capacity if not provided. 
            pos_retest_delay     (int)       : Number of days before a person who tests positive can retest. Independent of whether they really are positive. 
            neg_retest_delay     (int)       : Number of days before a person who tests negative can retest. Independent of whether they really are negative. 
            capacity             (any)       : Proportion of test seeking individuals who can be accomodated
            RAT_ind              (int)       : Specify which RAT surveillance system to interface with if supplied 
            analysis             (bool)      : Whether to call analysis functions or not
        '''
        super().__init__(initial_unit)
        self.sensitivity = sensitivity
        self.specificity = specificity 
        self.criteria = criteria
        self.seek_probs = seek_probs
        self.pos_retest_delay = pos_retest_delay  # TODO: Research this parameter
        self.neg_retest_delay = neg_retest_delay  # TODO: Research this parameter
        self.latency_shape = latency_shape
        self.latency_scale = latency_scale
        self.static_latency = static_latency
        self.capacity = capacity
        self.pop = pop
        self.analysis = analysis
        return


    def initialize(self, sim):
        '''
        Args: 
            sim         (sim): Simulation object. Do not initialize simulation with start/end day, use n_day to avoid potential error)
        '''
        super().initialize(sim)  # Unlock people object and set initialized as true 

        # Initialize count trackers
        self.pos_history = np.full(sim.pars['n_days'] + 1, self.initial_unit)  # The simulation always goes on for one extra day, hence +1
        self.neg_history = np.full(sim.pars['n_days'] + 1, self.initial_unit)

        # Initialize date trackers
        self.date_positive = np.full(sim.pars['pop_size'], -1)
        self.date_negative = np.full(sim.pars['pop_size'], -1)

        # Initialize overflow trackers 
        self.pos_overflow = 0
        self.neg_overflow = 0

        # Initialize other trackers
        self.retest_tracker = np.full(sim.pars['pop_size'], 0)  # Track whether an individual is allowed to retest 

        return


    def finalize(self, sim):
        '''
        Do any finalizing work here (e.g. resetting to save memory)
        '''
        super().finalize(sim)
        return


    def apply(self, sim, analysis=True): 
        cur = sim.t

        # Determine who is eligible to test
        eligible_module = sa.TestEligible(self, sim, criteria=self.criteria)
        eligible_dict = eligible_module.all_eligible

        # Determine who seeks a test
        seeker_module = sa.TestSeek(eligible_dict, seek_probs=self.seek_probs)
        seekers_dict = seeker_module.all_seekers

        # Compute capacity
        unique_seekers = set()
        for value in seekers_dict.values(): 
            unique_seekers.update(value)
        tests_left = int(self.capacity * len(unique_seekers))
        
        # Determine who is allocated a test
        tester_module = sa.TestAllocate(tests_left, eligible_dict, seekers_dict, eligible_module.criteria, seeker_module.seek_probs, allocation_type='proportional')
        test_uid = tester_module.all_testers  # Numpy array 

        if cur == 0: 
            setattr(self, 'test_uid_length', [])
        
        self.test_uid_length.append(len(test_uid))

        # Administer the tests
        test_results = st.basic_RAT_test(self, sim, test_uid, self.sensitivity, self.specificity)

        # Do some analyses
        if self.analysis: 
            # General analysis 
            # log_symptomatic(self, sim, eligible_dict)

            # Allocation analysis
            sa.TestEligibleGroups(sim, self, eligible_module.criteria, eligible_dict)
            sa.TestSeekGroups(sim, self, eligible_module.criteria, seekers_dict)
            sa.TestAllocateGroups(sim, self, test_uid, eligible_module.criteria, eligible_dict)

        if len(test_results) == 0:  # Nothing need to be updated if no tests were conducted
            return
        
        # Record the positive results with a result latency term (swab delay is not included in this term, as we are not delaying test administration)
        if not(self.static_latency is None): 
            test_latencies = np.full(len(test_results), self.static_latency)
            unique_latencies = [self.static_latency]
        else: 
            test_latencies = np.trunc(stats.gamma.rvs(a=self.latency_shape, scale=self.latency_scale, size=len(test_results)))
            unique_latencies = np.unique(test_latencies).astype(np.intc)

        for l in unique_latencies:  # This should be a short loop
            if cur + l <= sim.pars['n_days']: 

                # Track case counts
                self.pos_history[cur + l] += sum(test_results[test_latencies == l])
                self.neg_history[cur + l] += sum(~test_results[test_latencies == l])  # Evaluates ~test_results before boolean indexing

                # Track dates
                self.date_positive[test_uid[np.logical_and(test_results, test_latencies == l)]] = cur + l  # Pick people who returned positive result and was assigned latency l 
                self.date_negative[test_uid[np.logical_and(~test_results, test_latencies == l)]] = cur + l  # Pick people who returned negative result and was assigned latency l 

            else:  # Pool tests that are recorded after the hypothetical outbreak ends
                self.pos_overflow += sum(test_results[test_latencies == l])
                self.neg_overflow += sum(~test_results[test_latencies == l])

        return 



class ActiveRAT(Surveillance): 
    '''
        Args: 
            sensitivity          (float): PCR test sensitivity
            specificity          (float): PCR test specificity 
            frequency            (int)  : Number of days between individual testing, or between policy rotations
            population           (list) : The populations that undergo active surveillance. Options are "all", "household", "school", "workplace", "community"
            policy               (str)  : Strategy deployed for active surveillance. Options are "cross_sectional" or "longitudinal"
            participation        (list) : Percentage of each population that undergoes active surveillance, with some frequency 
            age_report           (str)  : Which age distribution to use (report1, report2, report3, etc.)
            enforced             (bool) : Whether testing should be adjusted for compliance
            compliance           (float): Percentage of testing individuals that successfully conduct the test and report it. Used in the context of self-reporting.
            f_compliance         (float): Compliance at the end of the outbreak. Used in the context of longitudinal sampling, where individuals may drop out over time.
            latency              (dict)  : Number of days after taking the test before test results are considered reported, with associated probabilities. 
        '''

    def __init__(self,
                initial_unit,
                sensitivity=0.7, 
                specificity=0.99,
                frequency=3, 
                population=['all'],
                policy='cross_sectional',
                participation={
                    'all': 0.5,
                },
                age_report=None,
                enforced=True,
                compliance=1.0,
                f_compliance=1.0,
                latency_shape=2,
                latency_scale=1,
                static_latency=None):

        super().__init__(initial_unit)
        self.sensitivity = sensitivity 
        self.specificity = specificity 
        self.frequency = frequency 
        self.population = population
        self.policy = policy
        self.participation = participation 
        self.age_report = age_report
        self.enforced = enforced
        self.compliance = compliance 
        self.f_compliance = f_compliance
        self.latency_shape = latency_shape
        self.latency_scale = latency_scale
        self.static_latency = static_latency

        fname = '/Users/Ritchie/Desktop/cv-test/covasim/covasim/data/age_distribution_reports.JSON'
        with open(fname, 'r') as f: 
            self.all_age_reports = json.load(f)


    def initialize(self, sim):
        '''
        Args: 
            sim         (sim): Simulation object. Do not initialize simulation with start/end day, use n_day to avoid potential error)
        '''
        super().initialize(sim)  # Unlock people object and set initialized as true 

        # Initialize count trackers
        self.surv_history = np.full(sim.pars['n_days'] + 1, self.initial_unit)  # Surveillance report history, the simulation goes on for one extra day, hence +1
        self.pos_history = np.full(sim.pars['n_days'] + 1, self.initial_unit)
        self.neg_history = np.full(sim.pars['n_days'] + 1, self.initial_unit)

        # Initialize date trackers
        self.date_positive = np.full(sim.pars['pop_size'], -1)
        self.date_negative = np.full(sim.pars['pop_size'], -1)

        # Initialize overflow trackers 
        self.pos_overflow = 0
        self.neg_overflow = 0
        return


    def finalize(self, sim):
        '''
        Do any finalizing work here (e.g. resetting to save memory)
        '''
        super().finalize(sim)
        return


    def apply(self, sim, analysis=True): 
        cur = sim.t

        # Allocate tests according to designated policy 
        tester_module = sa.RAT_allocation(self, sim)
        test_uid = tester_module.cur_testers  # np array 

        # Filter individuals through compliance
        if not(self.enforced): 
            if self.policy == 'cross_sectional':  # Static compliance
                keep_inds = np.random.choice(len(test_uid), int(len(test_uid) * self.compliance), replace=False)
                test_uid = test_uid[keep_inds]
            elif self.policy == 'longitudinal':  # Linearly decreasing compliance
                keep_proportion = self.compliance - ((self.compliance - self.f_compliance) / sim.pars['n_days']) * cur  # Linearly decreasing compliance down to f_compliance
                keep_inds = np.random.choice(len(test_uid), int(len(test_uid) * keep_proportion), replace=False)
                test_uid = test_uid[keep_inds]

        # Test selected individuals
        test_results = st.basic_RAT_test(self, sim, test_uid, self.sensitivity, self.specificity)

        # Generate test latencies
        if not(self.static_latency is None): 
            test_latencies = np.full(len(test_results), self.static_latency)
            unique_latencies = [self.static_latency]
        else: 
            test_latencies = np.trunc(stats.gamma.rvs(a=self.latency_shape, scale=self.latency_scale, size=len(test_results)))
            unique_latencies = np.unique(test_latencies).astype(np.intc)

        for l in unique_latencies: 
            if cur + l <= sim.pars['n_days']: 
                # Track case counts
                self.pos_history[cur + l] += sum(test_results[test_latencies == l])
                self.neg_history[cur + l] += sum(~test_results[test_latencies == l])  # First evaluates ~test_results then ~test_results[latency == l]

                # Track dates
                self.date_positive[test_uid[np.logical_and(test_results, test_latencies == l)]] = cur + l  # Pick people who returned positive result and was assigned latency l 
                self.date_negative[test_uid[np.logical_and(~test_results, test_latencies == l)]] = cur + l  # Pick people who returned negative result and was assigned latency l 

            else:  # Pool tests that are recorded after the hypothetical outbreak ends
                self.pos_overflow += sum(test_results[test_latencies == l])
                self.neg_overflow += sum(~test_results[test_latencies == l])
        
        return


