import numpy as np
import pandas as pd 
from scipy import stats
from pathosim import base_covid_testing as bct
from pathosim.test_pipeline import allocate as ac
from pathosim.test_pipeline import test_mechanism as tm

master_criteria = ['cont_conf', 
                   'cont_conf_quar',
                   'cont_vuln', 
                   'cont_sx', 
                   'cont_cs_sx', 
                   'cont_ncs_sx', 
                   'sx', 
                   'cs_sx', 
                   'ncs_sx', 
                   'sx_quar',
                   'sev_crit',
                   'pos_RAT', 
                   'pos_RAT_asx', 
                   'neg_RAT_sx', 
                   'work',
                   'swa_RAT',
                   'swa_PCR']

master_seekprobs = {
                'cont_conf': 0.8, 
                'cont_conf_quar': 0.8,
                'cont_vuln': 0.2,
                'cont_sx': None, 
                'cont_cs_sx': 0.3,
                'cont_ncs_sx': 0.1,
                'sx': None, 
                'cs_sx': 0.8,
                'ncs_sx': 0.4,
                'sx_quar': 0.8,
                'sev_crit': 0.8,
                'pos_RAT': 0.5,
                'pos_RAT_asx': 0.8,
                'neg_RAT_sx': 0.8,
                'neg_RAT': 0.4,
                'other': 0.4,
                'work': 1.0,
                'school': 1.0,
                'swa_RAT': 0.8,
                'swa_PCR': 0.4,
            }

default_PCR_disc_criteria = {'diagnostic': ['cont_conf', 'cs_sx', 'pos_RAT_asx', 'neg_RAT_sx', 'other'], 
                             'screening': ['work']}
default_RAT_disc_criteria = {'diagnostic': ['cont_conf', 'cont_vuln', 'cont_cs_sx', 'cont_ncs_sx', 'cs_sx', 'ncs_sx', 'neg_RAT', 'other'], 
                             'screening': ['work']}

default_PCR_disc_seekprobs = {'diagnostic': {c: master_seekprobs[c] for c in default_PCR_disc_criteria['diagnostic']},
                              'screening': {c: master_seekprobs[c] for c in default_PCR_disc_criteria['screening']}}
default_RAT_disc_seekprobs = {'diagnostic': {c: master_seekprobs[c] for c in default_RAT_disc_criteria['diagnostic']},
                              'screening': {c: master_seekprobs[c] for c in default_RAT_disc_criteria['screening']}}


class TestObj: 
    def __init__(self): 
        self.system = self.__class__.__name__
        self.initialized = False
        self.finalized = False
        return
    
    def __call__(self, *args, **kwargs):
        # Makes TestObj(sim) equivalent to TestObj.apply(sim)
        if not self.initialized:  # pragma: no cover
            errormsg = f'TestObj (label={self.system}, {type(self)}) has not been initialized'
            raise RuntimeError(errormsg)
        return self.apply(*args, **kwargs)

    def update_pars(self, pars_dict): 
        for param_name in pars_dict.keys():
            if hasattr(self, param_name): 
                setattr(self, pars_dict[param_name])
            else: 
                raise RuntimeError(f'Tried to assign {param_name} to testobj, but testobj does not have that attribute')

    def finalize(self):
        '''
        Finalize the TestObj object
        '''
        if self.finalized: # pragma: no cover
            raise RuntimeError('TestObj already finalized')  # Raise an error because finalizing multiple times has a high probability of producing incorrect results e.g. applying rescale factors twice
        self.finalized = True
        return


class PCR_disc(TestObj): 

    def __init__(self, 
                di_criteria=default_PCR_disc_criteria['diagnostic'], 
                di_seekprobs=None,
                sc_criteria=default_PCR_disc_criteria['screening'], 
                sc_seekprobs=None,
                crit_kwargs=None,
                sensitivity=0.85,
                specificity=1.00,
                LOD=3,  # Log base 10
                capacity=0.015,
                pos_retest_delay=90,  
                neg_retest_delay=14,  
                mean_latency=3,
                static_latency=None,
                RAT_ind=None,
                timeline=None,
                test_mode='basic',
                name='PCR'):
        '''
        PCR-based 'di'agnostic testing and 'sc'reening.

        Args: 
            system             (str)    : Test system name
            di_criteria        (list)   : List of eligibility criteria used for diagnostic testing
            di_seekprobs       (dict)   : Dictionary of test seeking probabilities (value) associated with each criteria (key)
            sc_criteria        (list)   : List of eligibiltiy criteria used for screening 
            sc_seekprobs       (dict)   : Dictionary of test seeking probabilities (value) associated with each criteria (key) 
            sensitivity        (float)  : Test sensitivity
            specificity        (float)  : Test specificity 
            capacity           (n/s)    : Test capacity as a proportion of total population (float or array)
            pos_retest_delay   (int)    : Number of days before an individual can retest, after receiving a positive test
            neg_retest_delay   (int)    : Number of days before an individual can retest, after receiving a negative test
            latency_shape      (int)    : Shape of gamma distribution for test reporting latency 
            latency_scale      (int)    : Scale of gamma distribution for test reporting latency 
            static_latency     (int)    : Specificy a static test reporting latency (optional)
            RAT_ind            (int)    : Index of RAT testing system, which also inherits from TestObj. Usually RAT_disc, can also be RAT_surv. 
            hybrid             (bool)   : Whether the testing system is a PCR-RAT hybrid 
            timeline           (dict)   : Parameter updates at different times 
        '''
        super().__init__()
        self.di_criteria = di_criteria
        self.di_seekprobs = di_seekprobs
        self.sc_criteria = sc_criteria
        self.sc_seekprobs = sc_seekprobs
        self.crit_kwargs = crit_kwargs
        self.sensitivity = sensitivity
        self.specificity = specificity
        self.LOD = LOD
        self.capacity = capacity 
        self.pos_retest_delay = pos_retest_delay 
        self.neg_retest_delay = neg_retest_delay 
        self.mean_latency = mean_latency
        self.latency_shape = 2  # fixed
        self.latency_scale = self.mean_latency / self.latency_shape
        self.static_latency = static_latency
        self.RAT_ind = RAT_ind  # For ease of use, we could remove this and hide it from user, force RAT_ind=0 always. 
        self.RAT = None
        self.hybrid = False
        self.timeline = timeline
        self.test_mode = test_mode
        self.name = name

        if self.di_seekprobs is None: 
            self.di_seekprobs = {c: master_seekprobs[c] for c in self.di_criteria} 
        if self.sc_seekprobs is None: 
            self.sc_seekprobs = {c: master_seekprobs[c] for c in self.sc_criteria} 

    # TODO: Decide which trackers are unnecessary for testing, but necessary for surveillance
    def initialize(self, sim): 
        self.diagnostic = bct.Diagnostic(self.di_criteria, self.di_seekprobs)
        self.screening = bct.Screening(self.sc_criteria, self.sc_seekprobs)

        # Initialize base test objects
        self.diagnostic.initialize()
        self.screening.initialize(sim.pars['pop_type'], sim.people)

        # Initialize individual trackers
        self.pos_agents = np.zeros((sim.pars['n_days'] + 1, sim.pars['pop_size']))  
        # Want something like: self.true_pos_agents, in order to calculate the positive predictve value

        # Initialize count trackers
        self.pos_history = np.zeros(sim.pars['n_days'] + 1)  # The simulation always goes on for one extra day, hence +1
        self.neg_history = np.zeros(sim.pars['n_days'] + 1)

        # Initialize date trackers
        self.date_positive = np.full(sim.pars['pop_size'], -1)  # For integration 
        self.date_negative = np.full(sim.pars['pop_size'], -1)
        self.date_pos_test = np.full(sim.pars['pop_size'], -1)  # For integration

        # Initialize overflow trackers
        self.pos_overflow = 0
        self.neg_overflow = 0

        # Initialize other trackers
        # May not need these first two trackers for PCR testing
        self.neg_was_sx = np.zeros(sim.pars['pop_size'])  # [Legend] 0: Never tested before, 1: Tested and was sx, 2: Tested and was asx
        self.pos_was_asx = np.zeros(sim.pars['pop_size'])
        self.retest_tracker = np.zeros(sim.pars['pop_size'])  # Track whether an individual is allowed to retest 

        # Initialize two-way interface with another testing system - assumed to be RAT currently, can it be PCR? 
        if not(self.RAT_ind is None): 
            self.RAT = sim.pars['testing'][self.RAT_ind] 
            self.RAT.PCR = self

            # Activate
            self.RAT.hybrid = True
            self.hybrid = True

        self.initialized = True
        

    def apply(self, sim): 
        '''
        Every day: 
            1) Get people eligible for diagnostic and screening
            2) Allocate tests to test-seekers
            3) Conduct testing
        '''
        if not(self.timeline is None): 
            update_parameters(self, sim.t)

            # Refresh diagnostic and screening objects with any updates to criteria/probabilities
            self.diagnostic.criteria = self.di_criteria 
            self.diagnostic.seekprobs = self.di_seekprobs
            self.screening.criteria = self.sc_criteria
            self.screening.seekprobs = self.sc_seekprobs

        # Update retest restrictions 
        self.retest_tracker -= 1

        di_eg, di_sk = self.diagnostic.apply(self, sim)
        sc_eg, sc_sk = self.screening.apply(self, sim)
        
        # Merge dictionaries
        all_eg = {**di_eg, **sc_eg}
        all_sk = {**di_sk, **sc_sk}

        # Allocate tests, conduct tests, record test results
        test_uids = ac.allocate_tests(sim, self, all_sk, master_seekprobs, mode='proportional')
        if self.test_mode == 'basic':
            test_results = tm.basic_test(self, sim, test_uids)
        elif self.test_mode == 'LOD':
            test_results = tm.LOD_test(self, sim, test_uids)
        else: 
            raise RuntimeError('Invalid test mode supplied. Options are "basic" and "LOD"')
        
        #TODO: Why are there still so many people receiving several diagnoses on the same day?
        # Update retest restrictions
        self.retest_tracker[test_uids[test_results]] = self.pos_retest_delay
        self.retest_tracker[test_uids[~test_results]] = self.neg_retest_delay

        # Record information - TODO: Optimize record_pipeline
        record_pipeline(self, sim, all_eg, all_sk, test_uids)
        record_test_results(self, sim, test_uids, test_results)  # Latencies get recorded here
        record_consumed(self, sim, len(test_uids))

        sim.results['new_tests'][sim.t] += len(test_uids)

        return

    def finalize(self, sim):
        '''
        Do any finalizing work here (e.g. resetting to save memory)
        '''
        super().finalize()
        return


class RAT_disc(TestObj): 

    def __init__(self, 
                di_criteria=default_RAT_disc_criteria['diagnostic'], 
                di_seekprobs=None,
                sc_criteria=default_RAT_disc_criteria['screening'], 
                sc_seekprobs=None,
                crit_kwargs=None,
                sensitivity=0.70,
                specificity=0.99,
                LOD=5,  # Log scale
                capacity=0.05,
                pos_retest_delay=14,  
                neg_retest_delay=0,  
                mean_latency = 0,
                static_latency=None,
                timeline=None,
                test_mode='basic',
                name='RAT'):
        '''
        RAT-based 'di'agnostic testing and 'sc'reening.
        '''
        super().__init__()
        self.di_criteria = di_criteria
        self.di_seekprobs = di_seekprobs
        self.sc_criteria = sc_criteria
        self.sc_seekprobs = sc_seekprobs
        self.crit_kwargs = crit_kwargs
        self.sensitivity = sensitivity
        self.specificity = specificity
        self.LOD = LOD
        self.capacity = capacity
        self.pos_retest_delay = pos_retest_delay
        self.neg_retest_delay = neg_retest_delay
        self.mean_latency = mean_latency
        self.latency_shape = 2
        self.latency_scale = self.mean_latency / self.latency_shape
        self.static_latency = static_latency
        self.PCR = None
        self.RAT = None
        self.hybrid = False
        self.timeline = timeline
        self.test_mode = test_mode
        self.name = name

        if self.di_seekprobs is None: 
            self.di_seekprobs = {c: master_seekprobs[c] for c in self.di_criteria} 
        if self.sc_seekprobs is None: 
            self.sc_seekprobs = {c: master_seekprobs[c] for c in self.sc_criteria} 

    def initialize(self, sim): 
        self.diagnostic = bct.Diagnostic(self.di_criteria, self.di_seekprobs)
        self.screening = bct.Screening(self.sc_criteria, self.sc_seekprobs)

        # Initialize base test objects
        self.diagnostic.initialize()
        self.screening.initialize(sim.pars['pop_type'], sim.people)

        # Initialize individual trackers
        self.pos_agents = np.zeros((sim.pars['n_days'] + 1, sim.pars['pop_size']))  

        # Initialize count trackers
        self.pos_history = np.zeros(sim.pars['n_days'] + 1)  # The simulation always goes on for one extra day, hence +1
        self.neg_history = np.zeros(sim.pars['n_days'] + 1)

        # Initialize date trackers
        self.date_positive = np.full(sim.pars['pop_size'], -1)  # For integration 
        self.date_negative = np.full(sim.pars['pop_size'], -1)
        self.date_pos_test = np.full(sim.pars['pop_size'], -1)  # For integration

        # Initialize overflow trackers
        self.pos_overflow = 0
        self.neg_overflow = 0

        # Initialize other trackers
        self.neg_was_sx = np.zeros(sim.pars['pop_size'])  # [Legend] 0: Never tested before, 1: Tested and was sx, 2: Tested and was asx
        self.pos_was_asx = np.zeros(sim.pars['pop_size'])
        self.retest_tracker = np.zeros(sim.pars['pop_size'])  # Track whether an individual is allowed to retest
        self.cons_days_neg_rat = np.zeros(sim.pars['pop_size'])  # Track the number of consecutive days individual has received negative RAT

        self.initialized = True

    def apply(self, sim): 
        '''
        Every day: 
            1) Get people eligible for diagnostic and screening
            2) Allocate tests to test-seekers
            3) Conduct testing
        '''
        if not(self.timeline is None): 
            update_parameters(self, sim.t)

            # Refresh diagnostic and screening objects with any updates to criteria/probabilities
            self.diagnostic.criteria = self.di_criteria
            self.diagnostic.seekprobs = self.di_seekprobs
            self.screening.criteria = self.sc_criteria
            self.screening.seekprobs = self.sc_seekprobs

        # Update retest restrictions 
        self.retest_tracker -= 1

        di_eg, di_sk = self.diagnostic.apply(self, sim)
        sc_eg, sc_sk = self.screening.apply(self, sim)
        
        # Merge dictionaries
        all_eg = {**di_eg, **sc_eg}
        all_sk = {**di_sk, **sc_sk}

        # Allocate tests, conduct tests, record test results
        test_uids = ac.allocate_tests(sim, self, all_sk, master_seekprobs, mode='proportional')
        if self.test_mode == 'basic':
            test_results = tm.basic_test(self, sim, test_uids)
        elif self.test_mode == 'LOD':
            test_results = tm.LOD_test(self, sim, test_uids)
        else: 
            raise RuntimeError('Invalid test mode supplied. Options are "basic" and "LOD"')

        #TODO: Why are there still so many people receiving several diagnoses on the same day?
        # Update retest restrictions
        self.retest_tracker[test_uids[test_results]] = self.pos_retest_delay
        self.retest_tracker[test_uids[~test_results]] = self.neg_retest_delay
        # print(self.retest_tracker)  # Confirm that boolean indexing is functioning correctly 

        # Record information
        record_pipeline(self, sim, all_eg, all_sk, test_uids)
        record_test_results(self, sim, test_uids, test_results)
        record_consumed(self, sim, len(test_uids))

        # Track consecutive days negative (assuming that an individual instantly knows test result) 
        tested_negative = np.zeros(sim.pars['pop_size'], dtype=bool) 
        tested_negative[test_uids[~test_results]] = True
        self.cons_days_neg_rat[tested_negative] += 1
        self.cons_days_neg_rat[~tested_negative] = 0  # reset everyone else

    def finalize(self, sim):
        '''
        Do any finalizing work here (e.g. resetting to save memory)
        '''
        super().finalize()
        return


class RAT_surv(TestObj): 

    def __init__(self,
                sensitivity=0.7, 
                specificity=0.99,
                LOD=5,
                frequency=7, 
                population=['all'],
                policy='cross_sectional',
                participation={
                    'all': 0.5,
                },
                age_report=None,
                enforced=True,
                i_compliance=1.0,
                f_compliance=1.0,
                mean_latency=0,
                static_latency=None,
                name='RAT'):
        '''
        RAT-based surveillance testing.
        '''
        super().__init__()
        self.sensitivity = sensitivity 
        self.specificity = specificity 
        self.LOD = LOD
        self.frequency = frequency 
        self.population = population
        self.policy = policy 
        self.participation = participation 
        self.age_report = age_report
        self.enforced = enforced 
        self.i_compliance = i_compliance 
        self.f_compliance = f_compliance 
        self.mean_latency = mean_latency
        self.latency_shape = 2,
        self.latency_scale = self.mean_latency / self.latency_shape
        self.static_latency = static_latency
        self.name = name
        
    def initialize(self, sim): 
        self.surveillance = bct.Surveillance()

        # Initialize base test objects
        self.surveillance.initialize()  # Currently this attribute is useless

        # Initialize individual trackers
        self.pos_agents = np.zeros((sim.pars['n_days'] + 1, sim.pars['pop_size']))  # Tracks the positive status of every individual over time

        # Initialize count trackers
        self.pos_history = np.zeros(sim.pars['n_days'] + 1)  # The simulation always goes on for one extra day, hence +1
        self.neg_history = np.zeros(sim.pars['n_days'] + 1)

        # Initialize date trackers
        self.date_positive = np.full(sim.pars['pop_size'], -1)  # For integration 
        self.date_negative = np.full(sim.pars['pop_size'], -1)
        self.date_pos_test = np.full(sim.pars['pop_size'], -1)  # For integration

        # Initialize overflow trackers
        self.pos_overflow = 0
        self.neg_overflow = 0

        # Initialize other trackers
        self.neg_was_sx = np.zeros(sim.pars['pop_size'])  # [Legend] 0: Never tested before, 1: Tested and was sx, 2: Tested and was asx
        self.pos_was_asx = np.zeros(sim.pars['pop_size'])

        self.initialized = True

    def apply(self, sim): 
         # Allocate tests, conduct tests, record test results
        test_uids = ac.allocate_tests(sim, self, mode=self.policy)
        test_results = tm.LOD_test(self, sim, test_uids)

        record_test_results(self, sim, test_uids, test_results)
        record_consumed(self, sim, len(test_uids))
        return

    def finalize(self, sim):
        '''
        Do any finalizing work here (e.g. resetting to save memory)
        '''
        super().finalize()
        return


# Utils
def record_test_results(testobj, sim, test_uids, test_results): 
    cur = sim.t
    
    # Record the positive results with a result latency term (swab delay is not included in this term, as we are not delaying test administration)
    if not(testobj.static_latency is None): 
        test_latencies = np.full(len(test_results), testobj.static_latency)
        unique_latencies = [testobj.static_latency]
    else: 
        test_latencies = np.trunc(stats.gamma.rvs(a=testobj.latency_shape, scale=testobj.latency_scale, size=len(test_results)))
        unique_latencies = np.unique(test_latencies).astype(np.intc)

    # Record date of tests that will eventually return positive (i.e., the date before latency is introduced)
    testobj.date_pos_test[test_uids[test_results]] = cur

    # Record dates that involve a latency - these dates are when test results become known to the individual
    for l in unique_latencies:  # This should be a short loop
        if cur + l <= sim.pars['n_days']:

            # Track case counts
            testobj.pos_history[cur + l] += sum(test_results[test_latencies == l])
            testobj.neg_history[cur + l] += sum(~test_results[test_latencies == l])  # Evaluates ~test_results before boolean indexing
            
            # Track dates
            testobj.date_positive[test_uids[np.logical_and(test_results, test_latencies == l)]] = cur + l  # Pick people who returned positive result and was assigned latency l 
            testobj.date_negative[test_uids[np.logical_and(~test_results, test_latencies == l)]] = cur + l  # Pick people who returned negative result and was assigned latency l ]

            # Track agent test results. Represents the date that a patient is officially diagnosed in a binary matrix.
            testobj.pos_agents[cur + l, test_uids[np.logical_and(test_results, test_latencies == l)]] += 1

        else:  # Pool tests that are recorded after the outbreak ends
            testobj.pos_overflow += sum(test_results[test_latencies == l])
            testobj.neg_overflow += sum(~test_results[test_latencies == l])
    return 

def record_pipeline(testobj, sim, all_eligible, all_seekers, test_uids):
    criteria = all_eligible.keys()
    recordEligibleGroups(testobj, sim, criteria, all_eligible)
    recordSeekGroups(testobj, sim, criteria, all_seekers)
    recordAllocateGroups(testobj, sim, criteria, all_eligible, test_uids)
    return 

def recordEligibleGroups(testobj, sim, criteria, all_eligible): 
    '''
    Get the number of people who meet each criteria each day. 
    '''
    if sim.t == 0: 
        setattr(testobj, 'crit_groups', dict())
        for criterion in criteria: 
            testobj.crit_groups[criterion] = []

    for criterion in testobj.crit_groups.keys():
        testobj.crit_groups[criterion].append(len(all_eligible[criterion])) 


def recordSeekGroups(testobj, sim, criteria, all_seekers): 
    '''
    For each testing criterion, get the number of test seekers belonging to the group. 
    '''
    if sim.t == 0: 
        setattr(testobj, 'seek_groups', dict())
        for criterion in criteria: 
            testobj.seek_groups[criterion] = []

    # Get all test seekers at the current state of the simulation
    for criterion in testobj.seek_groups.keys(): 
        testobj.seek_groups[criterion].append(len(all_seekers[criterion]))


def recordAllocateGroups(testobj, sim, criteria, all_eligible, testers): 
    '''
    Among those that received a PCR test, identify which criteria each person satisfied - this is their "reason" for testing. 
    '''
    if sim.t == 0: 
        setattr(testobj, 'alloc_groups', dict())
        for criterion in criteria: 
            testobj.alloc_groups[criterion] = []
    
    testers = set(testers)

    for criterion in testobj.alloc_groups.keys():
        testobj.alloc_groups[criterion].append(len(testers.intersection(all_eligible[criterion])))

    return


def record_consumed(testobj, sim, tests): 
    if sim.t == 0: 
        setattr(testobj, 'tests_consumed', [])

    testobj.tests_consumed.append(tests)


def update_parameters(testobj, t):
    '''
    Check if parameter updates have been scheduled for the current day, and update accordingly. 
    '''
    # Throw error if first element of keys is a string. (only numbers work)
    # First check if the dict is nonempty
    if len(testobj.timeline.keys()) > 0:
        if isinstance(list(testobj.timeline.keys())[0], str):
            raise TypeError('Timeline keys must be integers.')

    if t in testobj.timeline.keys(): 
        for param, param_val in testobj.timeline[t].items():
            if not(param in testobj.__dict__.keys()):
                raise RuntimeError(f'The {param} parameter of {testobj.__class__} does not exist, update attempt failed.')
            else: 
                setattr(testobj, param, param_val)
        print(f'Parameters for {testobj.__class__} have been updated!')


# Old function
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
