from . import people_categories as cvp


class TestScheduler:
    """
    Initialized once, in sim.py. 
    Call schedule_tests() once per timestep, ex in sim.py:step. 
    A wrapper function for pcr and rat models that splits people between them. 
    """

    def __init__(self, **kwargs):
        """
        __init__ requires pre-initialized pcr and/or rat models to be passed in.
        """

        self.make_pars() # Populate fields with default values. 
        self.update_pars(**kwargs)

        # Determine which regime we're in.
        if self.pcr_model == None:
            self.p_rat_given_t = 1
        if self.rat_model == None:
            self.p_rat_given_t = 0
        if self.pcr_model == None and self.rat_model == None:
            raise ValueError('No PCR or RAT model specified.')

    def make_pars(self):
        self.pcr_model = None
        self.rat_model = None
        self.p_rat_given_t = 0.8 # People prefer RAT with 80% chance. 
        self.t_seek_probs = { # TODO: Update. 
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
        } # Test seeking probabilities.
        self.pcr_eligiility_criteria = list(self.t_seek_probs.keys()) # Assume everyone eligible by default 

    def schedule_tests(self, sim):
        """
        Call once per timestep. 
        """
        people_cat_dict = cvp.PeopleCategories(sim.people).cat_dict
        test_seekers = TestSeekers(people_cat_dict, self.t_seek_probs)
        
        # Choose who goes to RAT, and who goes to PCR.
        rat_inds, pcr_inds = self.split_test_seekers(test_seekers)

        # Pipe some to RAT.
        if self.pcr_model is not None:
            self.pcr_model.schedule_tests(sim, pcr_inds)
        if self.rat_model is not None:
            self.rat_model.schedule_tests(sim, rat_inds)

    def split_test_seekers(self, test_seekers):
        """
        Args: 
            test_seekers: a dictionary of who seeks tests.
        Returns:
            rat_inds: indices of people who go to RAT.
            pcr_inds: indices of people who go to PCR.
        Body Notes:
            - p_rat_given_t must be 1 for people ineligibile for PCR.
            - For those who are eligible for both, toss a coin (p_rat_given_t i.e. given they're a test seeker)
            - Can go through all sets and choose randomly for new people encountered.
            - Would be faster for bitmask representation of set membership. O(<=n) because vectorised. 
        """

        pcr_crit = self.pcr_eligiility_criteria
        p_rat_given_t = self.p_rat_given_t
        raise NotImplementedError('split_test_seekers not implemented.')

    def update_pars(self, **kwargs):
        """
        Args: 
            **kwargs: keyword arguments to update parameters.
        Body:
            Update fields with parameters. New fields cannot be initialized in this way. 
        """
        for k, v in kwargs.items():
            # Check if k is already a field. 
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise ValueError('Invalid parameter: %s' % k)

class Test:
    def __init__():
        pass

class RAT(Test): # TODO: Rename Surveillance class to Test.
    """
    Basically the same as surveillance, but without the eligibility and test seeking steps. 
    """
    def __init__():
        pass
    def initialize():
        pass
    def finalize():
        pass
    def schedule_tests(sim, tester_inds):
        # eligibility() and test_seekers() are already handled by TestScheduler
        # Determine who is allocated a test
        # Administer the tests & implement delay
        # Update state variables sim.people.date_diagnosed and sim.people.date_pos_test. (See interventions:test_prob for an example)
        pass

class PCR(Test):
    def __init__():
        pass
    def initialize():
        pass
    def finalize():
        pass
    def schedule_tests(sim, tester_inds):
        # eligibility() and test_seekers() are already handled by TestScheduler
        # Determine who is allocated a test
        # Administer the tests & implement delay
        # Update state variables sim.people.date_diagnosed and sim.people.date_pos_test. (See interventions:test_prob for an example)
        pass

class TestSeekers:
    """
    Objective: trim down the dictionary of categories. 
    """
    def __init__(self, people_cats, seek_probs=None): 
            """
            People cats: dictionary of categories, with values being the indicies of people. 
            """
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

            self.all_seekers = self.get_all_seekers(people_cats)
    
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
