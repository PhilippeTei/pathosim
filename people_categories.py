from . import utils as cvu

class PeopleCategories:
    """
    OBJECTIVE: Categorise people, for example if they have a sw alert. 
    SKELE: Basic categories: 
        - Tested positive, i.e. union of RAT and PCR positive testers. 
        - Have Covid Specific Symptoms
        - Have SW alert
    TASKS: 
        - Refactor and populate with additional categories. 
        - complete testers() function
    """
    def __init__(self, sim, criteria=None):
        if criteria is None:
            self.criteria = ['pos_test', 
                             'cs_sx', 
                             'sw_alert']
        else:
            # May want to only specify a few relevant criteria. 
            self.criteria = criteria

        self.cat_dict = self.check_criteria(sim)

    def check_criteria(self, sim): 
        '''
        Given the criterion names, call respective criterion functions to identify who meets each criterion. 

        Args: 
            sim     (simulation)    : Initialized simulation oboject
            surv    (surveillance)  : The PCR surveillance system being used 
        '''
        meets_criteria = dict()  # Will be a dictionary of sets

        for name in self.criteria:  # Not sure if there's a cleaner way to do this.
            if name == 'pos_test': 
                meets_criteria['pos_test'] = self.testers(sim, mode='pos')
            if name == 'cs_sx': 
                meets_criteria['cs_sx'] = self.sx(sim, mode='specific')
            if name == 'sw_alert': 
                meets_criteria['sw_alert'] = self.sw_alerts(sim, mode='pos')
        
        return meets_criteria

    def testers(self, sim, mode='pos'):
        '''
        Get people who have tested, + we can apply different filters. 

        Args: 
            sim          (sim) : Simulation object 
            mode         (str) : Choose return option 

        Returns: 
            sw_alerted  (set) : Set of people who have tested, with the following filters: 
                                    If mode == 'pos', returns set of all people who received a positive result on that day. 
                                    If mode == 'neg', returns set of all people who received a negative result on that day.
                                    If mode == 'RAT_pos' returns set of all people who received a RAT positive.
                                    If mode == 'RAT_neg' returns set of all people who received a RAT negative.
                                    If mode == 'PCR_pos' returns set of all people who received a PCR positive.
                                    If mode == 'PCR_neg' returns set of all people who received a PCR negative.
        '''
        return_options = ['pos'] # all others are to be implemented.

        if not(mode in return_options): 
            raise RuntimeError('Return choice is invalid, the options are: "all", "specific, "nonspecific"')

        if mode == 'pos':
            return sim.people.date_diagnosed == sim.t

    def sw_alerts(self, sim, mode='pos'):
        '''
        Get people who have positive smartwatch alerts. 

        Args: 
            sim          (sim) : Simulation object 
            mode         (str) : Choose return option 

        Returns: 
            sw_alerted  (set) : Set of people who have smartwatches and ...
                                    If mode == 'pos', returns set of people who have positive alerts on the day. 
                                    If mode == 'neg', returns set of watch users without alerts on the day. 
        '''
        return_options = ['pos'] # neg is to be implemented. 

        if not(mode in return_options): 
            raise RuntimeError('Return choice is invalid, the options are: "all", "specific, "nonspecific"')

        alerted = sim.people.true('sw_alarmed')
        if mode == 'pos':
            return alerted

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