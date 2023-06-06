'''
Defines the People class and functions associated with making people and handling
the transitions between states (e.g., from susceptible to infected).
'''

#%% Imports
import numpy as np
import scipy.stats as stats
import sciris as sc
import scipy.stats as stats
from collections import defaultdict
from . import version as cvv
from . import utils as cvu
from . import defaults as cvd
from . import base as cvb
from . import plotting as cvplt
from . import immunity as cvi
from . import watches as cvw
from . import pathogens as pat

__all__ = ['People']

class People(cvb.BasePeople):
    '''
    A class to perform all the operations on the people -- usually not invoked directly.

    This class is usually created automatically by the sim. The only required input
    argument is the population size, but typically the full parameters dictionary
    will get passed instead since it will be needed before the People object is
    initialized. However, ages, contacts, etc. will need to be created separately --
    see ``cv.make_people()`` instead.

    Note that this class handles the mechanics of updating the actual people, while
    ``cv.BasePeople`` takes care of housekeeping (saving, loading, exporting, etc.).
    Please see the BasePeople class for additional methods.

    Args:
        pars (dict): the sim parameters, e.g. sim.pars -- alternatively, if a number, interpreted as pop_size
        strict (bool): whether or not to only create keys that are already in self.meta.person; otherwise, let any key be set
        kwargs (dict): the actual data, e.g. from a popdict, being specified

    **Examples**::

        ppl1 = cv.People(2000)

        sim = cv.Sim()
        ppl2 = cv.People(sim.pars)
    '''

    def __init__(self, pars, strict=True, **kwargs):
        
        # Temp variables meant for debugging.
        # print("DEBUG: TRACK ALERT TIMINGS")
        # self.tp_alarm_rel_infection_date = []
        # print("DEBUG: SMARTWATCH ALERT PROBS")
        # days = np.arange(-21, 22, 1)
        # # dict of day:0
        # self.sum_alert_on_day = {day:0 for day in days}
        # self.sum_people_on_day = {day:0 for day in days}

        # self.all_peoples_infection_days = np.empty(20000)
        ##### END DEBUG

        # Handle pars and population size
        self.set_pars(pars)
        self.version = cvv.__version__ # Store version info

        # Other initialization
        self.t = 0 # Keep current simulation time
        self._lock = False # Prevent further modification of keys
        self.meta = cvd.PeopleMeta() # Store list of keys and dtypes
        self.contacts = None
        self.init_contacts() # Initialize the contacts
        self.infection_log = [] # Record of infections - keys for ['source','target','date','layer']
        self.stratifications = None # Gets updated in sim.py : initialize() 

        # Set person properties -- all floats except for UID
        for key in self.meta.person:
            if key == 'uid':
                self[key] = np.arange(self.pars['pop_size'], dtype=cvd.default_int)
            elif key in ['n_infections', 'n_breakthroughs']:
                self[key] = np.zeros(self.pars['pop_size'], dtype=cvd.default_int)
            elif key in ['cons_days_in_quar', 'cons_days_neg_rat']:
                self[key] = np.zeros(self.pars['pop_size'], dtype=cvd.default_int)
            elif key in ['has_watch']:
                self[key] = np.full(self.pars['pop_size'], False, dtype=bool)
            elif key in ['income']:
                self[key] = np.zeros(self.pars['pop_size'], dtype=cvd.default_int)
            elif key in ['viral_load']:
                self[key] = np.zeros(self.pars['pop_size'], dtype=cvd.default_float)
            else:
                self[key] = np.full(self.pars['pop_size'], np.nan, dtype=cvd.default_float)

        # Set health states -- only susceptible is true by default -- booleans except exposed by variant which should return the variant that ind is exposed to
        for key in self.meta.states:
            val = (key in ['susceptible', 'naive']) # Default value is True for susceptible and naive, False otherwise
            self[key] = np.full(self.pars['pop_size'], val, dtype=bool)

        # Set variant states, which store info about which variant a person is exposed to
        for key in self.meta.variant_states:
            self[key] = np.full(self.pars['pop_size'], np.nan, dtype=cvd.default_float)
        for key in self.meta.by_variant_states:
            self[key] = np.full((self.pars['n_variants'], self.pars['pop_size']), False, dtype=bool)

        # Set immunity and antibody states
        for key in self.meta.imm_states:  # Everyone starts out with no immunity
            self[key] = np.zeros((self.pars['n_variants'], self.pars['pop_size']), dtype=cvd.default_float)
        for key in self.meta.nab_states:  # Everyone starts out with no antibodies
            dtype = cvd.default_int if key == 't_nab_event' else cvd.default_float
            self[key] = np.zeros(self.pars['pop_size'], dtype=dtype)
        for key in self.meta.vacc_states:
            self[key] = np.zeros(self.pars['pop_size'], dtype=cvd.default_int)

        # Set dates and durations -- both floats
        for key in self.meta.dates + self.meta.durs:
            self[key] = np.full(self.pars['pop_size'], np.nan, dtype=cvd.default_float)

        # Set dates for viral load profile -- floats
        for key in self.meta.vl_points:
            self[key] = np.full(self.pars['pop_size'], np.nan, dtype=cvd.default_float)

        # Store the dtypes used in a flat dict
        self._dtypes = {key:self[key].dtype for key in self.keys()} # Assign all to float by default
        if strict:
            self.lock() # If strict is true, stop further keys from being set (does not affect attributes)

        # Store flows to be computed during simulation
        self.init_flows()

        # Although we have called init(), we still need to call initialize()
        self.initialized = False 

        # Handle contacts, if supplied (note: they usually are)
        if 'contacts' in kwargs:
            self.add_contacts(kwargs.pop('contacts'))

        if 'workplaces' in kwargs:
            self.workplaces = kwargs.pop('workplaces') 
        if 'n_workplaces' in kwargs:
            self.n_workplaces = kwargs.pop('n_workplaces')


        # Handle all other values, e.g. age
        for key,value in kwargs.items():
            if strict:
                self.set(key, value)
            else:
                self[key] = value

        self._pending_quarantine = defaultdict(list)  # Internal cache to record people that need to be quarantined on each timestep {t:(person, quarantine_end_day)}

        return


    def init_watches(self, pars):
        '''
        Initialize the watches for the people.
        '''

        self.watches = cvw.Watches(pars)
        
        # Initialize the people-specific P(non-inf alert)
        n = len(self.uid)
        if self.watches.pars['use_variable_fpr']:
            self.watch_fpr = self.watches.pars['fpr_sampler'].sample(n=n)
        else:
            self.watch_fpr = np.ones(n) * self.watches.pars['mean_fpr']
        
        # Update the number of people who actually download the app and use it.
        # We model this by setting all other people to "not have watches". 
        usage_rate = self.watches.pars['usage_rate']
        mask = np.random.binomial(1, usage_rate, n) # Get a binomial mask
        self.has_watch = np.logical_and(self.has_watch, mask) # Update self.has_watch

        # Bad coding practice to init here, but much simpler. TODO: Verify.
        self.LEN_ALERT_HIST = 3
        self.alert_histories = np.zeros((len(self.uid), self.LEN_ALERT_HIST))
        return


    def init_flows(self):
        ''' Initialize flows to be zero '''
        self.flows = {key:0 for key in cvd.new_result_flows}
        self.flows_variant = {}
        for key in cvd.new_result_flows_by_variant:
            self.flows_variant[key] = np.zeros(self.pars['n_variants'], dtype=cvd.default_float)
        if self.pars['enable_smartwatches']:
            for key in cvd.new_result_flows_smartwatches:
                self.flows[key] = 0
        if self.pars['enable_multiregion']:
            targets = ['new_diagnoses', 'new_severe']
            regions = self.pars['multiregion']['rnames']
            for target in targets:
                for region in regions:
                    self.flows[f'{region}_{target}'] = 0
        return

    def initialize(self, sim_pars=None):
        ''' Perform initializations '''
        self.validate(sim_pars=sim_pars) # First, check that essential-to-match parameters match
        self.set_pars(sim_pars) # Replace the saved parameters with this simulation's
        self.set_prognoses()
        self.initialized = True
        return


    def set_prognoses(self):
        '''
        Set the prognoses for each person based on age during initialization. Need
        to reset the seed because viral loads are drawn stochastically.
        '''

        pars = self.pars # Shorten
        if 'prognoses' not in pars or 'rand_seed' not in pars:
            errormsg = 'This people object does not have the required parameters ("prognoses" and "rand_seed"). Create a sim (or parameters), then do e.g. people.set_pars(sim.pars).'
            raise sc.KeyNotFoundError(errormsg)

        def find_cutoff(age_cutoffs, age):
            '''
            Find which age bin each person belongs to -- e.g. with standard
            age bins 0, 10, 20, etc., ages [5, 12, 4, 58] would be mapped to
            indices [0, 1, 0, 5]. Age bins are not guaranteed to be uniform
            width, which is why this can't be done as an array operation.
            '''
            return np.nonzero(age_cutoffs <= age)[0][-1]  # Index of the age bin to use

        cvu.set_seed(pars['rand_seed'])

        progs = pars['prognoses'] # Shorten the name
        inds = np.fromiter((find_cutoff(progs['age_cutoffs'], this_age) for this_age in self.age), dtype=cvd.default_int, count=len(self)) # Convert ages to indices
        self.symp_prob[:]   = progs['symp_probs'][inds] # Probability of developing symptoms
        self.severe_prob[:] = progs['severe_probs'][inds]*progs['comorbidities'][inds] # Severe disease probability is modified by comorbidities
        self.crit_prob[:]   = progs['crit_probs'][inds] # Probability of developing critical disease
        self.death_prob[:]  = progs['death_probs'][inds] # Probability of death
        self.rel_sus[:]     = progs['sus_ORs'][inds]  # Default susceptibilities
        self.rel_trans[:]   = progs['trans_ORs'][inds] * cvu.sample(**self.pars['beta_dist'], size=len(inds))  # Default transmissibilities, with viral load drawn from a distribution

        return


    def check_alert_accuracy(self):
        '''
        Sorts alerts into true and false. Run before check_alerted() because check_alerted() will set alerted status to False in preparation for the next day.
        '''

        # Get those who were alerted and those who are infected
        watch_owners = self.true('has_watch') # list of inds of people with watches
        bool_mask_alerted = self.alerted[watch_owners]
        bool_mask_exposed = self.exposed[watch_owners]

        tp = np.sum(np.logical_and(bool_mask_alerted, bool_mask_exposed))
        fn = np.sum(np.logical_and(np.logical_not(bool_mask_alerted), bool_mask_exposed))
        tn = np.sum(np.logical_and(np.logical_not(bool_mask_alerted), np.logical_not(bool_mask_exposed)))
        fp = np.sum(np.logical_and(bool_mask_alerted, np.logical_not(bool_mask_exposed)))

        # print("DEBUG: SMARTWATCH ALERT PROBS 2")
        # days = np.arange(-21, 22, 1)
        # for day in days:
        #     # filled conditionally on watch ownership. 
        #     w_inf_day_i = self.all_peoples_infection_days == day
        #     # alerted and has watch and infected on day i. 
        #     w_inf_day_i_alerted = np.logical_and(w_inf_day_i, self.sw_alarmed)

        #     self.sum_people_on_day[day] += np.sum(w_inf_day_i)
        #     self.sum_alert_on_day[day] += np.sum(w_inf_day_i_alerted)
        # #### END DEBUG

        return tp, fn, tn, fp


    def update_states_pre(self, t):
        ''' Perform all state updates at the current timestep '''

        # Initialize
        self.t = t
        self.is_exp = self.true('exposed') # For storing the interim values since used in every subsequent calculation

        # Perform updates
        self.init_flows()
        self.flows['new_infectious']    += self.check_infectious() # For people who are exposed and not infectious, check if they begin being infectious
        self.flows['new_symptomatic']   += self.check_symptomatic()

        if self.pars['enable_multiregion']:
            self.flows['new_severe']        += self.check_severe_mr()
        else:
            self.flows['new_severe']        += self.check_severe()

        self.flows['new_critical']      += self.check_critical()
        self.flows['new_recoveries']    += self.check_recovery()
        new_deaths, new_known_deaths     = self.check_death()
        self.flows['new_deaths']        += new_deaths
        self.flows['new_known_deaths']  += new_known_deaths
        return


    def update_states_post(self):
        ''' Perform post-timestep updates '''

        if self.pars['enable_multiregion']:
            self.flows['new_diagnoses'] += self.check_diagnosed_mr()
        else:
            self.flows['new_diagnoses'] += self.check_diagnosed()
        self.flows['new_quarantined'] += self.check_quar()

        # Take care of smartwatch calculations
        if self.pars['enable_smartwatches']:
            tp, fn, tn, fp = self.check_alert_accuracy()
            self.flows['new_alerts_tp'] += tp
            self.flows['new_alerts_fn'] += fn
            self.flows['new_alerts_tn'] += tn
            self.flows['new_alerts_fp'] += fp
            self.flows['new_alerted'] += self.check_alerted()
            sw_quar, sw_i_quar = self.check_sw_quarantined() # Update the number of smartwatch users quarantined. 
            self.flows['new_Q_w'] += sw_quar
            self.flows['new_Q_w_i'] += sw_i_quar

        del self.is_exp  # Tidy up

        return


    def check_sw_quarantined(self):
        sw_quarantined_arr = np.logical_and(self.quarantined, self.has_watch)
        sw_quarantined_count = np.sum(sw_quarantined_arr)
        sw_correctly_quarantined_arr = np.logical_and(sw_quarantined_arr, self.exposed)
        sw_i_quarantined_count = sw_quarantined_count - np.sum(sw_correctly_quarantined_arr)

        return sw_quarantined_count, sw_i_quarantined_count


    def update_contacts(self):
        ''' Refresh dynamic contacts, e.g. community '''
        # Figure out if anything needs to be done -- e.g. {'h':False, 'c':True}
        for lkey, is_dynam in self.pars['dynam_layer'].items():
            if is_dynam:
                self.contacts[lkey].update(self)

        return self.contacts


    def schedule_behaviour(self, behaviour_pars):
        ''' Schedules quarantines on the basis of results or notifications received today '''

        # PEOPLE WHO BECAME SYMPTOMATIC TODAY
        if self.pars['enable_testobjs']:
            symptomatic_inds = cvu.true((self.date_symptomatic == self.t) | (self.date_symptomatic_ILI == self.t)) # Checks who became symptomatic today
        else:
            symptomatic_inds = cvu.true(self.date_symptomatic == self.t)
        key_to_use = self.get_key_to_use(behaviour_pars['symptom_quar_pars']) # Retrieves the latest policy
        symptomatic_who_quarantine = cvu.binomial_filter(behaviour_pars['symptom_quar_pars'][key_to_use][1], symptomatic_inds)
        self.schedule_quarantine(symptomatic_who_quarantine, start_date=self.t, period=behaviour_pars['symptom_quar_pars'][key_to_use][0]) # Quarantines indices for that long


        # HOUSEHOLD MEMBERS OF PEOPLE WHO BECAME SYMPTOMATIC TODAY
        # Find baseline list of household members
        if self.pars['enable_testobjs']:
            symptomatic_inds = cvu.true((self.date_symptomatic == self.t) | (self.date_symptomatic_ILI == self.t)) # Checks who became symptomatic today
        else:
            symptomatic_inds = cvu.true(self.date_symptomatic == self.t)
        hh_of_symptomatic_inds = self.contacts['h'].find_contacts(symptomatic_inds) # Check who they live with
        
        # Remove exempt household members
        eligible_hh_of_symptomatic_inds = [x for x in hh_of_symptomatic_inds if x not in cvu.true(self.t - self.date_symptomatic < 10)] # Remove household members who already became symptomatic in the last 10 days
        key_to_use = self.get_key_to_use(behaviour_pars['enable_hh_symp_exemption']) # Retrieves the latest policy related to EXEMPTIONS
        if behaviour_pars['enable_hh_symp_exemption'][key_to_use]:
            exempt = cvu.true((self.age < 18) | (self.vaccinated))
            eligible_hh_of_symptomatic_inds = [x for x in eligible_hh_of_symptomatic_inds if x not in exempt]

        # Quarantine eligible household members
        key_to_use = self.get_key_to_use(behaviour_pars['hh_symp_quar_pars']) # Retrieves the latest policy
        hh_of_symptomatic_who_quarantine = cvu.binomial_filter(behaviour_pars['hh_symp_quar_pars'][key_to_use][1], np.array(eligible_hh_of_symptomatic_inds))
        self.schedule_quarantine(hh_of_symptomatic_who_quarantine, start_date=self.t, period=behaviour_pars['hh_symp_quar_pars'][key_to_use][0]) # Quarantines indices for that long
        

        # PEOPLE WHO WERE NOTIFIED TODAY THEY WERE CLOSE CONTACTS
        contacted_inds = cvu.true(self.date_known_contact == self.t) # Checks who was contacted today
        key_to_use = self.get_key_to_use(behaviour_pars['enable_contact_exemption']) # Retrieves the latest policy related to EXEMPTIONS
        if behaviour_pars['enable_contact_exemption'][key_to_use]:
            exempt = cvu.true((self.age < 18) | (self.vaccinated))
            contacted_inds = [x for x in contacted_inds if x not in exempt]
        key_to_use = self.get_key_to_use(behaviour_pars['contact_quar_pars']) # Retrieves the latest policy
        contacted_who_quarantine = cvu.binomial_filter(behaviour_pars['contact_quar_pars'][key_to_use][1], np.array(contacted_inds))
        self.schedule_quarantine(contacted_who_quarantine, start_date=self.t, period=behaviour_pars['contact_quar_pars'][key_to_use][0]) # Quarantines indices for that long


        # PEOPLE WHO RECEIVED SMARTWATCH ALERTS
        if self.pars['enable_smartwatches']:

            # Retrieve the latest policy related to smartwatch behaviour
            key_to_use = self.get_key_to_use(behaviour_pars['smartwatch_behaviour']) # Retrieves the latest policy related to EXEMPTIONS

            # If we have enabled smartwatch behaviour on the present day
            if behaviour_pars['smartwatch_behaviour'][key_to_use][0]:

                # Get the number of alerts the person has received; if standardized, just consider current day; if historical, consider alert histories
                n_alerts = (self.alerted).astype(np.int64)
                if behaviour_pars['smartwatch_behaviour'][key_to_use][1] == 'hist':
                    n_alerts = np.sum(self.alert_histories, axis=1).astype(np.int64) * self.alerted # In the historical case, only count alerts if they received one today

                # Get the duration of quarantine
                quar_len = np.zeros(len(n_alerts))
                if behaviour_pars['smartwatch_behaviour'][key_to_use][1] == 'hist':
                    
                    # If considering historical behaviour, draw quarantine durations based on the number of alerts received
                    for i in range(self.LEN_ALERT_HIST):
                        indices = cvu.true(n_alerts == i+1)
                        n = self.pars['alert_behav_dist_hist'][i+1]['n']
                        a = int(self.pars['alert_behav_dist_hist'][i+1]['s'] * self.pars['alert_behav_dist_hist'][i+1]['p'])
                        b = int(self.pars['alert_behav_dist_hist'][i+1]['s'] * (1-self.pars['alert_behav_dist_hist'][i+1]['p']))
                        quar_len[indices] = stats.betabinom.rvs(n, a, b, size=len(indices))
                        
                else:
                    # If considering standardized behaviour, draw durations from the same distribution all the time
                    indices = cvu.true(n_alerts == 1)
                    n = self.pars['alert_behav_dist_stan']['n']
                    a = int(self.pars['alert_behav_dist_stan']['s'] * self.pars['alert_behav_dist_stan']['p'])
                    b = int(self.pars['alert_behav_dist_stan']['s'] * (1-self.pars['alert_behav_dist_stan']['p']))
                    quar_len[indices] = stats.betabinom.rvs(n, a, b, size=len(indices))

                # Quarantine people accordingly, batching them by number of days to quarantine; currently assumes maximum number of days is 5
                for i in np.arange(0,5):
                    indices = cvu.true(quar_len == i+1)
                    self.schedule_quarantine(indices, start_date=self.t, period=i+1)

        return

    
    def get_key_to_use(self, policy_dict):
        ''' Helper function used to determine what policy to use when scheduling behaviour '''

        keys = np.array(list(policy_dict.keys()))
        key_to_use = keys[cvu.true(keys <= self.t)[-1]]

        return key_to_use


    #%% Methods for updating state

    def check_inds(self, current, date, filter_inds=None):
        ''' Return indices for which the current state is false and which are assigned a date on or before the current date

        Args:
            current (array): list of boolean values that represent a current state
            date (array): list that contains either a date or a Nan
        '''
        if filter_inds is None:
            not_current = cvu.false(current)
        else:
            not_current = cvu.ifalsei(current, filter_inds)
        has_date = cvu.idefinedi(date, not_current)
        inds     = cvu.itrue(self.t >= date[has_date], has_date)
        return inds


    def check_infectious(self):
        ''' Check if they become infectious '''
        inds = self.check_inds(self.infectious, self.date_infectious, filter_inds=self.is_exp)
        self.infectious[inds] = True
        self.infectious_variant[inds] = self.exposed_variant[inds]

        if self.pars['enable_stratifications']:
            for strat, strat_pars in self.pars['stratification_pars'].items():
                metrics = strat_pars['metrics']

                if 'new_infectious' in metrics:
                    bracs = strat_pars['brackets']
                    
                    # update the value for the current day, for each bracket. 
                    for brac in bracs:
                        brac_name = "_".join([ str(brac[0]), str(brac[1])])

                        # Count the number of individuals meeting the criteria. 
                        if strat == 'age':
                            num_inds = np.sum( (self.age[inds] >= brac[0]) & (self.age[inds] < brac[1]) )
                        elif strat == 'income':
                            num_inds = np.sum( (self.income[inds] >= brac[0]) & (self.income[inds] < brac[1]) )
                        else:
                            raise ValueError(f"Stratification {strat} not recognized.")

                        self.stratifications[strat][brac_name]['new_infectious'][self.t] += num_inds

        for variant in range(self.pars['n_variants']):
            this_variant_inds = cvu.itrue(self.infectious_variant[inds] == variant, inds)
            n_this_variant_inds = len(this_variant_inds)
            self.flows_variant['new_infectious_by_variant'][variant] += n_this_variant_inds
            self.infectious_by_variant[variant, this_variant_inds] = True
        return len(inds)


    def check_symptomatic(self):
        ''' Check for new progressions to symptomatic '''
        inds = self.check_inds(self.symptomatic, self.date_symptomatic, filter_inds=self.is_exp)
        self.symptomatic[inds] = True
        return len(inds)

    def check_severe_mr(self):
        ''' Check for new progressions to severe. Update flows for each region. '''
        # Main result
        inds = self.check_inds(self.severe, self.date_severe, filter_inds=self.is_exp)
        
        #  STRATIFICATION
        if self.pars['enable_stratifications']:
            for strat, strat_pars in self.pars['stratification_pars'].items():
                metrics = strat_pars['metrics']

                if 'new_severe' in metrics:
                    bracs = strat_pars['brackets']
                    
                    # update the value for the current day, for each bracket.
                    for brac in bracs:
                        brac_name = "_".join([ str(brac[0]), str(brac[1])])

                        # Count the number of individuals meeting the criteria. 
                        if strat == 'age':
                            num_inds = np.sum( (self.age[inds] >= brac[0]) & (self.age[inds] < brac[1]) )
                        elif strat == 'income':
                            num_inds = np.sum( (self.income[inds] >= brac[0]) & (self.income[inds] < brac[1]) )
                        else:
                            raise ValueError(f"Stratification {strat} not recognized.")

                        self.stratifications[strat][brac_name]['new_severe'][self.t] += num_inds

        rnames = self.pars['multiregion']['rnames']
        rsizes = self.pars['multiregion']['rsizes']
        rstarts = self.pars['multiregion']['rstarts']

        for rname, rstart, rsize in zip(rnames, rstarts, rsizes):
            # Count the # indicies in the region. i.e. indicies greater than rstart and less than the end .
            num_inds = np.sum((inds >= rstart) & (inds < rstart+rsize))
            self.flows[f'{rname}_new_severe'] += num_inds

        self.severe[inds] = True

        # Assert that the sum over the regions = len(inds) -> Works. 
        # assert np.sum(self.flows[f'{rname}_new_severe'] for rname in rnames) == len(inds)

        return len(inds)

    def check_severe(self):
        ''' Check for new progressions to severe '''
        inds = self.check_inds(self.severe, self.date_severe, filter_inds=self.is_exp)
        
        #  STRATIFICATION
        if self.pars['enable_stratifications']:
            for strat, strat_pars in self.pars['stratification_pars'].items():
                metrics = strat_pars['metrics']

                if 'new_severe' in metrics:
                    bracs = strat_pars['brackets']
                    
                    # update the value for the current day, for each bracket. 
                    for brac in bracs:
                        brac_name = "_".join([ str(brac[0]), str(brac[1])])

                        # Count the number of individuals meeting the criteria. 
                        if strat == 'age':
                            num_inds = np.sum( (self.age[inds] >= brac[0]) & (self.age[inds] < brac[1]) )
                        elif strat == 'income':
                            num_inds = np.sum( (self.income[inds] >= brac[0]) & (self.income[inds] < brac[1]) )
                        else:
                            raise ValueError(f"Stratification {strat} not recognized.")

                        self.stratifications[strat][brac_name]['new_severe'][self.t] += num_inds
        
        self.severe[inds] = True
        return len(inds)


    def check_critical(self):
        ''' Check for new progressions to critical '''
        inds = self.check_inds(self.critical, self.date_critical, filter_inds=self.is_exp)
        self.critical[inds] = True
        return len(inds)


    def check_recovery(self, inds=None, filter_inds='is_exp'):
        '''
        Check for recovery.

        More complex than other functions to allow for recovery to be manually imposed
        for a specified set of indices.
        '''

        # Handle more flexible options for setting indices
        if filter_inds == 'is_exp':
            filter_inds = self.is_exp
        if inds is None:
            inds = self.check_inds(self.recovered, self.date_recovered, filter_inds=filter_inds)

        # Now reset all disease states
        self.exposed[inds]          = False
        self.infectious[inds]       = False
        self.symptomatic[inds]      = False
        self.severe[inds]           = False
        self.critical[inds]         = False
        self.recovered[inds]        = True
        self.recovered_variant[inds] = self.exposed_variant[inds]
        self.infectious_variant[inds] = np.nan
        self.exposed_variant[inds]    = np.nan
        self.exposed_by_variant[:, inds] = False
        self.infectious_by_variant[:, inds] = False

        # Handle immunity aspects
        if self.pars['use_waning']:

            # Reset additional states
            self.susceptible[inds] = True
            self.diagnosed[inds]   = False # Reset their diagnosis state because they might be reinfected

        # # Handle instances of false positive diagnoses; if someone is not exposed, is diagnosed, and has passed the required time in isolation, un-diagnose them
        # fp_diagnoses = cvu.true((~self.exposed) & (self.diagnosed) & (self.t - self.date_diagnosed >= 14))
        # self.diagnosed[fp_diagnoses] = False

        return len(inds)


    def check_death(self):
        ''' Check whether or not this person died on this timestep  '''
        inds = self.check_inds(self.dead, self.date_dead, filter_inds=self.is_exp)
    
        #  STRATIFICATION
        if self.pars['enable_stratifications']:
            for strat, strat_pars in self.pars['stratification_pars'].items():
                metrics = strat_pars['metrics']

                if 'new_deaths' in metrics:
                    bracs = strat_pars['brackets']
                    
                    # update the value for the current day, for each bracket. 
                    for brac in bracs:
                        brac_name = "_".join([ str(brac[0]), str(brac[1])])

                        # Count the number of individuals meeting the criteria. 
                        if strat == 'age':
                            num_inds = np.sum( (self.age[inds] >= brac[0]) & (self.age[inds] < brac[1]) )
                        elif strat == 'income':
                            num_inds = np.sum( (self.income[inds] >= brac[0]) & (self.income[inds] < brac[1]) )
                        else:
                            raise ValueError(f"Stratification {strat} not recognized.")

                        self.stratifications[strat][brac_name]['new_deaths'][self.t] += num_inds        
        
        self.dead[inds]             = True
        diag_inds = inds[self.diagnosed[inds]] # Check whether the person was diagnosed before dying
        self.known_dead[diag_inds]  = True
        self.susceptible[inds]      = False
        self.exposed[inds]          = False
        self.infectious[inds]       = False
        self.symptomatic[inds]      = False
        self.severe[inds]           = False
        self.critical[inds]         = False
        self.known_contact[inds]    = False
        self.quarantined[inds]      = False
        self.recovered[inds]        = False
        self.infectious_variant[inds] = np.nan
        self.exposed_variant[inds]    = np.nan
        self.recovered_variant[inds]  = np.nan

        return len(inds), len(diag_inds)


    def check_alerted(self):
        '''
        Check which people received alerts on this timestep and whether they were correct or incorrect.
        Reset alerted status to False for everyone to prepare for the next day.
        '''

        # Check who was alerted
        n_alerted = len(cvu.true(self.alerted))

        # Set the alerted status of everyone to false
        self.alerted[:] = False

        return n_alerted


    def check_diagnosed_mr(self):
        '''
        MULTI-REGION VERSION
        Check for new diagnoses. Since most data are reported with diagnoses on
        the date of the test, this function reports counts not for the number of
        people who received a positive test result on a day, but rather, the number
        of people who were tested on that day who are schedule to be diagnosed in
        the future.
        '''

        #### THIS DRIVES THE OUTPUT SO I'LL EDIT THIS TO MULTIREGION. 
        # Handle people who tested today who will be diagnosed in future (i.e., configure them to have have finished being tested)
        test_pos_inds = self.check_inds(self.diagnosed, self.date_pos_test, filter_inds=None) # Find people who are not diagnosed and have a date of a positive test that is today or earlier

        # Update per variant. (Written 31-05-23)
        for variant in range(self.pars['n_variants']):
            # Get the people who are exposed for this variant, and have a positive test. 
            # TODO: exposed_variant ok? 
            var_count = np.sum(self.exposed_by_variant[variant, test_pos_inds])
            self.flows_variant['new_diagnoses_by_variant'][variant] += var_count

        #### REGIONAL RESULTS
        rnames = self.pars['multiregion']['rnames']
        rsizes = self.pars['multiregion']['rsizes']
        rstarts = self.pars['multiregion']['rstarts']

        for rname, rstart, rsize in zip(rnames, rstarts, rsizes):
            # Count the # indicies in the region. i.e. indicies greater than rstart and less than the end .
            num_inds = np.sum((test_pos_inds >= rstart) & (test_pos_inds < rstart+rsize))
            self.flows[f'{rname}_new_diagnoses'] += num_inds

        self.date_pos_test[test_pos_inds] = np.nan # Clear date of having will-be-positive test
        
        inds = test_pos_inds
        #### OTHER STATIFICATIONS
        if self.pars['enable_stratifications']:
            for strat, strat_pars in self.pars['stratification_pars'].items():
                metrics = strat_pars['metrics']

                if 'new_diagnoses' in metrics:
                    bracs = strat_pars['brackets']
                    
                    # update the value for the current day, for each bracket. 
                    for brac in bracs:
                        brac_name = "_".join([ str(brac[0]), str(brac[1])])

                        # Count the number of individuals meeting the criteria. 
                        if strat == 'age':
                            num_inds = np.sum( (self.age[inds] >= brac[0]) & (self.age[inds] < brac[1]) )
                        elif strat == 'income':
                            num_inds = np.sum( (self.income[inds] >= brac[0]) & (self.income[inds] < brac[1]) )
                        else:
                            raise ValueError(f"Stratification {strat} not recognized.")

                        self.stratifications[strat][brac_name]['new_diagnoses'][self.t] += num_inds

        #### WEIRD DEFAULT COVASIM CODE
        # Handle people who were actually diagnosed today (i.e., set them as diagnosed and remove any of them that were quarantining from quarantine)
        diag_inds  = self.check_inds(self.diagnosed, self.date_diagnosed, filter_inds=None) # Find who are not diagnosed and have a date of diagnosis that is today or earlier
        self.diagnosed[diag_inds]   = True # Set these people to be diagnosed
        quarantined = cvu.itruei(self.quarantined, diag_inds) # Find individuals who were just diagnosed who are in quarantine
        self.date_end_quarantine[quarantined] = self.t # Set end quarantine date to match when the person left quarantine (and entered isolation)
        self.quarantined[diag_inds] = False # If you are diagnosed, you are isolated, not in quarantine

        return len(test_pos_inds)

    def check_diagnosed(self):
        '''
        Check for new diagnoses. Since most data are reported with diagnoses on
        the date of the test, this function reports counts not for the number of
        people who received a positive test result on a day, but rather, the number
        of people who were tested on that day who are schedule to be diagnosed in
        the future.
        '''

        # Handle people who tested today who will be diagnosed in future (i.e., configure them to have have finished being tested)
        test_pos_inds = self.check_inds(self.diagnosed, self.date_pos_test, filter_inds=None) # Find people who are not diagnosed and have a date of a positive test that is today or earlier
        self.date_pos_test[test_pos_inds] = np.nan # Clear date of having will-be-positive test

        # Update per variant. (Written 31-05-23)
        for variant in range(self.pars['n_variants']):
            # Get the people who are exposed for this variant, and have a positive test. 
            # TODO: exposed_variant ok? 
            var_count = np.sum(self.exposed_by_variant[variant, test_pos_inds])
            self.flows_variant['new_diagnoses_by_variant'][variant] += var_count

        # STRATIFICATIONS
        inds = test_pos_inds
        if self.pars['enable_stratifications']:
            for strat, strat_pars in self.pars['stratification_pars'].items():
                metrics = strat_pars['metrics']

                if 'new_diagnoses' in metrics:
                    bracs = strat_pars['brackets']
                    
                    # update the value for the current day, for each bracket. 
                    for brac in bracs:
                        brac_name = "_".join([ str(brac[0]), str(brac[1])])

                        # Count the number of individuals meeting the criteria. 
                        if strat == 'age':
                            num_inds = np.sum( (self.age[inds] >= brac[0]) & (self.age[inds] < brac[1]) )
                        elif strat == 'income':
                            num_inds = np.sum( (self.income[inds] >= brac[0]) & (self.income[inds] < brac[1]) )
                        else:
                            raise ValueError(f"Stratification {strat} not recognized.")

                        self.stratifications[strat][brac_name]['new_diagnoses'][self.t] += num_inds

        # Handle people who were actually diagnosed today (i.e., set them as diagnosed and remove any of them that were quarantining from quarantine)
        diag_inds  = self.check_inds(self.diagnosed, self.date_diagnosed, filter_inds=None) # Find who are not diagnosed and have a date of diagnosis that is today or earlier
        self.diagnosed[diag_inds]   = True # Set these people to be diagnosed
        quarantined = cvu.itruei(self.quarantined, diag_inds) # Find individuals who were just diagnosed who are in quarantine
        self.date_end_quarantine[quarantined] = self.t # Set end quarantine date to match when the person left quarantine (and entered isolation)
        self.quarantined[diag_inds] = False # If you are diagnosed, you are isolated, not in quarantine

        return len(test_pos_inds)


    def check_quar(self):
        ''' Update quarantine state '''

        n_quarantined = 0 # Number of people entering quarantine
        for ind,end_day in self._pending_quarantine[self.t]:
            if self.quarantined[ind]: # Update when quarantine should be finished (in case schedule_quarantine is called on someone already in quarantine)
                self.date_end_quarantine[ind] = max(self.date_end_quarantine[ind], end_day) # Extend quarantine if required
            elif not (self.dead[ind] or self.recovered[ind] or self.diagnosed[ind]): # Unclear whether recovered should be included here # elif not (self.dead[ind] or self.diagnosed[ind]):
                self.quarantined[ind] = True
                self.date_quarantined[ind] = self.t
                self.date_end_quarantine[ind] = end_day
                n_quarantined += 1

        # If someone on quarantine has reached the end of their quarantine, release them
        end_inds = self.check_inds(~self.quarantined, self.date_end_quarantine, filter_inds=None) # Note the double-negative here (~)
        self.quarantined[end_inds] = False # Release from quarantine

        # Remove people from quarantine if enable_behaviour is True and if the right conditions are met
        if self.pars['enable_behaviour']:
            
            # Remove people based on consecutive negative RAT
            key_to_use = self.get_key_to_use(self.pars['behaviour_pars']['enable_symp_early_end']) # Retrieves the latest policy on ending quarantine early
            if self.pars['behaviour_pars']['enable_symp_early_end'][key_to_use]:
                exempt = cvu.true(self.cons_days_neg_rat >= 2)
                self.quarantined[exempt] = False
                self.date_end_quarantine[exempt] = self.t

        # Update the counter for consecutive days in quarantine
        self.cons_days_in_quar[self.quarantined] += 1
        self.cons_days_in_quar[~self.quarantined] = 0

        return n_quarantined


    #%% Methods to make events occur (infection and diagnosis)

    def make_naive(self, inds, reset_vx=False):
        '''
        Make a set of people naive. This is used during dynamic resampling.

        Args:
            inds (array): list of people to make naive
            reset_vx (bool): whether to reset vaccine-derived immunity
        '''
        for key in self.meta.states:
            if key in ['susceptible', 'naive']:
                self[key][inds] = True
            else:
                if (key != 'vaccinated') or reset_vx: # Don't necessarily reset vaccination
                    self[key][inds] = False

        # Reset variant states
        for key in self.meta.variant_states:
            self[key][inds] = np.nan
        for key in self.meta.by_variant_states:
            self[key][:, inds] = False

        # Reset immunity and antibody states
        non_vx_inds = inds if reset_vx else inds[~self['vaccinated'][inds]]
        for key in self.meta.imm_states:
            self[key][:, non_vx_inds] = 0
        for key in self.meta.nab_states + self.meta.vacc_states:
            self[key][non_vx_inds] = 0

        # Reset dates
        for key in self.meta.dates + self.meta.durs:
            if (key != 'date_vaccinated') or reset_vx: # Don't necessarily reset vaccination
                self[key][inds] = np.nan

        return


    def make_nonnaive(self, inds, set_recovered=False, date_recovered=0):
        '''
        Make a set of people non-naive.

        This can be done either by setting only susceptible and naive states,
        or else by setting them as if they have been infected and recovered.
        '''
        self.make_naive(inds) # First make them naive and reset all other states

        # Make them non-naive
        for key in ['susceptible', 'naive']:
            self[key][inds] = False

        if set_recovered:
            self.date_recovered[inds] = date_recovered # Reset date recovered
            self.check_recovered(inds=inds, filter_inds=None) # Set recovered

        return



    def infect(self, inds, hosp_max=None, icu_max=None, source=None, layer=None, variant=0):
        '''
        Infect people and determine their eventual outcomes.

            * Every infected person can infect other people, regardless of whether they develop symptoms
            * Infected people that develop symptoms are disaggregated into mild vs. severe (=requires hospitalization) vs. critical (=requires ICU)
            * Every asymptomatic, mildly symptomatic, and severely symptomatic person recovers
            * Critical cases either recover or die
            * If the simulation is being run with waning, this method also sets/updates agents' neutralizing antibody levels

        Method also deduplicates input arrays in case one agent is infected many times
        and stores who infected whom in infection_log list.

        Args:
            inds     (array): array of people to infect
            hosp_max (bool):  whether or not there is an acute bed available for this person
            icu_max  (bool):  whether or not there is an ICU bed available for this person
            source   (array): source indices of the people who transmitted this infection (None if an importation or seed infection)
            layer    (str):   contact layer this infection was transmitted on
            variant  (int):   the variant people are being infected by

        Returns:
            count (int): number of people infected
        '''

        if len(inds) == 0:
            return 0

        # Remove duplicates
        inds, unique = np.unique(inds, return_index=True)
        if source is not None:
            source = source[unique]

        # Keep only susceptibles
        keep = self.susceptible[inds] # Unique indices in inds and source that are also susceptible
        inds = inds[keep]
        if source is not None:
            source = source[keep]

        # Deal with variant parameters
        variant_keys = ['rel_symp_prob', 'rel_severe_prob', 'rel_crit_prob', 'rel_death_prob']
        infect_pars = {k:self.pars[k] for k in variant_keys}
        variant_label = self.pars['variant_map'][variant]
        if variant:
            for k in variant_keys:
                infect_pars[k] *= self.pars['variant_pars'][variant_label][k]

        n_infections = len(inds)
        durpars      = self.pars['dur']

        # Retrieve those with a breakthrough infection (defined nabs)
        breakthrough_inds = inds[cvu.true(self.peak_nab[inds])]
        if len(breakthrough_inds):
            no_prior_breakthrough = (self.n_breakthroughs[breakthrough_inds] == 0) # We only adjust transmissibility for the first breakthrough
            new_breakthrough_inds = breakthrough_inds[no_prior_breakthrough]
            self.rel_trans[new_breakthrough_inds] *= self.pars['trans_redux']

        # Update states, variant info, and flows
        self.susceptible[inds]    = False
        self.naive[inds]          = False
        self.recovered[inds]      = False
        self.diagnosed[inds]      = False
        self.exposed[inds]        = True
        self.n_infections[inds]  += 1
        self.n_breakthroughs[breakthrough_inds] += 1
        self.exposed_variant[inds] = variant
        self.exposed_by_variant[variant, inds] = True
        self.flows['new_infections']   += len(inds)
        self.flows['new_reinfections'] += len(cvu.defined(self.date_recovered[inds])) # Record reinfections
        self.flows_variant['new_infections_by_variant'][variant] += len(inds)

        # Record transmissions
        for i, target in enumerate(inds):
            entry = dict(source=source[i] if source is not None else None, target=target, date=self.t, layer=layer, variant=variant_label)
            self.infection_log.append(entry)

        # Calculate how long before this person can infect other people
        self.dur_exp2inf[inds] = cvu.sample(**durpars['exp2inf'], size=n_infections)
        self.date_exposed[inds] = self.t
        self.date_infectious[inds] = self.dur_exp2inf[inds] + self.t

        # Reset all other dates
        for key in ['date_symptomatic', 'date_severe', 'date_critical', 'date_diagnosed', 'date_recovered']:
            self[key][inds] = np.nan

        # Use prognosis probabilities to determine what happens to them
        symp_probs = infect_pars['rel_symp_prob']*self.symp_prob[inds]*(1-self.symp_imm[variant, inds]) # Calculate their actual probability of being symptomatic
        is_symp = cvu.binomial_arr(symp_probs) # Determine if they develop symptoms
        symp_inds = inds[is_symp]
        asymp_inds = inds[~is_symp] # Asymptomatic
        self.flows_variant['new_symptomatic_by_variant'][variant] += len(symp_inds)

        # CASE 1: Asymptomatic: may infect others, but have no symptoms and do not die
        dur_asym2rec = cvu.sample(**durpars['asym2rec'], size=len(asymp_inds))
        self.date_recovered[asymp_inds] = self.date_infectious[asymp_inds] + dur_asym2rec  # Date they recover. Nx1.
        self.dur_disease[asymp_inds] = self.dur_exp2inf[asymp_inds] + dur_asym2rec  # Store how long this person had COVID-19. Nx1. 

        # CASE 2: Symptomatic: can either be mild, severe, or critical
        n_symp_inds = len(symp_inds)
        self.dur_inf2sym[symp_inds] = cvu.sample(**durpars['inf2sym'], size=n_symp_inds) # Store how long this person took to develop symptoms
        self.date_symptomatic[symp_inds] = self.date_infectious[symp_inds] + self.dur_inf2sym[symp_inds] # Date they become symptomatic
        sev_probs = infect_pars['rel_severe_prob'] * self.severe_prob[symp_inds]*(1-self.sev_imm[variant, symp_inds]) # Probability of these people being severe
        is_sev = cvu.binomial_arr(sev_probs) # See if they're a severe or mild case
        sev_inds = symp_inds[is_sev]
        mild_inds = symp_inds[~is_sev] # Not severe
        self.flows_variant['new_severe_by_variant'][variant] += len(sev_inds)

        # CASE 2.1: Mild symptoms, no hospitalization required and no probability of death
        dur_mild2rec = cvu.sample(**durpars['mild2rec'], size=len(mild_inds))
        self.date_recovered[mild_inds] = self.date_symptomatic[mild_inds] + dur_mild2rec  # Date they recover
        self.dur_disease[mild_inds] = self.dur_exp2inf[mild_inds] + self.dur_inf2sym[mild_inds] + dur_mild2rec  # Store how long this person had COVID-19

        # CASE 2.2: Severe cases: hospitalization required, may become critical
        self.dur_sym2sev[sev_inds] = cvu.sample(**durpars['sym2sev'], size=len(sev_inds)) # Store how long this person took to develop severe symptoms
        self.date_severe[sev_inds] = self.date_symptomatic[sev_inds] + self.dur_sym2sev[sev_inds]  # Date symptoms become severe
        crit_probs = infect_pars['rel_crit_prob'] * self.crit_prob[sev_inds] * (self.pars['no_hosp_factor'] if hosp_max else 1.) # Probability of these people becoming critical - higher if no beds available
        is_crit = cvu.binomial_arr(crit_probs)  # See if they're a critical case
        crit_inds = sev_inds[is_crit]
        non_crit_inds = sev_inds[~is_crit]

        # CASE 2.2.1 Not critical - they will recover
        dur_sev2rec = cvu.sample(**durpars['sev2rec'], size=len(non_crit_inds))
        self.date_recovered[non_crit_inds] = self.date_severe[non_crit_inds] + dur_sev2rec  # Date they recover
        self.dur_disease[non_crit_inds] = self.dur_exp2inf[non_crit_inds] + self.dur_inf2sym[non_crit_inds] + self.dur_sym2sev[non_crit_inds] + dur_sev2rec  # Store how long this person had COVID-19

        # CASE 2.2.2: Critical cases: ICU required, may die
        self.dur_sev2crit[crit_inds] = cvu.sample(**durpars['sev2crit'], size=len(crit_inds))
        self.date_critical[crit_inds] = self.date_severe[crit_inds] + self.dur_sev2crit[crit_inds]  # Date they become critical
        death_probs = infect_pars['rel_death_prob'] * self.death_prob[crit_inds] * (self.pars['no_icu_factor'] if icu_max else 1.)# Probability they'll die
        is_dead = cvu.binomial_arr(death_probs)  # Death outcome
        dead_inds = crit_inds[is_dead]
        alive_inds = crit_inds[~is_dead]

        # CASE 2.2.2.1: Did not die
        dur_crit2rec = cvu.sample(**durpars['crit2rec'], size=len(alive_inds))
        self.date_recovered[alive_inds] = self.date_critical[alive_inds] + dur_crit2rec # Date they recover
        self.dur_disease[alive_inds] = self.dur_exp2inf[alive_inds] + self.dur_inf2sym[alive_inds] + self.dur_sym2sev[alive_inds] + self.dur_sev2crit[alive_inds] + dur_crit2rec  # Store how long this person had COVID-19

        # CASE 2.2.2.2: Did die
        dur_crit2die = cvu.sample(**durpars['crit2die'], size=len(dead_inds))
        self.date_dead[dead_inds] = self.date_critical[dead_inds] + dur_crit2die # Date of death
        self.dur_disease[dead_inds] = self.dur_exp2inf[dead_inds] + self.dur_inf2sym[dead_inds] + self.dur_sym2sev[dead_inds] + self.dur_sev2crit[dead_inds] + dur_crit2die   # Store how long this person had COVID-19
        self.date_recovered[dead_inds] = np.nan # If they did die, remove them from recovered

        # HANDLE VIRAL LOAD CONTROL POINTS
        if self.pars['enable_vl']:
            # Get P_inf: where viral load crosses 10^6 cp/mL
            self.x_p_inf[inds] = self.dur_exp2inf[inds]
            self.y_p_inf[inds] = 6

            # Get P1: where viral load crosses 10^3 cp/mL; time difference obtained empirically through simulation
            self.x_p1[inds] = np.maximum(self.x_p_inf[inds] - (np.random.gamma(2, 0.35, size=len(inds)) + 0.25), 0)
            self.y_p1[inds] = 3

            # Get P2: where viral load peaks; time difference obtained empirically through simulation
            self.x_p2[inds] = self.x_p_inf[inds] + (np.random.gamma(3, 0.26, size=len(inds)) + 0.1)
            self.y_p2[inds] = ((self.y_p_inf[inds] - self.y_p1[inds])*(self.x_p2[inds] - self.x_p1[inds])/(self.x_p_inf[inds] - self.x_p1[inds])) + self.y_p1[inds]

            # Align P1, P_inf, and P2 to current time
            self.x_p1[inds] = self.x_p1[inds] + self.t
            self.x_p_inf[inds] = self.x_p_inf[inds] + self.t
            self.x_p2[inds] = self.x_p2[inds] + self.t

            # Get P3: where viral load drops below 10^6 cp/mL
            time_recovered = np.ones(len(self.date_recovered), dtype=cvd.default_float)*self.date_recovered # This is needed to make a copy
            inds_dead = ~np.isnan(self.date_dead)
            time_recovered[inds_dead] = self.date_dead[inds_dead]
            self.x_p3[inds] = np.maximum(time_recovered[inds], self.x_p2[inds])
            self.y_p3[inds] = 6

            # # For testing purposes
            # if self.t < self.pars['x_p1'].shape[1]:
            #     self.pars['x_p1'][:, self.t] = self.x_p1
            #     self.pars['x_p_inf'][:, self.t] = self.x_p_inf
            #     self.pars['x_p2'][:, self.t] = self.x_p2
            #     self.pars['x_p3'][:, self.t] = self.x_p3
            #     self.pars['y_p1'][:, self.t] = self.y_p1
            #     self.pars['y_p_inf'][:, self.t] = self.y_p_inf
            #     self.pars['y_p2'][:, self.t] = self.y_p2
            #     self.pars['y_p3'][:, self.t] = self.y_p3

        # Handle immunity aspects
        if self.pars['use_waning']:
            symp = dict(asymp=asymp_inds, mild=mild_inds, sev=sev_inds)
            cvi.update_peak_nab(self, inds, nab_pars=self.pars, symp=symp)

        return n_infections # For incrementing counters


    def test(self, inds, test_sensitivity=1.0, loss_prob=0.0, test_delay=0):
        '''
        Method to test people. Typically not to be called by the user directly;
        see the test_num() and test_prob() interventions.

        Args:
            inds: indices of who to test
            test_sensitivity (float): probability of a true positive
            loss_prob (float): probability of loss to follow-up
            test_delay (int): number of days before test results are ready
        '''

        inds = np.unique(inds)
        self.tested[inds] = True
        self.date_tested[inds] = self.t # Only keep the last time they tested

        is_infectious = cvu.itruei(self.infectious, inds)
        pos_test      = cvu.n_binomial(test_sensitivity, len(is_infectious))
        is_inf_pos    = is_infectious[pos_test]

        not_diagnosed = is_inf_pos[np.isnan(self.date_diagnosed[is_inf_pos])]
        not_lost      = cvu.n_binomial(1.0-loss_prob, len(not_diagnosed))
        final_inds    = not_diagnosed[not_lost]

        # Store the date the person will be diagnosed, as well as the date they took the test which will come back positive
        self.date_diagnosed[final_inds] = self.t + test_delay
        self.date_pos_test[final_inds] = self.t

        return final_inds


    def schedule_quarantine(self, inds, start_date=None, period=None):
        '''
        Schedule a quarantine. Typically not called by the user directly except
        via a custom intervention; see the contact_tracing() intervention instead.

        This function will create a request to quarantine a person on the start_date for
        a period of time. Whether they are on an existing quarantine that gets extended, or
        whether they are no longer eligible for quarantine, will be checked when the start_date
        is reached.

        Args:
            inds (int): indices of who to quarantine, specified by check_quar()
            start_date (int): day to begin quarantine (defaults to the current day, `sim.t`)
            period (int): quarantine duration (defaults to ``pars['quar_period']``)
        '''

        start_date = self.t if start_date is None else int(start_date)
        period = self.pars['quar_period'] if period is None else int(period)
        for ind in inds:
            self._pending_quarantine[start_date].append((ind, start_date + period))
        return


    #%% Analysis methods

    def plot(self, *args, **kwargs):
        '''
        Plot statistics of the population -- age distribution, numbers of contacts,
        and overall weight of contacts (number of contacts multiplied by beta per
        layer).

        Args:
            bins      (arr)   : age bins to use (default, 0-100 in one-year bins)
            width     (float) : bar width
            font_size (float) : size of font
            alpha     (float) : transparency of the plots
            fig_args  (dict)  : passed to pl.figure()
            axis_args (dict)  : passed to pl.subplots_adjust()
            plot_args (dict)  : passed to pl.plot()
            do_show   (bool)  : whether to show the plot
            fig       (fig)   : handle of existing figure to plot into
        '''
        fig = cvplt.plot_people(people=self, *args, **kwargs)
        return fig


    def story(self, uid, *args):
        '''
        Print out a short history of events in the life of the specified individual.

        Args:
            uid (int/list): the person or people whose story is being regaled
            args (list): these people will tell their stories too

        **Example**::

            sim = cv.Sim(pop_type='hybrid', verbose=0)
            sim.run()
            sim.people.story(12)
            sim.people.story(795)
        '''

        def label_lkey(lkey):
            ''' Friendly name for common layer keys '''
            if lkey.lower() == 'a':
                llabel = 'default contact'
            if lkey.lower() == 'h':
                llabel = 'household'
            elif lkey.lower() == 's':
                llabel = 'school'
            elif lkey.lower() == 'w':
                llabel = 'workplace'
            elif lkey.lower() == 'c':
                llabel = 'community'
            else:
                llabel = f'"{lkey}"'
            return llabel

        uids = sc.promotetolist(uid)
        uids.extend(args)

        for uid in uids:

            p = self[uid]
            sex = 'female' if p.sex == 0 else 'male'

            intro = f'\nThis is the story of {uid}, a {p.age:.0f} year old {sex}'

            if not p.susceptible:
                if np.isnan(p.date_symptomatic):
                    print(f'{intro}, who had asymptomatic COVID.')
                else:
                    print(f'{intro}, who had symptomatic COVID.')
            else:
                print(f'{intro}, who did not contract COVID.')

            total_contacts = 0
            no_contacts = []
            for lkey in p.contacts.keys():
                llabel = label_lkey(lkey)
                n_contacts = len(p.contacts[lkey])
                total_contacts += n_contacts
                if n_contacts:
                    print(f'{uid} is connected to {n_contacts} people in the {llabel} layer')
                else:
                    no_contacts.append(llabel)
            if len(no_contacts):
                nc_string = ', '.join(no_contacts)
                print(f'{uid} has no contacts in the {nc_string} layer(s)')
            print(f'{uid} has {total_contacts} contacts in total')

            events = []

            dates = {
                'date_critical'       : 'became critically ill and needed ICU care',
                'date_dead'           : 'died ☹',
                'date_diagnosed'      : 'was diagnosed with COVID',
                'date_end_quarantine' : 'ended quarantine',
                'date_infectious'     : 'became infectious',
                'date_known_contact'  : 'was notified they may have been exposed to COVID',
                'date_pos_test'       : 'took a positive test',
                'date_quarantined'    : 'entered quarantine',
                'date_recovered'      : 'recovered',
                'date_severe'         : 'developed severe symptoms and needed hospitalization',
                'date_symptomatic'    : 'became symptomatic',
                'date_tested'         : 'was tested for COVID',
                'date_vaccinated'     : 'was vaccinated against COVID',
            }

            for attribute, message in dates.items():
                date = getattr(p,attribute)
                if not np.isnan(date):
                    events.append((date, message))

            for infection in self.infection_log:
                lkey = infection['layer']
                llabel = label_lkey(lkey)
                if infection['target'] == uid:
                    if lkey:
                        events.append((infection['date'], f'was infected with COVID by {infection["source"]} via the {llabel} layer'))
                    else:
                        events.append((infection['date'], 'was infected with COVID as a seed infection'))

                if infection['source'] == uid:
                    x = len([a for a in self.infection_log if a['source'] == infection['target']])
                    events.append((infection['date'],f'gave COVID to {infection["target"]} via the {llabel} layer ({x} secondary infections)'))

            if len(events):
                for day, event in sorted(events, key=lambda x: x[0]):
                    print(f'On day {day:.0f}, {uid} {event}')
            else:
                print(f'Nothing happened to {uid} during the simulation.')
        return

