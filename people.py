'''
Defines the People class and functions associated with making people and handling
the transitions between states (e.g., from susceptible to infected).
'''

#%% Imports
from array import array
from heapq import merge
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
from . import stratify as strat
from . import pathogen_interactions as pint
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
        self.sim = None

        #For IgG conversion error
        self.IgG_conversion_slope = 0

        variantNums = []
        for p in self.pars['pathogens']:
            variantNums.append(len(p.variants)+1)
        max_variants = np.max(variantNums)

        # Set person properties -- all floats except for UID
        for key in self.meta.person:
            if key == 'uid':
                self[key] = np.arange(self.pars['pop_size'], dtype=cvd.default_int)
            elif key in ['n_infections', 'n_breakthroughs']:
                self[key] = np.zeros((self.pars['n_pathogens'],self.pars['pop_size']), dtype=cvd.default_int)
            elif key in ['cons_days_in_quar', 'cons_days_neg_rat']:
                self[key] = np.zeros(self.pars['pop_size'], dtype=cvd.default_int)
            elif key in ['has_watch']:
                self[key] = np.full(self.pars['pop_size'], False, dtype=bool)
            elif key in ['income']:
                self[key] = np.zeros(self.pars['pop_size'], dtype=cvd.default_int)
            elif key in ['viral_load']:
                self[key] = np.zeros((self.pars['n_pathogens'],self.pars['pop_size']), dtype=cvd.default_float)
            elif key in ['is_coinfected']:
                self[key] = np.full(self.pars['pop_size'], False,dtype=bool)
            elif key in ['symp_prob', 'severe_prob',  'crit_prob', 'death_prob','rel_trans', 'rel_sus', 'abs_symp_prob', 'abs_severe_prob', 'abs_crit_prob', 'abs_death_prob']:  
                self[key] = np.zeros((self.pars['n_pathogens'],self.pars['pop_size']), dtype=cvd.default_float)
            else:
                self[key] = np.full(self.pars['pop_size'], np.nan, dtype=cvd.default_float)
        
        self.provides_sample_prob = np.zeros(self.pars['pop_size'], dtype=cvd.default_float) #Keeps track of population sampling 
        
        self.IgG_level = np.zeros(self.pars['pop_size'], dtype=cvd.default_float) #Tracks IgG levels used for populaiton active surveillance

        # Set health states -- only susceptible is true by default -- booleans except exposed by variant which should return the variant that ind is exposed to   
        for key in self.meta.states:
            val = (key in ['susceptible', 'naive']) # Default value is True for susceptible and naive, False otherwise
            self[key] = np.full(self.pars['pop_size'], val, dtype=bool)
             
        #Set pathogen specific states    self[key][pathogen_index,personID]             #NEW 
        for key in self.meta.pathogen_states: 
            val = (key in ['p_susceptible', 'p_naive']) # Default value is True for susceptible and naive, False otherwise 
            self[key] = np.full((pars['n_pathogens'], self.pars['pop_size']), val, dtype = bool)

        #Set pathogen variant state, which store info about which variant a person is exposed to   #NEW
        for key in self.meta.pathogen_variants:
            opt = (key in ['p_exposed_variant','p_infectious_variant','p_recovered_variant'])

            if opt: 
                self[key] = np.full((pars['n_pathogens'],self.pars['pop_size']), np.nan, dtype=cvd.default_float) #self[key][pathogen_index,personID]
            else:  
                #change this to use n_variants for pathogen i!! 
                self[key] = (np.full((pars['n_pathogens'], max_variants, self.pars['pop_size']), False, dtype=bool)) #self[key)][variant_index][personID]

        for key in self.meta.imm_levels:
            if key not in ['was_previously_infected', 'curr_min']:
                self[key] = np.full((pars['n_pathogens'],self.pars['pop_size']), 0, dtype = cvd.default_float) #imm_level
                
        self['was_previously_infected'] = np.full((pars['n_pathogens'],self.pars['pop_size']), False, dtype = bool)
         
        self['curr_min'] = np.full((pars['n_pathogens'],self.pars['pop_size']), 0.005, dtype = cvd.default_float)
        # Set variant states, which store info about which variant a person is exposed to
     #   for key in self.meta.variant_states:
      #      self[key] = np.full(self.pars['pop_size'], np.nan, dtype=cvd.default_float)
      ##  for key in self.meta.by_variant_states:
        #    self[key] = np.full((self.pars['n_variants'][0], self.pars['pop_size']), False, dtype=bool)
#
        # Set immunity and antibody states
        for key in self.meta.imm_states:  # Everyone starts out with no immunity
            self[key] = np.zeros((self.pars['n_pathogens'], max_variants, self.pars['pop_size']), dtype=cvd.default_float)

        for key in self.meta.nab_states:  # Everyone starts out with no antibodies
            dtype = cvd.default_int if key == 't_nab_event' else cvd.default_float
            self[key] = np.zeros((self.pars['n_pathogens'],self.pars['pop_size']), dtype=dtype) #size is for all pathogens, but most indices wont be used if the param use_nab_framwork is False


        for key in self.meta.vacc_states:
            self[key] = np.zeros(self.pars['pop_size'], dtype=cvd.default_int)

        # Set dates and durations -- both floats
        for key in self.meta.pathogen_dates:
            self[key] = np.full((pars['n_pathogens'], self.pars['pop_size']), np.nan, dtype=cvd.default_float) #NEW

        for key in self.meta.durs:
            self[key] = np.full((pars['n_pathogens'],self.pars['pop_size']), np.nan, dtype=cvd.default_float) #Change

        for key in self.meta.dates:
            self[key] = np.full(self.pars['pop_size'], np.nan, dtype=cvd.default_float)

        # Set dates for viral load profile -- floats
        for key in self.meta.vl_points:
            self[key] = np.full((self.pars['n_pathogens'],self.pars['pop_size']), np.nan, dtype=cvd.default_float)

        
     
        self._dtypes = {key: self[key].dtype for key in self.keys()} # Assign all to float by default 

        # Store flows
        #  to be computed during simulation
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

    def init_IgG_conversion_slope(self): 
        '''
        Select an error value for the simulation. 
        '''
        return np.random.uniform(0.14, 2.42) #values taken from the upper and lower bounds of the trendline from literature


    def init_flows(self):
        ''' Initialize flows to be zero '''
        self.flows = {}
        self.flows_variant = {}
        
        self.flows = {key:0 for key in cvd.new_result_flows}
        for p in range(len(self.pars['pathogens'])):
            self.flows[p] = {key:0 for key in cvd.new_result_flows}
            self.flows_variant[p] = {}
            for key in cvd.new_result_flows_by_variant:
                self.flows_variant[p][key] = np.zeros(self.pars['pathogens'][p].n_variants, dtype=cvd.default_float)
            if self.pars['enable_smartwatches']:
                for key in cvd.new_result_flows_smartwatches:
                    self.flows[p][key] = 0
            if self.pars['enable_multiregion']:
                targets = ['new_diagnoses', 'new_severe']
                regions = self.pars['multiregion']['rnames']
                for target in targets:
                    for region in regions:
                        self.flows[p][f'{region}_{target}'] = 0
        
        return

    def initialize(self, sim_pars=None, sim = None):
        ''' Perform initializations '''
        self.validate(sim_pars=sim_pars) # First, check that essential-to-match parameters match
        self.set_pars(sim_pars) # Replace the saved parameters with this simulation's

        for i in range(sim_pars['n_pathogens']):
            self.set_prognoses(i)
        self.sim = sim
        self.initialized = True
        return


    def set_prognoses(self, pathogen):
        '''
        Set the prognoses for each person based on age during initialization. Need
        to reset the seed because viral loads are drawn stochastically.
        '''

        pars = self.pars # Shorten 
        def find_cutoff(age_cutoffs, age):
            '''
            Find which age bin each person belongs to -- e.g. with standard
            age bins 0, 10, 20, etc., ages [5, 12, 4, 58] would be mapped to
            indices [0, 1, 0, 5]. Age bins are not guaranteed to be uniform
            width, which is why this can't be done as an array operation.
            '''
            return np.nonzero(age_cutoffs <= age)[0][-1]  # Index of the age bin to use

        cvu.set_seed(pars['rand_seed'])

        progs = pars['pathogens'][pathogen].prognoses
        inds = np.fromiter((find_cutoff(progs['age_cutoffs'], this_age) for this_age in self.age), dtype=cvd.default_int, count=len(self)) # Convert ages to indices
        self.symp_prob[pathogen,:]   = progs['symp_probs'][inds] # Probability of developing symptoms
        self.severe_prob[pathogen,:] = progs['severe_probs'][inds]*progs['comorbidities'][inds] # Severe disease probability is modified by comorbidities
        self.crit_prob[pathogen,:]   = progs['crit_probs'][inds] # Probability of developing critical disease
        self.death_prob[pathogen,:]  = progs['death_probs'][inds] # Probability of death
        self.rel_sus[pathogen,:]     = progs['sus_ORs'][inds]  # Default susceptibilities 
        self.rel_trans[pathogen,:]   = progs['trans_ORs'][inds] * cvu.sample(**self.pars['pathogens'][pathogen].beta_dist, size=len(inds))  # Default transmissibilities, with viral load drawn from a distribution


        self.abs_symp_prob[pathogen,:]   = pars['pathogens'][pathogen].absolute_prognoses['symp_probs'][inds] # Probability of developing symptoms
        self.abs_severe_prob[pathogen,:] = pars['pathogens'][pathogen].absolute_prognoses['severe_probs'][inds]*pars['pathogens'][pathogen].absolute_prognoses['comorbidities'][inds] # Severe disease probability is modified by comorbidities
        self.abs_crit_prob[pathogen,:]   = pars['pathogens'][pathogen].absolute_prognoses['crit_probs'][inds] # Probability of developing critical disease
        self.abs_death_prob[pathogen,:]  = pars['pathogens'][pathogen].absolute_prognoses['death_probs'][inds] # Probability of death


        return


    def check_alert_accuracy(self):
        '''
        Sorts alerts into true and false. Run before check_alerted() because check_alerted() will set alerted status to False in preparation for the next day.
        '''

        # Get those who were alerted and those who are infected
        watch_owners = self.true('has_watch') # list of inds of people with watches
        bool_mask_alerted = self.alerted[watch_owners]
        bool_mask_exposed = self.exposed[watch_owners]
        
        tp = np.sum(strat.get_indices_to_track(self.sim,np.logical_and(bool_mask_alerted, bool_mask_exposed)))
        fn = np.sum(strat.get_indices_to_track(self.sim,np.logical_and(np.logical_not(bool_mask_alerted), bool_mask_exposed)))
        tn = np.sum(strat.get_indices_to_track(self.sim,np.logical_and(np.logical_not(bool_mask_alerted), np.logical_not(bool_mask_exposed))))
        fp = np.sum(strat.get_indices_to_track(self.sim,np.logical_and(bool_mask_alerted, np.logical_not(bool_mask_exposed))))

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
          
        self.is_exp = []
        
        self.init_flows()

        new_total_deaths = 0
        new_total_known_deaths = 0
        mexp = self.true('exposed')

        for pathogen in range(len(self.pars['pathogens'])):
            self.is_exp.append(self.true_with_index('p_exposed', pathogen)) # For storing the interim values since used in every subsequent calculation 
            # Perform updates
            self.flows[pathogen]['new_infectious']    += self.check_infectious(pathogen = pathogen) # For people who are exposed and not infectious, check if they begin being infectious
            self.flows[pathogen]['new_symptomatic']   += self.check_symptomatic(pathogen = pathogen)
             
            self.flows[pathogen]['new_severe']        += self.check_severe(pathogen = pathogen)

            self.flows[pathogen]['new_critical']      += self.check_critical(pathogen = pathogen)
            self.flows[pathogen]['new_recoveries']    += self.check_recovery(pathogen = pathogen)
            new_deaths, new_known_deaths     = self.check_death(pathogen = pathogen)
            new_total_deaths += new_deaths
            new_total_known_deaths += new_known_deaths
            self.flows[pathogen]['new_deaths']        += new_deaths
            self.flows[pathogen]['new_known_deaths']  += new_known_deaths
            
             
        #update flows for overall people state
         
        self.flows['new_infectious']    += len(strat.get_indices_to_track(self.sim, self.check_inds(self.infectious, self.date_infectious, filter_inds=mexp)))
        self.flows['new_symptomatic']   +=len(strat.get_indices_to_track(self.sim, self.check_inds(self.symptomatic, self.date_symptomatic, filter_inds=mexp)))
    
        self.flows['new_severe']        += len(strat.get_indices_to_track(self.sim,self.check_inds(self.severe, self.date_severe, filter_inds=mexp)))

        self.flows['new_critical']      += len(strat.get_indices_to_track(self.sim,self.check_inds(self.critical, self.date_critical, filter_inds=mexp)))
       
  
        self.flows['new_recoveries']    += len(strat.get_indices_to_track(self.sim,self.check_inds(self.recovered, self.date_recovered, filter_inds=mexp)))
        self.flows['new_deaths']        +=  new_total_deaths
        self.flows['new_known_deaths']  += new_total_known_deaths
        
        self.merge_states(self.pars['n_pathogens'], False, False) 
        return


    def update_states_post(self, pathogen = 0):
        ''' Perform post-timestep updates '''
         
        self.flows[pathogen]['new_diagnoses'] += len(strat.get_indices_to_track(self.sim, self.check_inds(self.diagnosed[pathogen], self.date_pos_test)))
        quar = self.check_quar()
        self.flows[pathogen]['new_quarantined'] += quar
        
        self.flows['new_diagnoses'] += self.check_diagnosed(pathogen = pathogen)

        if pathogen == 0:
            self.flows['new_quarantined'] += quar


        # Take care of smartwatch calculations
        if self.pars['enable_smartwatches']:
            tp, fn, tn, fp = self.check_alert_accuracy()
            self.flows[pathogen]['new_alerts_tp'] += tp
            self.flows[pathogen]['new_alerts_fn'] += fn
            self.flows[pathogen]['new_alerts_tn'] += tn
            self.flows[pathogen]['new_alerts_fp'] += fp
            self.flows[pathogen]['new_alerted'] += self.check_alerted()
            sw_quar, sw_i_quar = self.check_sw_quarantined() # Update the number of smartwatch users quarantined. 
            self.flows[pathogen]['new_Q_w'] += sw_quar
            self.flows[pathogen]['new_Q_w_i'] += sw_i_quar

        if pathogen == 0:
            del self.is_exp  # Tidy up  

        return


    def check_sw_quarantined(self):
        sw_quarantined_arr = np.logical_and(self.quarantined, self.has_watch)
        sw_quarantined_count = np.sum(strat.get_indices_to_track(self.sim, sw_quarantined_arr))
        sw_correctly_quarantined_arr = np.logical_and(sw_quarantined_arr, self.exposed)
        sw_i_quarantined_count = sw_quarantined_count - np.sum(strat.get_indices_to_track(self.sim, sw_correctly_quarantined_arr))

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
     
    #region
    def check_infectious(self, pathogen = 0):
        ''' Check if they become infectious with the pathogen index passed'''
        inds = self.check_inds(self.p_infectious[pathogen], self.date_p_infectious[pathogen], filter_inds=self.is_exp[pathogen])
        self.p_infectious[pathogen, inds] = True 
        self.p_infectious_variant[pathogen, inds] = self.p_exposed_variant[pathogen,inds]
        
        

        for variant in range(self.pars['pathogens'][pathogen].n_variants):
            this_variant_inds = cvu.itrue(self.p_infectious_variant[pathogen, inds] == variant, inds)
            n_this_variant_inds = len(strat.get_indices_to_track(self.sim, this_variant_inds))
            self.flows_variant[pathogen]['new_infectious_by_variant'][variant] += n_this_variant_inds

            #TODO change that array and update this 
            self.p_infectious_by_variant[pathogen,variant, this_variant_inds] = True 
        return len(strat.get_indices_to_track(self.sim, inds))


    def check_symptomatic(self, pathogen = 0):
        ''' Check for new progressions to symptomatic '''
        inds = self.check_inds(self.p_symptomatic[pathogen], self.date_p_symptomatic[pathogen], filter_inds=self.is_exp[pathogen])
        
        self.p_symptomatic[pathogen, inds] = True
        return len(strat.get_indices_to_track(self.sim, inds))

 
    def check_severe(self, pathogen = 0):
        ''' Check for new progressions to severe '''
        inds = self.check_inds(self.p_severe[pathogen], self.date_p_severe[pathogen], filter_inds=self.is_exp[pathogen])

        strat_inds = strat.get_indices_to_track(self.sim, inds)

        if self.pars['enable_multiregion']:
            rnames = self.pars['multiregion']['rnames']
            rsizes = self.pars['multiregion']['rsizes']
            rstarts = self.pars['multiregion']['rstarts']

            for rname, rstart, rsize in zip(rnames, rstarts, rsizes):
                # Count the # indicies in the region. i.e. indicies greater than rstart and less than the end .
                num_inds = np.sum((strat_inds >= rstart) & (strat_inds < rstart+rsize))
                self.flows[pathogen][f'{rname}_new_severe'] += num_inds
            

        self.p_severe[pathogen, inds] = True
        return len(strat_inds)


    def check_critical(self,pathogen = 0):
        ''' Check for new progressions to critical '''
        inds = self.check_inds(self.p_critical[pathogen], self.date_p_critical[pathogen], filter_inds=self.is_exp[pathogen])
        self.p_critical[pathogen,inds] = True
        return len(strat.get_indices_to_track(self.sim, inds))


    def check_recovery(self, pathogen = 0, inds=None, filter_inds='is_exp'):
        '''
        Check for recovery.

        More complex than other functions to allow for recovery to be manually imposed
        for a specified set of indices.
        '''

        # Handle more flexible options for setting indices
        if filter_inds == 'is_exp':
            filter_inds = self.is_exp[pathogen]
        if inds is None:
            inds = self.check_inds(self.p_recovered[pathogen], self.date_p_recovered[pathogen], filter_inds=filter_inds)
             
        self.is_coinfected[inds] = False
        self.p_exposed[pathogen,inds]          = False
        self.p_infectious[pathogen,inds]       = False
        self.p_symptomatic[pathogen,inds]      = False
        self.p_severe[pathogen,inds]           = False
        self.p_critical[pathogen,inds]         = False
        self.p_recovered[pathogen,inds]        = True 
           
        self.p_recovered_variant[pathogen,inds] = self.p_exposed_variant[pathogen,inds] 
        self.p_infectious_variant[pathogen,inds] = np.nan
        self.p_exposed_variant[pathogen,inds]    = np.nan
         
        self.p_exposed_by_variant[pathogen, :, inds] = False
        self.p_infectious_by_variant[pathogen, :, inds] = False

        
        #self.date_p_exposed[pathogen] = np.nan
        #self.date_p_infectious[pathogen] = np.nan
        #self.date_p_symptomatic[pathogen] = np.nan
        ##self.date_p_severe[pathogen] = np.nan
        #self.date_p_critical[pathogen] = np.nan 
        # Handle immunity aspects
        if self.pars['use_waning']:

            # Reset additional states 
            self.diagnosed[inds]   = False # Reset their diagnosis state because they might be reinfected
            self.p_susceptible[pathogen,inds] = True
            self.p_diagnosed[pathogen,inds]   = False # Reset their diagnosis state because they might be reinfected

        # # Handle instances of false positive diagnoses; if someone is not exposed, is diagnosed, and has passed the required time in isolation, un-diagnose them
        # fp_diagnoses = cvu.true((~self.exposed) & (self.diagnosed) & (self.t - self.date_diagnosed >= 14))
        # self.diagnosed[fp_diagnoses] = False
        
        return len(strat.get_indices_to_track(self.sim, inds))


    def check_death(self, pathogen= 0):
        ''' Check whether or not this person died on this timestep  '''
        inds = self.check_inds(self.p_dead[pathogen], self.date_p_dead[pathogen], filter_inds=self.is_exp[pathogen])
    
       
        self.dead[inds]             = True
         
        diag_inds = inds[self.p_diagnosed[pathogen, inds]] # Check whether the person was diagnosed before dying
        self.known_dead[diag_inds]  = True
        self.dead[inds]  = True
        self.susceptible[inds]      = False
        self.exposed[inds]          = False
        self.infectious[inds]       = False
        self.symptomatic[inds]      = False
        self.severe[inds]           = False
        self.critical[inds]         = False
        self.known_contact[inds]    = False
        self.quarantined[inds]      = False
        self.recovered[inds]        = False 
          
        self.p_dead[pathogen,inds]             = True 

        for i in range(len(self.pars['pathogens'])):
            self.p_susceptible[i,inds]      = False
            self.p_exposed[i,inds]          = False
            self.p_infectious[i,inds]       = False
            self.p_symptomatic[i,inds]      = False
            self.p_severe[i,inds]           = False
            self.p_critical[i,inds]         = False  
            self.p_recovered[i,inds]        = False
            self.p_infectious_variant[i,inds] = np.nan 
            self.p_exposed_variant[i,inds]    = np.nan
            self.p_recovered_variant[i,inds]  = np.nan 
            
        return len(strat.get_indices_to_track(self.sim, inds)), len(strat.get_indices_to_track(self.sim, diag_inds))


    def check_alerted(self):
        '''
        Check which people received alerts on this timestep and whether they were correct or incorrect.
        Reset alerted status to False for everyone to prepare for the next day.
        '''

        # Check who was alerted
        n_alerted =  len(strat.get_indices_to_track(self.sim, cvu.true(self.alerted)))
        #TODO potentially save for which pathogen the person is alerted to ?

        # Set the alerted status of everyone to false
        self.alerted[:] = False

        return n_alerted

     
    def check_diagnosed(self, pathogen = 0):
        '''
        Check for new diagnoses. Since most data are reported with diagnoses on
        the date of the test, this function reports counts not for the number of
        people who received a positive test result on a day, but rather, the number
        of people who were tested on that day who are schedule to be diagnosed in
        the future.
        '''

        # Handle people who tested today who will be diagnosed in future (i.e., configure them to have have finished being tested)
        test_pos_inds = self.check_inds(self.p_diagnosed[pathogen], self.date_pos_test, filter_inds=None) # Find people who are not diagnosed and have a date of a positive test that is today or earlier
        self.date_pos_test[test_pos_inds] = np.nan # Clear date of having will-be-positive test
        
        test_pos_inds_strat = strat.get_indices_to_track(self.sim, test_pos_inds)
        # Update per variant. (Written 31-05-23)
        for variant in range(self.pars['pathogens'][pathogen].n_variants):
            # Get the people who are exposed for this variant, and have a positive test. 
            # TODO: exposed_variant ok? 
            var_count = np.sum(self.p_exposed_by_variant[pathogen, variant, test_pos_inds_strat])
            self.flows_variant[pathogen]['new_diagnoses_by_variant'][variant] += var_count


        
        if self.pars['enable_multiregion']:
            #### REGIONAL RESULTS
            rnames = self.pars['multiregion']['rnames']
            rsizes = self.pars['multiregion']['rsizes']
            rstarts = self.pars['multiregion']['rstarts']

            for rname, rstart, rsize in zip(rnames, rstarts, rsizes): 
                # Count the # indicies in the region. i.e. indicies greater than rstart and less than the end .
                num_inds = np.sum((test_pos_inds_strat >= rstart) & (test_pos_inds_strat < rstart+rsize))
                self.flows[pathogen][f'{rname}_new_diagnoses'] += num_inds
                
       
        # Handle people who were actually diagnosed today (i.e., set them as diagnosed and remove any of them that were quarantining from quarantine)
        diag_inds  = self.check_inds(self.diagnosed, self.date_diagnosed, filter_inds=None) # Find who are not diagnosed and have a date of diagnosis that is today or earlier
        self.diagnosed[diag_inds]   = True # Set these people to be diagnosed
        self.p_diagnosed[pathogen, diag_inds]   = True # Set these people to be diagnosed
        quarantined = cvu.itruei(self.quarantined, diag_inds) # Find individuals who were just diagnosed who are in quarantine
        self.date_end_quarantine[quarantined] = self.t # Set end quarantine date to match when the person left quarantine (and entered isolation)
        self.quarantined[diag_inds] = False # If you are diagnosed, you are isolated, not in quarantine

        return len(test_pos_inds_strat)


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

                if self.sim.enable_stratifications:
                    if ind in self.sim.stratification_indices:
                        n_quarantined += 1
                else: 
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

    #endregion
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

        for i in range(self.pars['n_pathogens']):
            pathogen_index = i

            
            for key in self.meta.pathogen_states:
                if key in ['p_susceptible', 'p_naive']:
                      self[key][pathogen_index, inds] = True
                else:
                    if (key != 'vaccinated') or reset_vx: # Don't necessarily reset vaccination
                        self[key][pathogen_index,inds] = False 

            for key in self.meta.pathogen_variants:
                if key in ['p_exposed_variant','p_infectious_variant','p_recovered_variant']:
                    self[key][pathogen_index,inds] = np.nan
                else:
                    self[key][pathogen_index,:, inds] = False

                 
            # Reset immunity and antibody states TODO change to multipathogen 
            non_vx_inds = inds if reset_vx else inds[~self['vaccinated'][inds]]
            for key in self.meta.imm_states:
                self[key][pathogen_index,:, non_vx_inds] = 0

            for key in self.meta.nab_states:
                self[key][pathogen_index, non_vx_inds] = 0

            # Reset dates
            for key in self.meta.dates:
                if (key != 'date_vaccinated') or reset_vx: # Don't necessarily reset vaccination
                    self[key][inds] = np.nan
        
            for key in self.meta.durs:
                if (key != 'date_vaccinated') or reset_vx: # Don't necessarily reset vaccination
                    self[key][pathogen_index,inds] = np.nan
                
            for key in self.meta.pathogen_dates:
                self[key][pathogen_index,inds] = np.nan
            
        for key in self.meta.vacc_states: 
            self[key][non_vx_inds] = 0
        return


    def make_nonnaive(self, inds, pathogen_index = 0,set_recovered=False, date_recovered=0):
        '''
        Make a set of people non-naive.

        This can be done either by setting only susceptible and naive states,
        or else by setting them as if they have been infected and recovered.
        '''


        self.make_naive(inds,pathogen_index) # First make them naive and reset all other states

        for i in range(self.pars['n_pathogen']):
            pathogen_index = i
            # Make them non-naive
            for key in ['susceptible', 'naive']:
                self[key][inds] = False
                self[f'p_{key}'][pathogen_index,inds] = False

            if set_recovered:
                self.date_recovered[inds] = date_recovered # Reset date recovered
                self.date_p_recovered[pathogen_index,inds] = date_recovered # Reset date recovered
                self.check_recovery(inds=inds, filter_inds=None, pathogen = pathogen_index) # Set recovered

        return



    def infect(self, inds, hosp_max=None, icu_max=None, source=None, layer=None, variant=0, pathogen_index = 0):
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
         
        inds = np.intersect1d(inds, np.logical_not(self.is_coinfected).nonzero()[0]) #keep inds susceptible and not co-infected 
        # Remove duplicates
        inds, unique = np.unique(inds, return_index=True)


        if source is not None:
            source = source[unique]
            
        # Keep only susceptibles 
        keep = self.p_susceptible[pathogen_index,inds] # Unique indices in inds and source that are also susceptible  
        inds = inds[keep]
        
        if source is not None:
            source = source[keep]
        
             
        #severity_mult_by_ind and dur_mult_by_ind will contain the multiplier of severity and duration
        #for everyone which depends on if they are co-infected with a pathogen or not
        Msev_by_ind = np.ones(self.pars['pop_size'])  
        co_infected_pathogen_inds = np.full(self.pars['pop_size'], -1, dtype = int) #Indices of pathogens co-infecting individuals, -1 if person is not infected with any pathogen 
        Mdur_by_ind = np.ones(self.pars['pop_size']) 

        for i in range(self.pars['n_pathogens']):
            if i == pathogen_index:
                assert len(np.intersect1d(inds,self.p_exposed[i].nonzero()[0])) == 0
                continue
            indices_infected = self.p_exposed[i].nonzero()[0]  #indices with pathogen index i
            indices_infected = np.intersect1d(indices_infected, inds)
              
            co_infected_pathogen_inds[indices_infected] = i #They are co-infected with pathogens[i]

            Msev_by_ind[indices_infected] = self.sim['Msev'][i, pathogen_index] #Set the indices with people infected with i to have the Msev multiplier at i
            Mdur_by_ind[indices_infected] = self.sim['Mdur'][pathogen_index, i] #Set the indices with people infected with i to have the Mdur multiplier at i
              
        indices_use_alpha = np.where(co_infected_pathogen_inds != -1)[0] #Indices of people that are infected
        
        self.sim.results['co-infections'][self.t] = len(indices_use_alpha) 
        self.is_coinfected[indices_use_alpha] = True

        '''
        for i in range(self.pars['n_pathogens']): 
            if i != pathogen_index: 
                for j in range(self.sim.pathogens[i].n_variants):
                    indices_with_variant = self.p_exposed_by_variant[i,j].nonzero()[0]  
                    indices_with_variant = np.intersect1d(indices_with_variant, inds)
                    for ind in indices_with_variant:
                        self.recalculate_disease_trajectory(ind, i, j, Msev_by_ind, Mdur_by_ind, pathogen_index)'''

 
        # Deal with variant parameters
        variant_keys = ['rel_symp_prob', 'rel_severe_prob', 'rel_crit_prob', 'rel_death_prob']
        infect_pars = {
           'rel_symp_prob':self.pars['pathogens'][pathogen_index].rel_symp_prob,
           'rel_severe_prob':self.pars['pathogens'][pathogen_index].rel_severe_prob,
           'rel_crit_prob':self.pars['pathogens'][pathogen_index].rel_crit_prob,
           'rel_death_prob':self.pars['pathogens'][pathogen_index].rel_death_prob}
        
        variant_label =self.pars['pathogens'][pathogen_index].get_variants_labels()[variant]
        if variant!=0:
            for k in variant_keys:
                infect_pars[k] *= self.pars['pathogens'][pathogen_index].variants[variant-1].p[k] 

        n_infections = len(inds)
        durpars      = self.pars['pathogens'][pathogen_index].dur
        
        # Retrieve those with a breakthrough infection (defined nabs)
        if self.pars['pathogens'][pathogen_index].use_nab_framework:
            breakthrough_inds = inds[cvu.true(self.peak_nab[pathogen_index,inds])]
            if len(breakthrough_inds):
                no_prior_breakthrough = (self.n_breakthroughs[pathogen_index, breakthrough_inds] == 0) # We only adjust transmissibility for the first breakthrough
                new_breakthrough_inds = breakthrough_inds[no_prior_breakthrough]
                self.rel_trans[pathogen_index, new_breakthrough_inds] *= self.pars['pathogens'][pathogen_index].trans_redux
        else:
            breakthrough_inds = []

        # Update states, variant info, and flows 
        self.p_susceptible[pathogen_index, inds]    = False
        self.p_naive[pathogen_index, inds]          = False
        self.p_recovered[pathogen_index, inds]      = False
        self.p_diagnosed[pathogen_index, inds]      = False
        self.p_exposed[pathogen_index, inds]        = True

        self.n_infections[pathogen_index, inds]  += 1
        self.n_breakthroughs[pathogen_index, breakthrough_inds] += 1
        
        self.p_exposed_variant[pathogen_index, inds] = variant
        self.p_exposed_by_variant[pathogen_index,variant, inds] = True
          
        self.flows[pathogen_index]['new_infections']   += len(strat.get_indices_to_track(self.sim, inds))
        self.flows[pathogen_index]['new_reinfections'] += len(strat.get_indices_to_track(self.sim, cvu.defined(self.date_p_recovered[pathogen_index, inds]))) # Record reinfections
         
        self.flows['new_infections']   += len(strat.get_indices_to_track(self.sim, inds))
        self.flows['new_reinfections'] += len(strat.get_indices_to_track(self.sim, cvu.defined(self.date_p_recovered[pathogen_index, inds]))) # Record reinfections

        self.flows_variant[pathogen_index]['new_infections_by_variant'][variant] += len(strat.get_indices_to_track(self.sim, inds))
        
        # Record transmissions
        for i, target in enumerate(inds):
            entry = dict(source=source[i] if source is not None else None, target=target, date=self.t, layer=layer, variant=variant_label)
            self.infection_log.append(entry)
            
        # Calculate how long before this person can infect other people
        self.dur_exp2inf[pathogen_index, inds] = cvu.sample(**durpars['exp2inf'], size=n_infections) * Mdur_by_ind[inds]
        
        self.date_p_exposed[pathogen_index, inds] = self.t
        self.date_p_infectious[pathogen_index, inds] = self.dur_exp2inf[pathogen_index, inds] + self.t
         
        #set peak immunity  date
        if not self.pars['pathogens'][pathogen_index].use_nab_framework:
            self['t_peak_imm'][pathogen_index, inds] = self.t + self.pars['pathogens'][pathogen_index].imm_days_to_peak
            inds_zeros = np.intersect1d(np.where(self['imm_level'][pathogen_index] <= 0.001), inds)
            self['imm_level'][pathogen_index, inds_zeros] = 0.005 

        # Reset all other dates  
       
        for key in ['date_p_symptomatic', 'date_p_severe', 'date_p_critical', 'date_p_diagnosed', 'date_p_recovered']:
            self[key][pathogen_index, inds] = np.nan


        # Use prognosis probabilities to determine what happens to them
        symp_prop_coinf_mult = np.ones(self.pars['pop_size'])
         
        #Scale symp_probs based on co-infection
        for i in indices_use_alpha:
            assert co_infected_pathogen_inds[i] != -1 
            symp_prop_coinf_mult[i] = pint.get_disease_traj_alpha(self.abs_symp_prob[co_infected_pathogen_inds[i],i], self.abs_symp_prob[pathogen_index,i], Msev_by_ind[i])
             

        symp_probs = infect_pars['rel_symp_prob']*self.symp_prob[pathogen_index,inds]*(1-self.symp_imm[pathogen_index,variant, inds])*symp_prop_coinf_mult[inds] # Calculate their actual probability of being symptomatic
        
        symp_probs = np.clip(symp_probs, 0, 1)

        is_symp = cvu.binomial_arr(symp_probs) # Determine if they develop symptoms
           
        symp_inds = inds[is_symp]
        asymp_inds = inds[~is_symp] # Asymptomatic
        self.flows_variant[pathogen_index]['new_symptomatic_by_variant'][variant] += len(strat.get_indices_to_track(self.sim, symp_inds))

        # CASE 1: Asymptomatic: may infect others, but have no symptoms and do not die
        dur_asym2rec = cvu.sample(**durpars['asym2rec'], size=len(asymp_inds)) * Mdur_by_ind[asymp_inds]
        
        self.date_p_recovered[pathogen_index, asymp_inds] = self.date_p_infectious[pathogen_index, asymp_inds] + dur_asym2rec  # Date they recover. PxN.
        self.date_p_symptomatic[pathogen_index, asymp_inds] = np.nan
        self.dur_disease[pathogen_index, asymp_inds] = self.dur_exp2inf[pathogen_index, asymp_inds] + dur_asym2rec  # Store how long this person had COVID-19. Nx1.  
        
        # CASE 2: Symptomatic: can either be mild, severe, or critical
        n_symp_inds = len(symp_inds)
        self.dur_inf2sym[pathogen_index, symp_inds] = cvu.sample(**durpars['inf2sym'], size=n_symp_inds)* Mdur_by_ind[symp_inds]# Store how long this person took to develop symptoms  
        self.date_p_symptomatic[pathogen_index, symp_inds] = self.date_p_infectious[pathogen_index, symp_inds] + self.dur_inf2sym[pathogen_index, symp_inds] # Date they become symptomatic
         
        sev_prop_coinf_mult = np.ones(self.pars['pop_size'])
         
        #Scale symp_probs based on co-infection
        for i in np.intersect1d(indices_use_alpha, symp_inds):
            assert co_infected_pathogen_inds[i] != -1 
            sev_prop_coinf_mult[i] = pint.get_disease_traj_alpha(self.abs_severe_prob[co_infected_pathogen_inds[i],i], self.abs_severe_prob[pathogen_index,i], Msev_by_ind[i]) 
             

        sev_probs = infect_pars['rel_severe_prob'] * self.severe_prob[pathogen_index, symp_inds]*(1-self.sev_imm[pathogen_index, variant, symp_inds]) * sev_prop_coinf_mult[symp_inds] # Probability of these people being severe
        
        sev_probs = np.clip(sev_probs, 0, 1)
        is_sev = cvu.binomial_arr(sev_probs) # See if they're a severe or mild case
        sev_inds = symp_inds[is_sev]
        mild_inds = symp_inds[~is_sev] # Not severe
        self.flows_variant[pathogen_index]['new_severe_by_variant'][variant] += len(strat.get_indices_to_track(self.sim, sev_inds))
        self.date_p_severe[pathogen_index,mild_inds] = np.nan
        # CASE 2.1: Mild symptoms, no hospitalization required and no probability of death
        dur_mild2rec = cvu.sample(**durpars['mild2rec'], size=len(mild_inds))* Mdur_by_ind[mild_inds]
        self.date_p_recovered[pathogen_index, mild_inds] = self.date_p_symptomatic[pathogen_index,mild_inds] + dur_mild2rec  # Date they recover
        self.dur_disease[pathogen_index, mild_inds] = self.dur_exp2inf[pathogen_index, mild_inds] + self.dur_inf2sym[pathogen_index,mild_inds] + dur_mild2rec  # Store how long this person had COVID-19

        # CASE 2.2: Severe cases: hospitalization required, may become critical

        crit_prop_coinf_mult = np.ones(self.pars['pop_size'])
         
        #Scale symp_probs based on co-infection
        for i in np.intersect1d(indices_use_alpha, sev_inds):
            assert co_infected_pathogen_inds[i] != -1 
            crit_prop_coinf_mult[i] = pint.get_disease_traj_alpha(self.abs_crit_prob[co_infected_pathogen_inds[i],i], self.abs_crit_prob[pathogen_index,i], Msev_by_ind[i])
            
        self.dur_sym2sev[pathogen_index,sev_inds] = cvu.sample(**durpars['sym2sev'], size=len(sev_inds)) * Mdur_by_ind[sev_inds]# Store how long this person took to develop severe symptoms
        self.date_p_severe[pathogen_index,sev_inds] = self.date_p_symptomatic[pathogen_index, sev_inds] + self.dur_sym2sev[pathogen_index, sev_inds]  # Date symptoms become severe
        crit_probs = infect_pars['rel_crit_prob'] * self.crit_prob[pathogen_index, sev_inds] * (self.pars['no_hosp_factor'] if hosp_max else 1.) * crit_prop_coinf_mult[sev_inds] # Probability of these people becoming critical - higher if no beds available
        crit_probs = np.clip(crit_probs, 0, 1)
        is_crit = cvu.binomial_arr(crit_probs)  # See if they're a critical case
        crit_inds = sev_inds[is_crit]
        non_crit_inds = sev_inds[~is_crit]
        self.date_p_critical[pathogen_index, non_crit_inds] = np.nan
        # CASE 2.2.1 Not critical - they will recover
        dur_sev2rec = cvu.sample(**durpars['sev2rec'], size=len(non_crit_inds)) * Mdur_by_ind[non_crit_inds]
        self.date_p_recovered[pathogen_index, non_crit_inds] = self.date_p_severe[pathogen_index, non_crit_inds] + dur_sev2rec  # Date they recover
        self.dur_disease[pathogen_index, non_crit_inds] = self.dur_exp2inf[pathogen_index,non_crit_inds] + self.dur_inf2sym[pathogen_index,non_crit_inds] + self.dur_sym2sev[pathogen_index,non_crit_inds] + dur_sev2rec  # Store how long this person had COVID-19

        # CASE 2.2.2: Critical cases: ICU required, may die

        death_prop_coinf_mult = np.ones(self.pars['pop_size'])
         
        #Scale symp_probs based on co-infection
        for i in np.intersect1d(indices_use_alpha, crit_inds):
            assert co_infected_pathogen_inds[i] != -1 
            death_prop_coinf_mult[i] = pint.get_disease_traj_alpha(self.abs_death_prob[co_infected_pathogen_inds[i],i], self.abs_death_prob[pathogen_index,i], Msev_by_ind[i])
            

        self.dur_sev2crit[pathogen_index, crit_inds] = cvu.sample(**durpars['sev2crit'], size=len(crit_inds)) * Mdur_by_ind[crit_inds]
        self.date_p_critical[pathogen_index, crit_inds] = self.date_p_severe[pathogen_index, crit_inds] + self.dur_sev2crit[pathogen_index, crit_inds]  # Date they become critical
        death_probs = infect_pars['rel_death_prob'] * self.death_prob[pathogen_index, crit_inds] * (self.pars['no_icu_factor'] if icu_max else 1.) * death_prop_coinf_mult[crit_inds]# Probability they'll die
        death_probs = np.clip(death_probs, 0, 1)
        is_dead = cvu.binomial_arr(death_probs)  # Death outcome
        dead_inds = crit_inds[is_dead]
        alive_inds = crit_inds[~is_dead]
         
        self.date_p_dead[pathogen_index, alive_inds] = np.nan

        # CASE 2.2.2.1: Did not die
        dur_crit2rec = cvu.sample(**durpars['crit2rec'], size=len(alive_inds)) * Mdur_by_ind[alive_inds]
        self.date_p_recovered[pathogen_index, alive_inds] = self.date_p_critical[pathogen_index,alive_inds] + dur_crit2rec # Date they recover
        self.dur_disease[pathogen_index, alive_inds] = self.dur_exp2inf[pathogen_index, alive_inds] + self.dur_inf2sym[pathogen_index,alive_inds] + self.dur_sym2sev[pathogen_index, alive_inds] + self.dur_sev2crit[pathogen_index, alive_inds] + dur_crit2rec  # Store how long this person had COVID-19

        # CASE 2.2.2.2: Did die
        dur_crit2die = cvu.sample(**durpars['crit2die'], size=len(dead_inds)) * Mdur_by_ind[dead_inds]
        self.date_p_dead[pathogen_index, dead_inds] = self.date_p_critical[pathogen_index, dead_inds] + dur_crit2die # Date of death 
        self.dur_disease[pathogen_index, dead_inds] = self.dur_exp2inf[pathogen_index,dead_inds] + self.dur_inf2sym[pathogen_index,dead_inds] + self.dur_sym2sev[pathogen_index,dead_inds] + self.dur_sev2crit[pathogen_index,dead_inds] + dur_crit2die   # Store how long this person had COVID-19
       
        self.date_p_recovered[pathogen_index,dead_inds] = np.nan # If they did die, remove them from recovered

        

        if not self.pars['pathogens'][pathogen_index].use_nab_framework:
            rec_inds = list(set(inds) - set(dead_inds))
            self['t_min_imm'][pathogen_index, rec_inds] = self.date_p_recovered[pathogen_index, rec_inds] + self.pars['pathogens'][pathogen_index].imm_days_to_final_value
            self['t_min_imm'][pathogen_index, dead_inds] = self.date_p_dead[pathogen_index, dead_inds] + self.pars['pathogens'][pathogen_index].imm_days_to_final_value

        # HANDLE VIRAL LOAD CONTROL POINTS
       
        # Get P_inf: where viral load crosses 10^6 cp/mL
        self.x_p_inf[pathogen_index,inds] = self.dur_exp2inf[pathogen_index, inds]
        self.y_p_inf[pathogen_index,inds] = 6

        # Get P1: where viral load crosses 10^3 cp/mL; time difference obtained empirically through simulation
        self.x_p1[pathogen_index,inds] = np.maximum(self.x_p_inf[pathogen_index,inds] - (np.random.gamma(2, 0.35, size=len(inds)) + 0.25), 0)
        self.y_p1[pathogen_index,inds] = 3

        # Get P2: where viral load peaks; time difference obtained empirically through simulation
        self.x_p2[pathogen_index,inds] = self.x_p_inf[pathogen_index,inds] + (np.random.gamma(3, 0.26, size=len(inds)) + 0.1)
        self.y_p2[pathogen_index,inds] = ((self.y_p_inf[pathogen_index,inds] - self.y_p1[pathogen_index,inds])*(self.x_p2[pathogen_index,inds] - self.x_p1[pathogen_index,inds])/(self.x_p_inf[pathogen_index,inds] - self.x_p1[pathogen_index,inds])) + self.y_p1[pathogen_index,inds]

        # Align P1, P_inf, and P2 to current time
        self.x_p1[pathogen_index,inds] = self.x_p1[pathogen_index,inds] + self.t
        self.x_p_inf[pathogen_index,inds] = self.x_p_inf[pathogen_index,inds] + self.t
        self.x_p2[pathogen_index,inds] = self.x_p2[pathogen_index,inds] + self.t

        # Get P3: where viral load drops below 10^6 cp/mL
        time_recovered = np.ones(len(self.date_p_recovered[pathogen_index]), dtype=cvd.default_float)*self.date_p_recovered[pathogen_index] # This is needed to make a copy
        inds_dead = ~np.isnan(self.date_p_dead[pathogen_index])
        time_recovered[inds_dead] = self.date_p_dead[pathogen_index,inds_dead]
        self.x_p3[pathogen_index,inds] = np.maximum(time_recovered[inds], self.x_p2[pathogen_index,inds])
        self.y_p3[pathogen_index,inds] = 6
          
        #self.sim.results['co-infected_deaths'][self.t] += len(self.is_coinfected[dead_inds].nonzero()[0])
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
        if self.pars['pathogens'][pathogen_index].use_nab_framework: 
            symp = dict(asymp=asymp_inds, mild=mild_inds, sev=sev_inds)
            cvi.update_peak_nab(self, inds, nab_pars=self.pars['pathogens'][pathogen_index], symp=symp, pathogen = pathogen_index)
            
        self.merge_states(self.pars['n_pathogens'], False, True)

        return n_infections # For incrementing counters
     
    def recalculate_disease_trajectory(self, index, pathogen_index, variant, Msev_by_ind, Mdur_by_ind, pathogen_coinfecting):
    
        #print(self.t - self.date_p_exposed[pathogen_index, index])
        #self.print_disease_trajectory(index, pathogen_index)
        #identify which states are in the past
        #Infectious Symptomatic Severe Critical Dead
        # Deal with variant parameters
        variant_keys = ['rel_symp_prob', 'rel_severe_prob', 'rel_crit_prob', 'rel_death_prob']
        infect_pars = {
           'rel_symp_prob':self.pars['pathogens'][pathogen_index].rel_symp_prob,
           'rel_severe_prob':self.pars['pathogens'][pathogen_index].rel_severe_prob,
           'rel_crit_prob':self.pars['pathogens'][pathogen_index].rel_crit_prob,
           'rel_death_prob':self.pars['pathogens'][pathogen_index].rel_death_prob}
         
        if variant!=0:
            for k in variant_keys:
                infect_pars[k] *= self.pars['pathogens'][pathogen_index].variants[variant-1].p[k] 


        alpha = (self.t - self.date_p_exposed[pathogen_index, index])/(self.date_p_recovered[pathogen_index, index] - self.date_p_exposed[pathogen_index, index])
        lbd = (Mdur_by_ind[index] - alpha)/(1-alpha)
         
        durpars      = self.pars['pathogens'][pathogen_index].dur
           
        # Calculate how long before this person can infect other people

        #Update duration spent in each state
        if self.date_p_exposed[pathogen_index, index] + self.dur_exp2inf[pathogen_index, index] > self.t: #where in exposed state still
            
            if self.t + lbd*(self.date_p_infectious[pathogen_index, index] - self.t) <= self.date_p_infectious[pathogen_index, index]: 
                self.dur_exp2inf[pathogen_index, index] = self.t - self.date_p_exposed[pathogen_index, index]
                self.date_p_exposed[pathogen_index,index] = self.t+1
            else:
                self.dur_exp2inf[pathogen_index, index] = (self.t - self.date_p_exposed[pathogen_index, index]) + lbd*(self.date_p_infectious[pathogen_index, index]-self.t)
                self.date_p_infectious[pathogen_index, index] = self.date_p_exposed[pathogen_index, index] +  self.dur_exp2inf[pathogen_index, index]
         
         
        #Scale symp_probs based on co-infection 

        is_symp = not np.isnan(self.date_p_symptomatic[pathogen_index, index])
        if (not is_symp) or (self.date_p_symptomatic[pathogen_index, index] > self.t): #If state was orinally not supposed to be entered or the state's date is in the future
            
            if (Msev_by_ind[index] >= 1 and not is_symp) or (Msev_by_ind[index] <= 1 and (not np.isnan(self.date_p_symptomatic[pathogen_index, index]))): 
                symp_prob_mult = pint.get_disease_traj_alpha(self.abs_symp_prob[pathogen_index,index], self.abs_symp_prob[pathogen_coinfecting,index], Msev_by_ind[index]) 
              
                symp_probs = infect_pars['rel_symp_prob']*self.symp_prob[pathogen_index,index]*(0.9)*symp_prob_mult #placeholder 0.9 to assume very low symp_imm protection originally
                symp_probs = symp_probs if symp_probs <=1 else 1

                prev = is_symp
                is_symp = len(cvu.binomial_arr([symp_probs])!= 0) # Determine if they develop symptoms
                if is_symp and not prev:
                    self.flows_variant[pathogen_index]['new_symptomatic_by_variant'][variant] += len(strat.get_indices_to_track(self.sim, [index]))

        if not is_symp:
            # CASE 1: Asymptomatic: may infect others, but have no symptoms and do not die 
            dur_asym2rec = cvu.sample(**durpars['asym2rec'], size=1) * lbd 
            date_rec = self.date_p_infectious[pathogen_index, index] + dur_asym2rec 

            if date_rec <= self.t:
                date_rec = self.t + 1
                dur_asym2rec = self.t - self.date_p_infectious[pathogen_index, index]

            self.date_p_recovered[pathogen_index, index] = date_rec  # Date they recover. PxN. 
            self.dur_disease[pathogen_index, index] = self.dur_exp2inf[pathogen_index, index] + dur_asym2rec  # Store how long this person had COVID-19. Nx1.  
            self.date_p_dead[pathogen_index, index] = np.nan
            self.date_p_symptomatic[pathogen_index, index] = np.nan
            self.date_p_severe[pathogen_index, index] = np.nan
            self.date_p_critical[pathogen_index, index] = np.nan
        if is_symp:
            
            # CASE 2: Symptomatic: can either be mild, severe, or critical
            if np.isnan(self.date_p_symptomatic[pathogen_index, index]):
                #Generate date_p_dymptomatic and durs normally
                self.dur_inf2sym[pathogen_index, index] = cvu.sample(**durpars['inf2sym'], size=1) * lbd   
                self.date_p_symptomatic[pathogen_index, index] = self.date_p_infectious[pathogen_index, index] + self.dur_inf2sym[pathogen_index, index]  
            else:
                if self.date_p_infectious[pathogen_index, index] <= self.t and (self.date_p_infectious[pathogen_index, index] + self.dur_inf2sym[pathogen_index, index]) > self.t:
                    if self.t + lbd*(self.date_p_symptomatic[pathogen_index, index] - self.t) <= self.date_p_symptomatic[pathogen_index, index]:
                        self.dur_inf2sym[pathogen_index, index] = self.t - self.date_p_infectious[pathogen_index, index]
                        self.date_p_symptomatic[pathogen_index, index] = self.t + 1
                    else:
                        self.dur_inf2sym[pathogen_index, index] = (self.t - self.date_p_infectious[pathogen_index, index]) + lbd*(self.date_p_symptomatic[pathogen_index, index]-self.t)
                        self.date_p_symptomatic[pathogen_index, index] = self.date_p_infectious[pathogen_index, index] +  self.dur_inf2sym[pathogen_index, index]
            
                          
            #Scale symp_probs based on co-infection 

            is_sev = not np.isnan(self.date_p_severe[pathogen_index, index])
            prev_is_sev = is_sev
            if (not is_sev) or (self.date_p_severe[pathogen_index, index] > self.t):

                if (Msev_by_ind[index] >= 1 and not is_sev) or (Msev_by_ind[index] <= 1 and is_sev):
                    sev_prob_mult = pint.get_disease_traj_alpha(self.abs_severe_prob[pathogen_index,index], self.abs_severe_prob[pathogen_coinfecting,index], Msev_by_ind[index]) 
                    sev_probs = infect_pars['rel_severe_prob'] * self.severe_prob[pathogen_index, index]*(1-self.sev_imm[pathogen_index, variant, index]) * sev_prob_mult # Probability of these people being severe
                    sev_probs = sev_probs if sev_probs < 1 else 1 
                    is_sev = len(cvu.binomial_arr([sev_probs])) != 0 # See if they're a severe or mild case

           
            if not is_sev:
                #mild
                dur_mild2rec = cvu.sample(**durpars['mild2rec'], size=1) * lbd 
                date_rec = self.date_p_symptomatic[pathogen_index, index] + dur_mild2rec 

                if date_rec <= self.t:
                    date_rec = self.t + 1
                    dur_mild2rec = self.t - self.date_p_symptomatic[pathogen_index, index]

                self.date_p_dead[pathogen_index, index] = np.nan
                self.date_p_severe[pathogen_index, index] = np.nan
                self.date_p_critical[pathogen_index, index] = np.nan
                self.date_p_recovered[pathogen_index, index] = date_rec  # Date they recover. PxN. 
                self.dur_disease[pathogen_index, index] = self.dur_exp2inf[pathogen_index, index] + self.dur_inf2sym[pathogen_index,index] + dur_mild2rec  # Store how long this person had COVID-19. Nx1.  
                 
            if is_sev:
                
                if prev_is_sev != is_sev:
                    self.flows_variant[pathogen_index]['new_severe_by_variant'][variant] += len(strat.get_indices_to_track(self.sim, [index]))
                # CASE 2.2: Severe cases: hospitalization required, may become critical

                if np.isnan(self.date_p_severe[pathogen_index, index]):
                    #Generate date_p_dymptomatic and durs normally
                    self.dur_sym2sev[pathogen_index, index] = cvu.sample(**durpars['sym2sev'], size=1) * lbd   
                    self.date_p_severe[pathogen_index, index] = self.date_p_symptomatic[pathogen_index, index] + self.dur_sym2sev[pathogen_index, index]  
                else:
                    if self.date_p_symptomatic[pathogen_index, index] <= self.t and (self.date_p_symptomatic[pathogen_index, index] + self.dur_sym2sev[pathogen_index, index]  ) > self.t:
                        if self.t + lbd*(self.date_p_severe[pathogen_index, index] - self.t) <= self.date_p_severe[pathogen_index, index]:
                            self.dur_sym2sev[pathogen_index, index] = self.t - self.date_p_symptomatic[pathogen_index, index]
                            self.date_p_severe[pathogen_index, index] = self.t + 1
                        else:
                            self.dur_sym2sev[pathogen_index, index] = (self.t - self.date_p_symptomatic[pathogen_index, index]) + lbd*(self.date_p_severe[pathogen_index, index]-self.t)
                            self.date_p_severe[pathogen_index, index] = self.date_p_symptomatic[pathogen_index, index] +  self.dur_sym2sev[pathogen_index, index]
                             
                             
                is_crit = not np.isnan(self.date_p_critical[pathogen_index, index])

                if (not is_crit) or (self.date_p_critical[pathogen_index, index] > self.t):  
                    if (Msev_by_ind[index] >= 1 and not is_crit) or (Msev_by_ind[index] <= 1 and is_crit): 
                        crit_prob_mult = pint.get_disease_traj_alpha(self.abs_crit_prob[pathogen_index,index], self.abs_crit_prob[pathogen_coinfecting,index], Msev_by_ind[index]) 
             
                        crit_probs = infect_pars['rel_crit_prob'] * self.crit_prob[pathogen_index, index] * crit_prob_mult # TODO change symp_immm
                        crit_probs = crit_probs if crit_probs <=1 else 1

                        is_crit = len(cvu.binomial_arr([crit_probs])!= 0) # Determine if they develop symptoms 
                   
                if not is_crit:
                      
                    # CASE 2.2.1 Not critical - they will recover
                    dur_sev2rec = cvu.sample(**durpars['sev2rec'], size=1) * lbd
                    date_rec = self.date_p_severe[pathogen_index, index] + dur_sev2rec
                    
                    if date_rec <= self.t:
                        date_rec = self.t + 1
                        dur_sev2rec = self.t - self.date_p_severe[pathogen_index, index]
                    self.date_p_dead[pathogen_index, index] = np.nan
                    self.date_p_recovered[pathogen_index, index] = date_rec # Date they recover
                    self.dur_disease[pathogen_index, index] = self.dur_exp2inf[pathogen_index,index] + self.dur_inf2sym[pathogen_index,index] + self.dur_sym2sev[pathogen_index,index] + dur_sev2rec  # Store how long this person had COVID-19
                    self.date_p_critical[pathogen_index, index] = np.nan

                if is_crit: 
                     
                    if np.isnan(self.date_p_critical[pathogen_index, index]):
                    #Generate date_p_dymptomatic and durs normally
                        self.dur_sev2crit[pathogen_index, index] = cvu.sample(**durpars['sev2crit'], size=1) * lbd   
                        self.date_p_critical[pathogen_index, index] = self.date_p_severe[pathogen_index, index] + self.dur_sev2crit[pathogen_index, index]  
                    else:
                        if self.date_p_severe[pathogen_index, index] <= self.t and (self.date_p_severe[pathogen_index, index] + self.dur_sev2crit[pathogen_index, index]) > self.t:
                            if self.t + lbd*(self.date_p_critical[pathogen_index, index] - self.t) <= self.date_p_critical[pathogen_index, index]:
                                self.dur_sev2crit[pathogen_index, index] = self.t - self.date_p_severe[pathogen_index, index]
                                self.date_p_critical[pathogen_index, index] = self.t + 2
                            else:
                                self.dur_sev2crit[pathogen_index, index] = (self.t - self.date_p_severe[pathogen_index, index]) + lbd*(self.date_p_critical[pathogen_index, index]-self.t)
                                self.date_p_critical[pathogen_index, index] = self.date_p_severe[pathogen_index, index] +  self.dur_sev2crit[pathogen_index, index]
                          
                    # CASE 2.2.2: Critical cases: ICU required, may die 
                    is_dead = not np.isnan(self.date_p_dead[pathogen_index, index])

                    if (not is_dead) or (self.date_p_dead[pathogen_index, index] > self.t):  
                        if (Msev_by_ind[index] >= 1 and not is_dead) or (Msev_by_ind[index] <= 1 and is_dead): 
                            dead_prob_mult = pint.get_disease_traj_alpha(self.abs_death_prob[pathogen_index,index], self.abs_death_prob[pathogen_coinfecting,index], Msev_by_ind[index]) 
             
                            death_probs = infect_pars['rel_death_prob'] * self.death_prob[pathogen_index, index] * dead_prob_mult  
                            death_probs = death_probs if death_probs <=1 else 1

                            is_dead = len(cvu.binomial_arr([death_probs])!= 0) # Determine if they develop symptoms 
                      
                    if not is_dead:
                             
                        dur_crit2rec = cvu.sample(**durpars['crit2rec'], size=1) * lbd
                        date_rec =  self.date_p_critical[pathogen_index,index] + dur_crit2rec

                        if date_rec <= self.t:
                            date_rec = self.t + 1
                            dur_crit2rec = self.t - self.date_p_critical[pathogen_index, index]
                        self.date_p_dead[pathogen_index, index] = np.nan
                        self.date_p_recovered[pathogen_index, index] = date_rec  
                        self.dur_disease[pathogen_index, index] = self.dur_exp2inf[pathogen_index, index] + self.dur_inf2sym[pathogen_index,index] + self.dur_sym2sev[pathogen_index, index] + self.dur_sev2crit[pathogen_index, index] + dur_crit2rec  # Store how long this person had COVID-19

                    else:
                        
                        dur_crit2die = cvu.sample(**durpars['crit2die'], size=1) * lbd
                        date_dead =  self.date_p_critical[pathogen_index, index] + dur_crit2die

                        if date_dead <= self.t:
                            date_dead = self.t + 1
                            dur_crit2die = self.t - self.date_p_critical[pathogen_index, index]
                              
                        self.date_p_dead[pathogen_index, index] = date_dead
                        self.dur_disease[pathogen_index, index] = self.dur_exp2inf[pathogen_index,index] + self.dur_inf2sym[pathogen_index,index] + self.dur_sym2sev[pathogen_index,index] + self.dur_sev2crit[pathogen_index,index] + dur_crit2die   # Store how long this person had COVID-19
                        self.date_p_recovered[pathogen_index,index] = np.nan # If they did die, remove them from recovered

        
        #self.print_disease_trajectory(index, pathogen_index)
        
        # HANDLE VIRAL LOAD CONTROL POINTS
       
        # Get P_inf: where viral load crosses 10^6 cp/mL
        self.x_p_inf[pathogen_index,index] = self.dur_exp2inf[pathogen_index, index]
        self.y_p_inf[pathogen_index,index] = 6

        # Get P1: where viral load crosses 10^3 cp/mL; time difference obtained empirically through simulation
        self.x_p1[pathogen_index,index] = np.maximum(self.x_p_inf[pathogen_index,index] - (np.random.gamma(2, 0.35, size=1) + 0.25), 0)
        self.y_p1[pathogen_index,index] = 3

        # Get P2: where viral load peaks; time difference obtained empirically through simulation
        self.x_p2[pathogen_index,index] = self.x_p_inf[pathogen_index,index] + (np.random.gamma(3, 0.26, size=1) + 0.1)
        self.y_p2[pathogen_index,index] = ((self.y_p_inf[pathogen_index,index] - self.y_p1[pathogen_index,index])*(self.x_p2[pathogen_index,index] - self.x_p1[pathogen_index,index])/(self.x_p_inf[pathogen_index,index] - self.x_p1[pathogen_index,index])) + self.y_p1[pathogen_index,index]

        # Align P1, P_inf, and P2 to current time
        self.x_p1[pathogen_index,index] = self.x_p1[pathogen_index,index] + self.t
        self.x_p_inf[pathogen_index,index] = self.x_p_inf[pathogen_index,index] + self.t
        self.x_p2[pathogen_index,index] = self.x_p2[pathogen_index,index] + self.t

        # Get P3: where viral load drops below 10^6 cp/mL
        time_recovered = self.date_p_recovered[pathogen_index, index]# This is needed to make a copy
        is_dead = np.isnan(self.date_p_dead[pathogen_index, index])
        if is_dead:
            time_recovered= self.date_p_dead[pathogen_index,index]
        self.x_p3[pathogen_index,index] = np.maximum(time_recovered, self.x_p2[pathogen_index,index])
        self.y_p3[pathogen_index,index] = 6
        return

    def print_disease_trajectory(self, i, p):
        print("-> DISEASE TRAJECTORY FOR PERSON ", i, "PATHOGEN", self.sim.pathogens[p].label)
        print("exposed", self.date_p_exposed[p, i])
        print("infectious", self.date_p_infectious[p, i]) 
        print("symptomatic", self.date_p_symptomatic[p, i]) 
        print("severe", self.date_p_severe[p, i]) 
        print("critical", self.date_p_critical[p, i]) 
        print("dead", self.date_p_dead[p, i])
        print("recovered", self.date_p_recovered[p, i])



    def merge_states(self, n_pathogens, debug, merge_dates):
        '''
        Method to merge the pathogen disease states into an overral state for the person, used in interventions... 
        Not to be called by the user directly; 
        Called at the end of each sim.step() function and after infect()
        '''
        states_to_merge_with_OR = [
            'susceptible', 
            'exposed',
            'infectious',
            'symptomatic',
            'severe',
            'critical',  
            'naive',
            'dead'
            ]
        states_to_merge_with_AND = [ 
            ]

        if(debug):
            states_to_merge_with_OR = [
            'susceptible'] 
            states_to_merge_with_AND = [
                'recovered',
                ]

        self['recovered'].fill(False)
        for state in states_to_merge_with_OR:
            self[state].fill(False)

        for state in states_to_merge_with_AND:
            self[state].fill(True)
             
        #merge states
        for i in range(n_pathogens):
            
            #Set states with logical OR (states are false by default, then we set True when for one of the pathogen the state is false for a person)
            for state in states_to_merge_with_OR:
                inds_with_state = self[f'p_{state}'][i].nonzero()[0]  
                self[f'{state}'][inds_with_state] = True 
                
            #Set states with logical AND (states are true by default, then we set False when for one of the pathogen the state is false for a person)
            for state in states_to_merge_with_AND: 
                inds_with_state = np.logical_not(self[f'p_{state}'][i]).nonzero()[0]
                self[f'{state}'][inds_with_state] = False 
        

        #Overall recovered
        for i in range(n_pathogens):
            inds_rec = self[f'p_recovered'][i].nonzero()[0]  
            self['recovered'][inds_rec] = True

        inds_exposed = self['exposed'].nonzero()[0]  
        self['recovered'][inds_exposed] = False

        if merge_dates:
            for state in states_to_merge_with_OR+states_to_merge_with_AND:
                self.set_overall_state_date(state, n_pathogens, True if state in states_to_merge_with_OR else False) 
            
            self.set_overall_state_date('recovered', n_pathogens, False) 


                 
  
    def set_overall_state_date(self, state_key, n_pathogens, useMin = True):   
        v = np.copy(self[f'date_p_{state_key}'][0])
        for i in range(n_pathogens): 
            if(useMin): 
                v = cvu.custom_np_fmin(self[f'date_p_{state_key}'][i], v, size = self.pars['pop_size'])
            else: 
                v = cvu.custom_np_fmax(self[f'date_p_{state_key}'][i], v, size = self.pars['pop_size'])
          
        self[f'date_{state_key}'] = v
             



    def test(self, inds, test_sensitivity=1.0, loss_prob=0.0, test_delay=0, pathogen = 0):
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
        
        self.p_tested[pathogen,inds] = True
        self.date_p_tested[pathogen,inds] = self.t # Only keep the last time they tested

        is_infectious = cvu.itruei(self.p_infectious[pathogen], inds)
        pos_test      = cvu.n_binomial(test_sensitivity, len(is_infectious))
        is_inf_pos    = is_infectious[pos_test]

        not_diagnosed = is_inf_pos[np.isnan(self.date_diagnosed[is_inf_pos])]
        not_lost      = cvu.n_binomial(1.0-loss_prob, len(not_diagnosed))
        final_inds    = not_diagnosed[not_lost]

        # Store the date the person will be diagnosed, as well as the date they took the test which will come back positive
        self.date_diagnosed[final_inds] = self.t + test_delay
        self.date_p_diagnosed[pathogen, final_inds] = self.t + test_delay
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
                    print(f'{intro}, who had asymptomatic infection.')
                else:
                    print(f'{intro}, who had symptomatic infection.')
            else:
                print(f'{intro}, who did not contract infection.')

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
                'date_diagnosed'      : 'was diagnosed with a pathogen',
                'date_end_quarantine' : 'ended quarantine',
                'date_infectious'     : 'became infectious',
                'date_known_contact'  : 'was notified they may have been exposed to a pathogen',
                'date_pos_test'       : 'took a positive test',
                'date_quarantined'    : 'entered quarantine',
                'date_recovered'      : 'recovered',
                'date_severe'         : 'developed severe symptoms and needed hospitalization',
                'date_symptomatic'    : 'became symptomatic',
                'date_tested'         : 'was tested for pathogen',
                'date_vaccinated'     : 'was vaccinated against pathogen',
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
                        events.append((infection['date'], f'was infected with a pathogen by {infection["source"]} via the {llabel} layer'))
                    else:
                        events.append((infection['date'], 'was infected with a pathogen as a seed infection'))

                if infection['source'] == uid:
                    x = len([a for a in self.infection_log if a['source'] == infection['target']])
                    events.append((infection['date'],f'gave a pathogen to {infection["target"]} via the {llabel} layer ({x} secondary infections)'))

            if len(events):
                for day, event in sorted(events, key=lambda x: x[0]):
                    print(f'On day {day:.0f}, {uid} {event}')
            else:
                print(f'Nothing happened to {uid} during the simulation.')
        return

