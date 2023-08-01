'''
Set the defaults across each of the different files.

To change the default precision from 32 bit (default) to 64 bit, use::

    cv.options.set(precision=64)
'''

from genericpath import samefile
import numpy as np
import numba as nb
import sciris as sc
from .settings import options as cvo # To set options

# Specify all externally visible functions this file defines -- other things are available as e.g. cv.defaults.default_int
__all__ = ['default_float', 'default_int', 'get_default_colors', 'get_default_plots']


#%% Specify what data types to use

result_float = np.float64 # Always use float64 for results, for simplicity
if cvo.precision == 32:
    default_float = np.float32
    default_int   = np.int32
    nbfloat       = nb.float32
    nbint         = nb.int32
elif cvo.precision == 64: # pragma: no cover
    default_float = np.float64
    default_int   = np.int64
    nbfloat       = nb.float64
    nbint         = nb.int64
else:
    raise NotImplementedError(f'Precision must be either 32 bit or 64 bit, not {cvo.precision}')


#%% Define all properties of people

class PeopleMeta(sc.prettyobj):
    ''' For storing all the keys relating to a person and people '''

    def __init__(self):

        # Set the properties of a person
        self.person = [
            'uid',              # Int
            'age',              # Float
            'sex',              # Float
            'symp_prob',        # Float
            'severe_prob',      # Float
            'crit_prob',        # Float
            'death_prob',       # Float
            'abs_symp_prob',        # Float
            'abs_severe_prob',      # Float
            'abs_crit_prob',        # Float
            'abs_death_prob',       # Float
            'rel_trans',        # Float
            'rel_sus',          # Float
            'viral_load',       # Float
            'n_infections',     # Int
            'n_breakthroughs',  # Int
            'watch_fpr',        # Float
            'has_watch',        # Bool
            'income',           # Int
            'cons_days_in_quar',    # Int
            'cons_days_neg_rat',    # Int
            'is_coinfected'
        ]

        # Set the states that a person can be in: these are all booleans per person -- used in people.py
        self.states = [
            'susceptible',
            'naive',#remove
            'exposed',
            'infectious',
            'symptomatic',
            'severe',
            'critical',
            'tested',
            'diagnosed',
            'recovered',
            'known_dead',
            'dead',
            'known_contact',
            'quarantined',
            'vaccinated',
            'alerted',       # Used to denote that a smartwatch alert was received on the current day; date_alerted refers to the date last alerted
        ]

        #Each field would be initialized as an matrix NxK where N is the number of pathogens in the simulation, K is the number of agens in the simulation  
        self.pathogen_states = [
            'p_susceptible',    #susceptible with a specific pathogen
            'p_naive',          #naive for a specific pathogen
            'p_exposed',        #exposed by a specific pathogen
            'p_infectious',     #infectious with a specific pathogen
            'p_symptomatic',    #symptomatic with a specific pathogen
            'p_severe',         #severe symptoms due to a specific pathogen
            'p_critical',       #in critical condition caused by specific pathogen
            'p_tested',         #tested for a specific pathogen
            'p_diagnosed',      #diagnosed for a specific pathogen
            'p_recovered',      #recovered after infectious with a specific pathogen
            'p_dead',           #dead due to a specific pathogen
        ]
        
        #Information on variants of pathogens (variant states).
        self.pathogen_variants =[
            'p_exposed_variant',    #exposed to variant: matrix NxK, int
            'p_infectious_variant', #infectious with variant: matrix NxK, int
            'p_recovered_variant',  #recovered from variant: matrix NxK, int

            'p_exposed_by_variant',     #Records of exposed to a variant: matrix NxNkxK
            'p_infectious_by_variant',  #Records of infectiou with variant: matrix NxNkxK
]

         
        # Immune states, by pathogen by variant 
        self.imm_states = [
            'sus_imm',  # Float, by pathogen, by variant, by people,(Matrix NxNkxp, where N is the number of pathogens, and Nk the number of variants of N, whrere p is the population size)
            'symp_imm', # Float, by pathogen, by variant, by people,(Matrix NxNkxp, where N is the number of pathogens, and Nk the number of variants of N, whrere p is the population size)
            'sev_imm',  # Float, by pathogen, by variant, by people,(Matrix NxNkxp, where N is the number of pathogens, and Nk the number of variants of N, whrere p is the population size)
        ]

        # Neutralizing antibody states, 1D array of length: number of pathogens in simulation
        self.nab_states = [
            'peak_nab',    # Float, peak neutralization titre relative to convalescent plasma
            'nab',         # Float, current neutralization titre relative to convalescent plasma
            't_nab_event', # Int, time since nab-conferring event
        ]

        self.imm_levels = [
            'imm_level',
            't_peak_imm',
            't_min_imm', 
            'curr_min',
            ]

        # Additional vaccination states
        self.vacc_states = [
            'doses',          # Number of doses given per person
            'vaccine_source', # index of vaccine that individual received
        ]

        # Set the dates various events took place: these are floats per person -- used in people.py
        self.dates = [f'date_{state}' for state in self.states] # Convert each state into a date
        self.dates.append('date_pos_test') # Store the date when a person tested which will come back positive
        self.dates.append('date_end_quarantine') # Store the date when a person comes out of quarantine

        
        self.pathogen_dates = [f'date_{state}' for state in self.pathogen_states] # Convert each state into a date, arrays of NxP where N is num of pathogen and P is num of people
        #self.dates.append('date_pos_test') # Store the date when a person tested which will come back positive
        #self.dates.append('date_end_quarantine') # Store the date when a person comes out of quarantine

        # Duration of different states: these are floats per person -- used in people.py
       
       #each field is a 1D array of length: number of pathogens in simulation
        self.durs = [
            'dur_exp2inf',
            'dur_inf2sym',
            'dur_sym2sev',
            'dur_sev2crit',
            'dur_disease',
        ]

        # Timing of control points for viral load; each field is a 1D array of length: number of pathogens in simulation
        self.vl_points = [
            'x_p_inf',
            'y_p_inf',
            'x_p1',
            'y_p1',
            'x_p2',
            'y_p2',
            'x_p3',
            'y_p3',
        ]

        self.population_sampling = [
            'provides_sample_prob', #float 
            'IgG_level' #float
        ]

        
        self.all_states = self.person + self.states + self.pathogen_states + self.pathogen_variants+ self.imm_states + self.nab_states + self.imm_levels+ self.vacc_states + self.dates+ self.pathogen_dates + self.durs + self.vl_points + self.population_sampling

        # Validate 
        self.state_types = ['person', 'states','pathogen_states', 'pathogen_variants', 'imm_states',
                            'nab_states', 'vacc_states', 'dates', 'pathogen_dates', 'durs', 'all_states', 'imm_levels', 'population_sampling']
        for state_type in self.state_types:
            states = getattr(self, state_type)
            n_states        = len(states)
            n_unique_states = len(set(states))
            if n_states != n_unique_states: # pragma: no cover
                errormsg = f'In {state_type}, only {n_unique_states} of {n_states} state names are unique'
                raise ValueError(errormsg)

        return



#%% Define other defaults

# A subset of the above states are used for results
result_stocks = {
    'susceptible': 'Number susceptible',
    'exposed':     'Number exposed',
    'infectious':  'Number infectious',
    'symptomatic': 'Number symptomatic',
    'severe':      'Number of severe cases',
    'critical':    'Number of critical cases',
    'recovered':   'Number recovered',
    'dead':        'Number dead',
    'diagnosed':   'Number of confirmed cases',
    'known_dead':  'Number of confirmed deaths',
    'quarantined': 'Number in quarantine',
    'vaccinated':  'Number of people vaccinated',
}

result_stocks_by_variant = {
    'exposed_by_variant':    'Number exposed by variant',
    'infectious_by_variant': 'Number infectious by variant',
}

# The types of result that are counted as flows -- used in sim.py; value is the label suffix
result_flows = {
    'infections':   'infections',
    'reinfections': 'reinfections',
    'infectious':   'infectious',
    'symptomatic':  'symptomatic cases',
    'severe':       'severe cases',
    'critical':     'critical cases',
    'recoveries':   'recoveries',
    'deaths':       'deaths',
    'tests':        'tests',
    'diagnoses':    'diagnoses',
    'known_deaths': 'known deaths',
    'quarantined':  'quarantines started',
    'doses':        'vaccine doses',
    'vaccinated':   'vaccinated people'
}

result_flows_by_variant = {
    'infections_by_variant':  'infections by variant',
    'symptomatic_by_variant': 'symptomatic by variant',
    'severe_by_variant':      'severe by variant',
    'infectious_by_variant':  'infectious by variant',
    'diagnoses_by_variant':   'diagnoses by variant'
}

result_flows_smartwatches = {'alerted': 'smartwatch alerts',
                             'alerts_tp': 'true positive smartwatch alerts',
                             'alerts_fn': 'false negative smartwatch alerts',
                             'alerts_tn': 'true negative smartwatch alerts',
                             'alerts_fp': 'false positives smartwatch alerts',
                             'Q_w_i': 'smartwatch users incorrectly quarantined',
                             'Q_w': 'smartwatch users quarantined'}

result_imm = {
    'pop_nabs':       'Population average nabs',
    'pop_protection': 'Population average protective immunity',
    'pop_imm': 'Population average immunity'
}

# Define new and cumulative flows
new_result_flows = [f'new_{key}' for key in result_flows.keys()]
cum_result_flows = [f'cum_{key}' for key in result_flows.keys()]
new_result_flows_by_variant = [f'new_{key}' for key in result_flows_by_variant.keys()]
cum_result_flows_by_variant = [f'cum_{key}' for key in result_flows_by_variant.keys()]
new_result_flows_smartwatches = [f'new_{key}' for key in result_flows_smartwatches.keys()]
cum_result_flows_smartwatches = [f'cum_{key}' for key in result_flows_smartwatches.keys()]

# Parameters that can vary by variant
variant_pars = [
    'rel_beta',
    'rel_symp_prob',
    'rel_severe_prob',
    'rel_crit_prob',
    'rel_death_prob',
]

# Immunity is broken down according to 3 axes, as listed here
immunity_axes = ['sus', 'symp', 'sev']

# Immunity protection also varies depending on your infection history
immunity_sources = [
    'asymptomatic',
    'mild',
    'severe',
]

# Default age data, based on Seattle 2018 census data -- used in population.py
default_age_data = np.array([
    [ 0,  4, 0.0605],
    [ 5,  9, 0.0607],
    [10, 14, 0.0566],
    [15, 19, 0.0557],
    [20, 24, 0.0612],
    [25, 29, 0.0843],
    [30, 34, 0.0848],
    [35, 39, 0.0764],
    [40, 44, 0.0697],
    [45, 49, 0.0701],
    [50, 54, 0.0681],
    [55, 59, 0.0653],
    [60, 64, 0.0591],
    [65, 69, 0.0453],
    [70, 74, 0.0312],
    [75, 79, 0.02016], # Calculated based on 0.0504 total for >=75
    [80, 84, 0.01344],
    [85, 89, 0.01008],
    [90, 99, 0.00672],
])


def get_default_colors():
    '''
    Specify plot colors -- used in sim.py.

    NB, includes duplicates since stocks and flows are named differently.
    '''
    c = sc.objdict()
    c.susceptible           = '#4d771e'
    c.exposed               = '#c78f65'
    c.exposed_by_variant    = '#c75649'
    c.infectious            = '#e45226'
    c.infectious_by_variant = c.infectious
    c.infections            = '#b62413'
    c.reinfections          = '#732e26'
    c.infections_by_variant = '#b62413'
    c.tests                 = '#aaa8ff'
    c.diagnoses             = '#5f5cd2'
    c.diagnosed             = c.diagnoses
    c.quarantined           = '#5c399c'
    c.doses                 = c.quarantined # TODO: new color
    c.vaccinated            = c.quarantined
    c.recoveries            = '#9e1149'
    c.recovered             = c.recoveries
    c.symptomatic           = '#c1ad71'
    c.symptomatic_by_variant= c.symptomatic
    c.severe                = '#c1981d'
    c.severe_by_variant     = c.severe
    c.diagnoses_by_variant  = '#5f5cd2'
    c.critical              = '#b86113'
    c.deaths                = '#000000'
    c.dead                  = c.deaths
    c.known_dead            = c.deaths
    c.known_deaths          = c.deaths
    c.default               = '#000000'
    c.pop_nabs              = '#32733d'
    c.pop_protection        = '#9e1149'
    c.pop_symp_protection   = '#b86113'
    c.alerted               = '#c55d28' # Desert orange
    c.alerts_tp             = '#ffff00' # Yellow
    c.alerts_fn             = '#ff0000' # Red
    c.alerts_tn             = '#00ff00' # Green
    c.alerts_fp             = '#f000ff' # Magenta
    return c


# Define the 'overview plots', i.e. the most useful set of plots to explore different aspects of a simulation
overview_plots = [
    'cum_infections',
    'cum_severe',
    'cum_critical',
    'cum_deaths',
    'cum_known_deaths',
    'cum_diagnoses',
    'new_infections',
    'new_severe',
    'new_critical',
    'new_deaths',
    'new_diagnoses',
    'n_infectious',
    'n_severe',
    'n_critical',
    'n_susceptible',
    'new_tests',
    'n_symptomatic',
    'new_quarantined',
    'n_quarantined',
    'new_doses',
    'new_vaccinated',
    'cum_vaccinated',
    'cum_doses',
    'test_yield',
    'r_eff',
]

overview_variant_plots = [
    'cum_infections_by_variant',
    'new_infections_by_variant',
    'n_infectious_by_variant',
    'cum_reinfections',
    'new_reinfections',
    'pop_nabs',
    'pop_protection',
    'pop_symp_protection',
]

def get_default_plots(which='default', kind='sim', sim=None):
    '''
    Specify which quantities to plot; used in sim.py.

    Args:
        which (str): either 'default' or 'overview'
    '''
    which = str(which).lower() # To make comparisons easier

    # Check that kind makes sense
    sim_kind   = 'sim'
    scens_kind = 'scens'
    kindmap = {
        None:      sim_kind,
        'sim':     sim_kind,
        'default': sim_kind,
        'msim':    scens_kind,
        'scen':    scens_kind,
        'scens':   scens_kind,
    }
    if kind not in kindmap.keys():
        errormsg = f'Expecting "sim" or "scens", not "{kind}"'
        raise ValueError(errormsg)
    else:
        is_sim = kindmap[kind] == sim_kind

    # Default plots -- different for sims and scenarios
    if which in ['none', 'default']:

        if is_sim:
            plots = sc.odict({
                'Total counts': [
                    'cum_infections',
                    'n_infectious',
                    'cum_diagnoses',
                ],
                'Daily counts': [
                    'new_infections',
                    'new_diagnoses',
                ],
                'Health outcomes': [
                    'cum_severe',
                    'cum_critical',
                    'cum_deaths',
                    'cum_known_deaths',
                ],
            })

        else: # pragma: no cover
            plots = sc.odict({
                'Cumulative infections': [
                    'cum_infections',
                ],
                'New infections per day': [
                    'new_infections',
                ],
                'Cumulative deaths': [
                    'cum_deaths',
                ],
            })

    # Show an overview
    elif which == 'overview': # pragma: no cover
        plots = sc.dcp(overview_plots)

    # Plot absolutely everything
    elif which == 'all': # pragma: no cover
        plots = sim.result_keys('all')

    # Show an overview plus variants
    elif 'overview' in which and 'variant' in which: # pragma: no cover
        plots = sc.dcp(overview_plots) + sc.dcp(overview_variant_plots)

    # Show default but with variants
    elif which.startswith('variant'): # pragma: no cover
        if is_sim:
            plots = sc.odict({
                'Cumulative infections by variant': [
                    'cum_infections_by_variant',
                ],
                'New infections by variant': [
                    'new_infections_by_variant',
                ],
                'Health outcomes': [
                    'cum_severe',
                    'cum_critical',
                    'cum_deaths',
                ],
            })

        else: # pragma: no cover
            plots = sc.odict({
                    'Cumulative infections by variant': [
                        'cum_infections_by_variant',
                    ],
                    'New infections by variant': [
                        'new_infections_by_variant',
                    ],
                    'New diagnoses': [
                        'new_diagnoses',
                    ],
                    'Cumulative deaths': [
                        'cum_deaths',
                    ],
            })

    # Plot SEIR compartments
    elif which == 'seir': # pragma: no cover
        plots = [
            'n_susceptible',
            'n_preinfectious',
            'n_infectious',
            'n_removed',
        ],

    else: # pragma: no cover
        errormsg = f'The choice which="{which}" is not supported: choices are "default", "overview", "all", "variant", "overview-variant", or "seir", along with any result key (see sim.results_keys(\'all\') for options)'
        raise ValueError(errormsg)

    return plots
