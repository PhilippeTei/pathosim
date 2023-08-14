'''
Set the parameters for Covasim.
'''

import numpy as np
import sciris as sc
from .settings import options as cvo # For setting global options
from . import misc as cvm
from . import defaults as cvd
from . import pathogens as pat

__all__ = ['make_pars', 'reset_layer_pars', 'get_vaccine_choices',
            'get_vaccine_variant_pars', 'get_vaccine_dose_pars']


def make_pars(version=None, **kwargs):
    '''
    Create the parameters for the simulation. Typically, this function is used
    internally rather than called by the user; e.g. typical use would be to do
    sim = cv.Sim() and then inspect sim.pars, rather than calling this function
    directly.

    Args:
        set_prognoses (bool): whether or not to create prognoses (else, added when the population is created)
        prog_by_age   (bool): whether or not to use age-based severity, mortality etc.
        kwargs        (dict): any additional kwargs are interpreted as parameter names
        version       (str):  if supplied, use parameters from this Covasim version

    Returns:
        pars (dict): the parameters of the simulation
    '''
    pars = {}

    #Multi-pathogen
    #pars['pathogens'] = pat.Pathogen('SARS-CoV-2')
    pars['n_pathogens'] = 1 #Number of pathogens circulating in the simulation
    pars['pathogens'] = []
    #----------------- ADDITIONAL MODULES (not in baseline multi-pathogen sim)-----------------#
    # Multi-region
    pars['enable_multiregion'] = False
    pars['multiregion'] = None

    # Calculate stratified statistics
    pars['enable_stratifications'] = False
    pars['stratification_pars'] = None 

    # Smartwatch par dictionaries
    pars['enable_smartwatches'] = False
    pars['smartwatch_pars']     = None

    # Testing par dictionaries. Note: This is used to toggle symptoms and bkg ILI, since we only need symptoms and bkg ILI when we need testing. 
    pars['enable_testobjs'] = False 
    pars['testobjs'] = None  # This naming convention is possibly ambiguous - do NOT supply actually test objects (do that for pars['testing'])

    #------------------------------------------------------------------------------------------#
    #Pathogen-pathogen interaction
    pars['Mtrans']= None # Mtrans[i,j]=P(transmit Pi | co�infected with Pj) / P(transmit Pi at baseline)
    pars['Miimm'] = None # P(infection with Pi | co�infected with Pj) / P(infection with Pi at baseline)
    pars['Mcimm'] = None # contribution of current immunity to Pj to immunity to Pi
    pars['Mdur'] =  None # (duration of Pi | co�infected with Pj) / (duration of Pi at baseline)
    pars['Msev'] =  None # (severity of Pi | co�infected with Pj) / (severity of Pi at baseline)
   
    #--------------------------- POPULATION & CONTACTS PARAMETERS------------------------------#
    # Population parameters
    pars['pop_size']     = 20e3     # Number of agents, i.e., people susceptible to SARS-CoV-2
    pars['pop_type']     = 'random' # What type of population data to use -- 'random' (fastest), 'synthpops' (best), 'hybrid' (compromise)
    pars['location']     = None     # What location to load data from -- default Seattle
    pars['pop']          = None     # Detailed synthetic population -- this is supplied for access to workplace information for surveillance

    # Simulation parameters
    pars['start_day']  = '2020-03-01' # Start day of the simulation
    pars['end_day']    = None         # End day of the simulation
    pars['n_days']     = 60           # Number of days to run, if end_day isn't specified
    pars['rand_seed']  = 1            # Random seed, if None, don't reset
    pars['verbose']    = cvo.verbose  # Whether or not to display information during the run -- options are 0 (silent), 0.1 (some; default), 1 (default), 2 (everything)

    # Rescaling parameters
    pars['pop_scale']         = 1    # Factor by which to scale the population -- e.g. pop_scale=10 with pop_size=100e3 means a population of 1 million
    pars['scaled_pop']        = None # The total scaled population, i.e. the number of agents times the scale factor
    pars['rescale']           = True # Enable dynamic rescaling of the population -- starts with pop_scale=1 and scales up dynamically as the epidemic grows
    pars['rescale_threshold'] = 0.05 # Fraction susceptible population that will trigger rescaling if rescaling
    pars['rescale_factor']    = 1.2  # Factor by which the population is rescaled on each step
    pars['frac_susceptible']  = 1.0  # What proportion of the population is susceptible to infection
    
    # Network parameters, generally initialized after the population has been constructed
    pars['contacts']        = None  # The number of contacts per layer; set by reset_layer_pars() below
    pars['dynam_layer']     = None  # Which layers are dynamic; set by reset_layer_pars() below
    pars['beta_layer']      = None  # Transmissibility per layer; set by reset_layer_pars() below

     
    #-------------------------------------------------------------------------------------------#
    pars['use_waning']   = True # Whether to use dynamically calculated immunity
     
    # Symptoms: symptom prevalence within the population
    pars['prev_COVID'] = None
    pars['prev_ILI'] = None

    # Efficacy of protection measures
    pars['iso_factor']   = None # Multiply beta by this factor for diagnosed cases to represent isolation; set by reset_layer_pars() below
    pars['quar_factor']  = None # Quarantine multiplier on transmissibility and susceptibility; set by reset_layer_pars() below
    pars['quar_period']  = 14   # Number of days to quarantine for; assumption based on standard policies

    # Parameters that govern actions people take in response to a result or notification
    pars['enable_behaviour'] = False                                 # Flag to use this module
    pars['behaviour_pars'] = {'symptom_quar_pars': {0: (10, 0)},     # Symptomatic individuals quarantine for v1 days with probability v2 in key:(v1,v2) where key is the largest key such that key <= t
                              'contact_quar_pars': {0: (10, 0)},     # Close contacts quarantine for v1 days with probability v2 in key:(v1,v2) where key is the largest key such that key <= t
                              'hh_symp_quar_pars': {0: (10, 0)},     # People with symptomatic household members quarantine for v1 days with probability v2 in key:(v1,v2) where key is the largest key such that key <= t
                              'enable_symp_early_end': {0: False},      # If enabled, those who have tested negative on RAT for more than two consecutive days can leave quarantine
                              'enable_hh_symp_exemption': {0: False},   # If enabled, household members that are younger than 18 OR vaccinated are exempt
                              'enable_contact_exemption': {0: False},   # If enabled, contacts that are younger than 18 OR vaccinated are exempt
                              'smartwatch_behaviour': {0: (False, 'stan')}} # Quarantine according to a standardized ('stan') or historical ('hist') model of behaviour applies to alerted smartwatch users (True/False) as of the most recent key
    pars['alert_behav_dist_stan'] = {'n':5, 'p':0.2, 's':10}
    pars['alert_behav_dist_hist'] = {1: {'n':5, 'p':0.05, 's':20},
                                     2: {'n':5, 'p':0.30, 's':20},
                                     3: {'n':5, 'p':0.50, 's':20}}

    # Events and interventions
    pars['interventions'] = []   # The interventions present in this simulation; populated by the user
    pars['surveillance'] = []    # The surveillance systems present in this simualtion; populated by the user
    pars['testing'] = []         # The testing systems present in this simualtion; populated by the user. These can be populated externally, or internally by supplying a parameter dictionary. 
    pars['analyzers']     = []   # Custom analysis functions; populated by the user
    pars['timelimit']     = None # Time limit for the simulation (seconds)
    pars['stopping_func'] = None # A function to call to stop the sim partway through

    # Testing
    pars['seek_multiplier'] = 1  # Uniform multiplier on test seeking probabilties
    pars['symp_prev_multiplier'] = 1  # Uniform multiplier on prevalence of different symptoms in the population

    # Health system parameters
    pars['n_beds_hosp']    = None # The number of hospital (adult acute care) beds available for severely ill patients (default is no constraint)
    pars['n_beds_icu']     = None # The number of ICU beds available for critically ill patients (default is no constraint)
    pars['no_hosp_factor'] = 2.0  # Multiplier for how much more likely severely ill people are to become critical if no hospital beds are available
    pars['no_icu_factor']  = 2.0  # Multiplier for how much more likely critically ill people are to die if no ICU beds are available

    # Handle vaccine and variant parameters
    pars['vaccine_pars'] = {} # Vaccines that are being used; populated during initialization
    pars['vaccine_map']  = {} #Reverse mapping from number to vaccine key 

    #NEW PARAMETERS FOR SURVEILLANCE ARCHITECTURE (CK_SV)

    # Surveilance parameters
    pars['enable_surveillance'] = False
    pars['surveillance_test_size'] = None
    pars['surveillance_test_percent'] = None
    pars['surveillance_viral_threshold'] = 10 #to be calibrated
    pars['surveillance_num_threshold'] = 1
    pars['surveillance_time_to_confirmation'] = 0 #to be calibrated

    #contact-based surveillance
    pars['enable_contact_testing'] = False
    pars['contact_percentile_lower'] = 75
    pars['contact_percentile_upper'] = 100
    pars['contact_test_frequency'] = 1
    pars['contact_test_start_date'] = None
    pars['contact_test_end_date'] = pars['end_day']
    
    #severity-based surveillance
    pars['enable_severity_testing'] = False
    pars['severity_test_frequency'] = 1
    pars['severity_test_start_date'] = None
    pars['severity_test_end_date'] = pars['end_day']
    pars['severity_symp_prob_threshold'] = 0.8

    #age-based surveillance
    pars['enable_age_testing'] = False
    pars['surveillance_age_lower'] = 0
    pars['surveillance_age_upper'] = 100
    pars['age_test_frequency'] = 1
    pars['age_test_start_date'] = None
    pars['age_test_end_date'] = pars['end_day']

    #random testing
    pars['enable_random_testing'] = False
    pars['random_test_frequency'] = 1
    pars['random_test_start_date'] = None
    pars['random_test_end_date'] = pars['end_day']

    #syndromic surveillance
    pars['enable_syndromic_testing'] = False
    pars['hospital_capacity_percent'] = 1
    pars['syndromic_test_percent'] = 0.25*0.3 #to be calibrated_lowprio
    pars['syndromic_days_to_test'] = 0
    pars['syndromic_time_to_confirmation'] = 0

    #metagenomics parameters
    pars['rna_depletion_enabled'] = False
    pars['pooling_enabled'] = False
    pars['pool_size'] = 1

    # Update with any supplied parameter values and generate things that need to be generated
    pars.update(kwargs)
    reset_layer_pars(pars)
 
    return pars


# Define which parameters need to be specified as a dictionary by layer -- define here so it's available at the module level for sim.py
layer_pars = ['beta_layer', 'contacts', 'dynam_layer', 'iso_factor', 'quar_factor']


def reset_layer_pars(pars, layer_keys=None, force=False):
    '''
    Helper function to set layer-specific parameters. If layer keys are not provided,
    then set them based on the population type. This function is not usually called
    directly by the user, although it can sometimes be used to fix layer key mismatches
    (i.e. if the contact layers in the population do not match the parameters). More
    commonly, however, mismatches need to be fixed explicitly.

    Args:
        pars (dict): the parameters dictionary
        layer_keys (list): the layer keys of the population, if available
        force (bool): reset the parameters even if they already exist
    '''

    # Specify defaults for random -- layer 'a' for 'all'
    layer_defaults = {}
    layer_defaults['random'] = dict(
        beta_layer  = dict(a=1.0), # Default beta
        contacts    = dict(a=20),  # Default number of contacts
        dynam_layer = dict(a=0),   # Do not use dynamic layers by default
        iso_factor  = dict(a=0.2), # Assumed isolation factor
        quar_factor = dict(a=0.3), # Assumed quarantine factor
    )

    # Specify defaults for hybrid -- household, school, work, and community layers (h, s, w, c)
    layer_defaults['hybrid'] = dict(
        beta_layer  = dict(h=3.0, s=0.6, w=0.6, c=0.3),  # Per-population beta weights; relative; in part based on Table S14 of https://science.sciencemag.org/content/sci/suppl/2020/04/28/science.abb8001.DC1/abb8001_Zhang_SM.pdf
        contacts    = dict(h=2.0, s=20,  w=16,  c=20),   # Number of contacts per person per day, estimated
        dynam_layer = dict(h=0,   s=0,   w=0,   c=0),    # Which layers are dynamic -- none by default
        iso_factor  = dict(h=0.3, s=0.1, w=0.1, c=0.1),  # Multiply beta by this factor for people in isolation
        quar_factor = dict(h=0.6, s=0.2, w=0.2, c=0.2),  # Multiply beta by this factor for people in quarantine
    )

    # Specify defaults for SynthPops -- same as hybrid except for LTCF layer (l)
    l_pars = dict(beta_layer=1.5, contacts=10, dynam_layer=0, iso_factor=0.2, quar_factor=0.3)
    layer_defaults['synthpops'] = sc.dcp(layer_defaults['hybrid'])
    for key,val in l_pars.items():
        layer_defaults['synthpops'][key]['l'] = val # Add LTCF layer params, since this isn't in hybrid. 

    # Behaviour module is the same except for quarantine factors. 
    layer_defaults['behaviour_module'] = sc.dcp(layer_defaults['synthpops'])
    b_quar_factors = dict(h=0.5, s=0.1, w=0.1, c=0.1, l=0.5)
    for key, val in b_quar_factors.items():
        layer_defaults['behaviour_module']['quar_factor'][key] = val

    # Choose the parameter defaults based on the population type, and get the layer keys
    try:
        defaults = layer_defaults[pars['pop_type']]
    except Exception as E:
        errormsg = f'Cannot load defaults for population type "{pars["pop_type"]}": must be hybrid, random, or synthpops'
        raise ValueError(errormsg) from E
    default_layer_keys = list(defaults['beta_layer'].keys()) # All layers should be the same, but use beta_layer for convenience

    # Actually set the parameters
    for pkey in layer_pars:
        par = {} # Initialize this parameter
        default_val = layer_defaults['random'][pkey]['a'] # Get the default value for this parameter

        # If forcing, we overwrite any existing parameter values
        if force:
            par_dict = defaults[pkey] # Just use defaults
        else:
            par_dict = sc.mergedicts(defaults[pkey], pars.get(pkey, None)) # Use user-supplied parameters if available, else default

        # Figure out what the layer keys for this parameter are (may be different between parameters)
        if layer_keys:
            par_layer_keys = layer_keys # Use supplied layer keys
        else:
            par_layer_keys = list(sc.odict.fromkeys(default_layer_keys + list(par_dict.keys())))  # If not supplied, use the defaults, plus any extra from the par_dict; adapted from https://www.askpython.com/python/remove-duplicate-elements-from-list-python

        # Construct this parameter, layer by layer
        for lkey in par_layer_keys: # Loop over layers
            par[lkey] = par_dict.get(lkey, default_val) # Get the value for this layer if available, else use the default for random
        pars[pkey] = par # Save this parameter to the dictionary

    return
 
 
def get_vaccine_choices():
    '''
    Define valid pre-defined vaccine names
    '''
    # List of choices currently available: new ones can be added to the list along with their aliases
    choices = {
        'default': ['default', None],
        'pfizer':  ['pfizer', 'biontech', 'pfizer-biontech', 'pf', 'pfz', 'pz', 'bnt162b2', 'comirnaty'],
        'moderna': ['moderna', 'md', 'spikevax'],
        'novavax': ['novavax', 'nova', 'covovax', 'nvx', 'nv'],
        'az':      ['astrazeneca', 'az', 'covishield', 'oxford', 'vaxzevria'],
        'jj':      ['jnj', 'johnson & johnson', 'janssen', 'jj'],
        'sinovac': ['sinovac', 'coronavac'],
        'sinopharm': ['sinopharm']
    }
    mapping = {name:key for key,synonyms in choices.items() for name in synonyms} # Flip from key:value to value:key
    return choices, mapping


def _get_from_pars(pars, default=False, key=None, defaultkey='default'):
    ''' Helper function to get the right output from vaccine and variant functions '''

    # If a string was provided, interpret it as a key and swap
    if isinstance(default, str):
        key, default = default, key

    # Handle output
    if key is not None:
        try:
            return pars[key]
        except Exception as E:
            errormsg = f'Key "{key}" not found; choices are: {sc.strjoin(pars.keys())}'
            raise sc.KeyNotFoundError(errormsg) from E
    elif default:
        return pars[defaultkey]
    else:
        return pars

    
def get_vaccine_variant_pars(default=False, vaccine=None):
    '''
    Define the effectiveness of each vaccine against each variant
    '''
    pars = dict(

        default = dict(
            wild  = 1.0,
            alpha = 1.0,
            beta  = 1.0,
            gamma = 1.0,
            delta = 1.0,
        ),

        pfizer = dict(
            wild  = 1.0,
            alpha = 1/2.0, # https://www.nejm.org/doi/full/10.1056/nejmc2100362
            beta  = 1/10.3, # https://www.nejm.org/doi/full/10.1056/nejmc2100362
            gamma = 1/6.7, # https://www.nejm.org/doi/full/10.1056/nejmc2100362
            delta = 1/2.9, # https://www.researchsquare.com/article/rs-637724/v1
            omicron = 26/802 # https://www.science.org/doi/full/10.1126/sciimmunol.abo2202
        ),

        moderna = dict(
            wild  = 1.0,
            alpha = 1/1.8,
            beta  = 1/4.5,
            gamma = 1/8.6, # https://www.nejm.org/doi/full/10.1056/nejmc2100362
            delta = 1/2.9,  # https://www.researchsquare.com/article/rs-637724/v1
            omicron = 58/1997 # https://www.science.org/doi/full/10.1126/sciimmunol.abo2202
        ),

        az = dict(
            wild  = 1.0,
            alpha = 1/2.3,
            beta  = 1/9,
            gamma = 1/2.9,
            delta = 1/6.2,  # https://www.researchsquare.com/article/rs-637724/v1
            omicron = 20/411  # https://www.science.org/doi/full/10.1126/sciimmunol.abo2202
        ),

        jj = dict(
            wild  = 1.0,
            alpha = 1.0,
            beta  = 1/3.6,  # https://www.biorxiv.org/content/10.1101/2021.07.01.450707v1.full.pdf
            gamma = 1/3.4,  # https://www.biorxiv.org/content/10.1101/2021.07.01.450707v1.full.pdf
            delta = 1/1.6,  # https://www.biorxiv.org/content/10.1101/2021.07.01.450707v1.full.pdf
            omicron = 15/163 # https://www.science.org/doi/full/10.1126/sciimmunol.abo2202
        ),

        novavax = dict( # Data from https://ir.novavax.com/news-releases/news-release-details/novavax-covid-19-vaccine-demonstrates-893-efficacy-uk-phase-3
            wild  = 1.0,
            alpha = 1/1.12,
            beta  = 1/4.7,
            gamma = 1/8.6, # Assumption, no data available yet
            delta = 1/6.2, # Assumption, no data available yet
        ),

        sinovac = dict(
            wild  = 1.0,
            alpha = 1/1.12,
            beta  = 1/4.7,
            gamma = 1/8.6, # Assumption, no data available yet
            delta = 1/6.2, # Assumption, no data available yet
        ),

        sinopharm = dict(
            wild  = 1.0,
            alpha = 1/1.12,
            beta  = 1/4.7,
            gamma = 1/8.6, # Assumption, no data available yet
            delta = 1/6.2, # Assumption, no data available yet
        )
    )

    return _get_from_pars(pars, default=default, key=vaccine)


def get_vaccine_dose_pars(default=False, vaccine=None):
    '''
    Define the parameters for each vaccine
    '''

    pars = dict(

        default = dict(
            nab_init  = dict(dist='normal', par1=0, par2=2), # Initial distribution of NAbs
            nab_boost = 2, # Factor by which a dose increases existing NABs
            doses     = 1, # Number of doses for this vaccine
            interval  = None, # Interval between doses
        ),

        pfizer = dict(
            nab_init  = dict(dist='normal', par1=-1, par2=2),
            nab_boost = 4,
            doses     = 2,
            interval  = 21,
        ),

        moderna = dict(
            nab_init  = dict(dist='normal', par1=-1, par2=2),
            nab_boost = 8,
            doses     = 2,
            interval  = 28,
        ),

        az = dict(
            nab_init  = dict(dist='normal', par1=-1.5, par2=2),
            nab_boost = 2,
            doses     = 2,
            interval  = 21,
        ),

        jj = dict(
            nab_init  = dict(dist='normal', par1=1, par2=2),
            nab_boost = 3,
            doses     = 1,
            interval  = None,
        ),

        novavax = dict(
            nab_init  = dict(dist='normal', par1=-0.9, par2=2),
            nab_boost = 3,
            doses     = 2,
            interval  = 21,
        ),

        sinovac = dict(
            nab_init  = dict(dist='normal', par1=-2, par2=2),
            nab_boost = 2,
            doses     = 2,
            interval  = 14,
        ),

        sinopharm = dict(
            nab_init  = dict(dist='normal', par1=-1, par2=2),
            nab_boost = 2,
            doses     = 2,
            interval  = 21,
        )
    )

    return _get_from_pars(pars, default, key=vaccine)