'''
Defines functions for making the population.
'''

#%% Imports
import numpy as np # Needed for a few things not provided by pl
import sciris as sc
from . import requirements as cvreq
from . import utils as cvu
from . import misc as cvm
from . import base as cvb
from . import data as cvdata
from . import defaults as cvd
from . import parameters as cvpar
from . import people as cvppl


# Specify all externally visible functions this file defines
__all__ = ['make_people', 'make_randpop', 'make_random_contacts',
           'make_microstructured_contacts']


def make_people(sim, popdict=None, workplaces = None, n_workplaces = None, reset=False, verbose=None, **kwargs):
    '''
    Make the actual people for the simulation.

    Usually called via ``sim.initialize()``. While in theory this function can be
    called directly by the user, usually it's better to call ``cv.People()`` directly.

    Args:
        sim             (Sim)  : the simulation object; population parameters are taken from the sim object
        popdict         (any)  : either None, pop.popdict, popfile, or People object
        reset           (bool) : whether to force population creation even if self.popdict/self.people exists
        verbose         (bool) : level of detail to print
        workplaces      (array): array of workplaces provided by behaviour module
        n_workplaces    (int)  : number of workplaces provided by behaviour module
        kwargs          (dict) : passed to make_randpop() or make_synthpop()

    Returns:
        people (People): people
    '''

     
    # Set inputs and defaults@
    pop_size = int(sim['pop_size']) # Shorten
    pop_type = sim['pop_type'] # Shorten
    if verbose is None:
        verbose = sim['verbose'] 


    #Handle the possible popdict input types
    if popdict is not None:  
        if sim.popfile is not None or (isinstance(popdict, cvppl.People)): 
            workplaces = popdict.workplaces;
            n_workplaces = popdict.n_workplaces; 
    
        elif sim.popfile is None: 
            popdict = parse_bhPopdict_to_ifPopdict(popdict, sim)
     
        sim.pop_type = 'behaviour_module'
    else:   
        popdict = make_randpop(sim, **kwargs)  #create random population if no population input
        sim.pop_type = 'random'
            
     
    validate_popdict(popdict, sim.pars, verbose=verbose)
     
    if sim.pars['enable_smartwatches']:
        people = cvppl.People(sim.pars, uid=popdict['uid'], age=popdict['age'], sex=popdict['sex'], contacts=popdict['contacts'], income=popdict['income'], has_watch=popdict['has_watch'], workplaces = workplaces, n_workplaces = n_workplaces) # List for storing the people
    else:
        people = cvppl.People(sim.pars, uid=popdict['uid'], age=popdict['age'], sex=popdict['sex'], contacts=popdict['contacts'], workplaces = workplaces, n_workplaces = n_workplaces)
        
    sc.printv(f'Created {pop_size} people, average age {people.age.mean():0.2f} years', 2, verbose)

    return people

 

def validate_popdict(popdict, pars, verbose=True):
    '''
    Check that the popdict is the correct type, has the correct keys, and has
    the correct length
    '''

    # Check it's the right type
    try:
        popdict.keys() # Although not used directly, this is used in the error message below, and is a good proxy for a dict-like object
    except Exception as E:
        errormsg = f'The popdict should be a dictionary or cv.People object, but instead is {type(popdict)}'
        raise TypeError(errormsg) from E

    # Check keys and lengths
    required_keys = ['uid', 'age', 'sex']
    popdict_keys = popdict.keys()
    pop_size = pars['pop_size']
    for key in required_keys:

        if key not in popdict_keys:
            errormsg = f'Could not find required key "{key}" in popdict; available keys are: {sc.strjoin(popdict.keys())}'
            sc.KeyNotFoundError(errormsg)

        actual_size = len(popdict[key])
        if actual_size != pop_size:
            errormsg = f'Could not use supplied popdict since key {key} has length {actual_size}, but all keys must have length {pop_size}'
            raise ValueError(errormsg)

        isnan = np.isnan(popdict[key]).sum()
        if isnan:
            errormsg = f'Population not fully created: {isnan:,} NaNs found in {key}. This can be caused by calling cv.People() instead of cv.make_people().'
            raise ValueError(errormsg)

    if ('contacts' not in popdict_keys) and (not hasattr(popdict, 'contacts')) and verbose:
        warnmsg = 'No contacts found. Please remember to add contacts before running the simulation.'
        cvm.warn(warnmsg)

    return


def make_randpop(pars, use_age_data=True, use_household_data=True, sex_ratio=0.5, **kwargs):
    '''
    Make a random population, with contacts.

    This function returns a "popdict" dictionary, which has the following (required) keys:

        - uid: an array of (usually consecutive) integers of length N, uniquely identifying each agent
        - age: an array of floats of length N, the age in years of each agent
        - sex: an array of integers of length N (not currently used, so does not have to be binary)
        - contacts: list of length N listing the contacts; see make_random_contacts() for details
        - layer_keys: a list of strings representing the different contact layers in the population; see make_random_contacts() for details

    Args:
        pars (dict): the parameter dictionary or simulation object
        use_age_data (bool): whether to use location-specific age data
        use_household_data (bool): whether to use location-specific household size data
        sex_ratio (float): proportion of the population that is male (not currently used)
        kwargs (dict): passed to contact creation method (e.g., make_hybrid_contacts)

    Returns:
        popdict (dict): a dictionary representing the population, with the following keys for a population of N agents with M contacts between them:
    '''

    pop_size = int(pars['pop_size']) # Number of people

    # Load age data and household demographics based on 2018 Seattle demographics by default, or country if available
    age_data = cvd.default_age_data
    location = pars['location']
    if location is not None:
        if pars['verbose']:
            print(f'Loading location-specific data for "{location}"')
        if use_age_data:
            try:
                age_data = cvdata.get_age_distribution(location)
            except ValueError as E:
                warnmsg = f'Could not load age data for requested location "{location}" ({str(E)}), using default'
                cvm.warn(warnmsg)
        if use_household_data:
            try:
                household_size = cvdata.get_household_size(location)
                if 'h' in pars['contacts']:
                    pars['contacts']['h'] = household_size - 1 # Subtract 1 because e.g. each person in a 3-person household has 2 contacts
                elif pars['verbose']:
                    keystr = ', '.join(list(pars['contacts'].keys()))
                    warnmsg = f'Not loading household size for "{location}" since no "h" key; keys are "{keystr}". Try "hybrid" population type?'
                    cvm.warn(warnmsg)
            except ValueError as E:
                if pars['verbose']>1: # These don't exist for many locations, so skip the warning by default
                    warnmsg = f'Could not load household size data for requested location "{location}" ({str(E)}), using default'
                    cvm.warn(warnmsg)

    # Handle sexes and ages
    uids           = np.arange(pop_size, dtype=cvd.default_int)
    sexes          = np.random.binomial(1, sex_ratio, pop_size)
    age_data_min   = age_data[:,0]
    age_data_max   = age_data[:,1] + 1 # Since actually e.g. 69.999
    age_data_range = age_data_max - age_data_min
    age_data_prob  = age_data[:,2]
    age_data_prob /= age_data_prob.sum() # Ensure it sums to 1
    age_bins       = cvu.n_multinomial(age_data_prob, pop_size) # Choose age bins
    ages           = age_data_min[age_bins] + age_data_range[age_bins]*np.random.random(pop_size) # Uniformly distribute within this age bin

    # Store output
    popdict = {}
    popdict['uid'] = uids
    popdict['age'] = ages
    popdict['sex'] = sexes

    # Actually create the contacts
    contacts = dict()
    for lkey,n in pars['contacts'].items():
        contacts[lkey] = make_random_contacts(pop_size, n, **kwargs)
   

    popdict['contacts']   = contacts
    popdict['layer_keys'] = list(pars['contacts'].keys())

    return popdict


def _tidy_edgelist(p1, p2, mapping):
    ''' Helper function to convert lists to arrays and optionally map arrays '''
    p1 = np.array(p1, dtype=cvd.default_int)
    p2 = np.array(p2, dtype=cvd.default_int)
    if mapping is not None:
        mapping = np.array(mapping, dtype=cvd.default_int)
        p1 = mapping[p1]
        p2 = mapping[p2]
    output = dict(p1=p1, p2=p2)
    return output


def make_random_contacts(pop_size, n, overshoot=1.2, dispersion=None, mapping=None, **kwargs):
    '''
    Make random static contacts for a single layer as an edgelist.

    Args:
        pop_size   (int)   : number of agents to create contacts between (N)
        n          (int) : the average number of contacts per person for this layer
        overshoot  (float) : to avoid needing to take multiple Poisson draws
        dispersion (float) : if not None, use a negative binomial distribution with this dispersion parameter instead of Poisson to make the contacts
        mapping    (array) : optionally map the generated indices onto new indices

    Returns:
        Dictionary of two arrays defining UIDs of the edgelist (sources and targets)

    New in 3.1.1: optimized and updated arguments.
    '''

    # Preprocessing
    pop_size = int(pop_size) # Number of people
    p1 = [] # Initialize the "sources"
    p2 = [] # Initialize the "targets"

    # Precalculate contacts
    n_all_contacts  = int(pop_size*n*overshoot) # The overshoot is used so we won't run out of contacts if the Poisson draws happen to be higher than the expected value
    all_contacts    = cvu.choose_r(max_n=pop_size, n=n_all_contacts) # Choose people at random
    if dispersion is None:
        p_count = cvu.n_poisson(n, pop_size) # Draw the number of Poisson contacts for this person
    else:
        p_count = cvu.n_neg_binomial(rate=n, dispersion=dispersion, n=pop_size) # Or, from a negative binomial
    p_count = np.array((p_count/2.0).round(), dtype=cvd.default_int)

    # Make contacts
    count = 0
    for p in range(pop_size):
        n_contacts = p_count[p]
        these_contacts = all_contacts[count:count+n_contacts] # Assign people
        count += n_contacts
        p1.extend([p]*n_contacts)
        p2.extend(these_contacts)

    # Tidy up
    output = _tidy_edgelist(p1, p2, mapping)

    return output


def make_microstructured_contacts(pop_size, cluster_size, mapping=None):
    '''
    Create microstructured contacts -- i.e. for households.

    Args:
        pop_size (int): total number of people
        cluster_size (int): the average size of each cluster (Poisson-sampled)

    New in version 3.1.1: optimized updated arguments.
    '''

    # Preprocessing -- same as above
    pop_size = int(pop_size) # Number of people
    p1 = [] # Initialize the "sources"
    p2 = [] # Initialize the "targets"

    # Initialize
    n_remaining = pop_size # Make clusters - each person belongs to one cluster. TODO: cluster = family? 

    # Loop over the clusters
    cluster_id = -1
    while n_remaining > 0:
        cluster_id += 1 # Assign cluster id
        this_cluster =  cvu.poisson(cluster_size)  # Sample the cluster size
        if this_cluster > n_remaining:
            this_cluster = n_remaining

        # Indices of people in this cluster
        cluster_indices = (pop_size-n_remaining) + np.arange(this_cluster)
        for source in cluster_indices: # Add symmetric pairwise contacts in each cluster
            targets = set()
            for target in cluster_indices:
                if target > source:
                    targets.add(target)
            p1.extend([source]*len(targets))
            p2.extend(list(targets))

        n_remaining -= this_cluster

    # Tidy up
    output = _tidy_edgelist(p1, p2, mapping)

    return output

 
def parse_bhPopdict_to_ifPopdict(bhPopdict, sim):
     
    '''
    Converts behaviour popdict into infection popdict
    Arguments: behaviour popdict (pop.popdict)
    '''
    if sim.pars['enable_multiregion']:
        try:
            import behaviour as bh # Optional import
        except ModuleNotFoundError as E: # pragma: no cover
            errormsg = 'Please install the behaviour module first' # Also caught in make_people()
            raise ModuleNotFoundError(errormsg) from E
        layer_mapping = None
        # Handle layer mapping
        default_layer_mapping = {'H':'h', 'S':'s', 'W':'w', 'C':'c', 'LTCF':'l'} # Remap keys from old names to new names
        layer_mapping = sc.mergedicts(default_layer_mapping, layer_mapping)
        people = bh.people.make_people(popdict=bhPopdict, pop_type='behaviour_module') # Note: veeery slow. Just a parsing function; contacts and stuff already made. 
         
        people.contacts = cvb.Contacts(**people.contacts)
        people['layer_keys'] = list(layer_mapping.values())
        return people

    bhPopdict = sc.dcp(bhPopdict)

    # Handle layer mapping
    layer_mapping = {'H':'h', 'S':'s', 'W':'w', 'C':'c', 'LTCF':'l'}

    # Create the basic lists
    pop_size = len(bhPopdict)
    uids, ages, sexes, have_watches, contacts, incomes = [], [], [], [], [], []
    for uid,person in bhPopdict.items():
        uids.append(uid)
        ages.append(person['age'])
        have_watches.append(person['has_watch'])
        sexes.append(person['sex'])
        incomes.append(person['hhincome'])

    # Replace contact UIDs with ints
    uid_mapping = {uid:u for u,uid in enumerate(uids)}
    for uid in uids:
        iid = uid_mapping[uid] # Integer UID
        person = bhPopdict.pop(uid)
        uid_contacts = sc.dcp(person['contacts'])
        int_contacts = {}
        for spkey in uid_contacts.keys():
            try:
                lkey = layer_mapping[spkey] # Map the SynthPops key into a Covasim layer key
            except KeyError: # pragma: no cover
                errormsg = f'Could not find key "{spkey}" in layer mapping "{layer_mapping}"'
                raise sc.KeyNotFoundError(errormsg)
            int_contacts[lkey] = []
            for cid in uid_contacts[spkey]: # Contact ID
                icid = uid_mapping[cid] # Integer contact ID
                if icid>iid: # Don't add duplicate contacts
                    int_contacts[lkey].append(icid)
            int_contacts[lkey] = np.array(int_contacts[lkey], dtype= np.int64)
        contacts.append(int_contacts)

    # Finalize
    popdict = {}
    popdict['uid']        = np.array(list(uid_mapping.values()), dtype= np.int64)
    popdict['age']        = np.array(ages)
    popdict['sex']        = np.array(sexes)
    popdict['has_watch']  = np.array(have_watches)
    popdict['contacts']   = sc.dcp(contacts)
    popdict['layer_keys'] = list(layer_mapping.values())
    popdict['income']    = np.array(incomes)

    return popdict