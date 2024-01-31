'''
Defines classes and methods for calculating immunity
'''

import numpy as np
import numba as nb
import sciris as sc
from . import utils as cvu
from . import defaults as cvd
from . import parameters as cvpar
from . import interventions as cvi
 
__all__ = []

#%% Neutralizing antibody methods


def update_peak_nab(people, inds, nab_pars, symp=None, pathogen = 0):
    '''
    Update peak NAb level

    This function updates the peak NAb level for individuals when a NAb event occurs.
        - individuals that already have NAbs from a previous vaccination/infection have their NAb level boosted;
        - individuals without prior NAbs are assigned an initial level drawn from a distribution. This level
            depends on whether the NAbs are from a natural infection (and if so, on the infection's severity)
            or from a vaccination (and if so, on the type of vaccine).

    Args:
        people: A people object
        inds: Array of people indices
        nab_pars: Parameters from which to draw values for quantities like ['nab_init'] - either sim pars (for natural immunity) or vaccine pars
        symp: either None (if NAbs are vaccine-derived), or a dictionary keyed by 'asymp', 'mild', and 'sev' giving the indices of people with each of those symptoms
        pathogen: index of the pathogen that we are modifying immunity of

    Returns: None
    '''

    # Extract parameters and indices
    pars = people.pars
    has_nabs = people.nab[pathogen, inds] > 0
    no_prior_nab_inds = inds[~has_nabs]
    prior_nab_inds = inds[has_nabs]

    # 1) Individuals that already have NAbs from a previous vaccination/infection have their NAb level boosted
    if len(prior_nab_inds):
        boost_factor = nab_pars.nab_boost
        people.peak_nab[pathogen,prior_nab_inds] *= boost_factor


    # 2) Individuals without prior NAbs are assigned an initial level drawn from a distribution.
    if len(no_prior_nab_inds):

        # Firstly, ensure that we don't try to apply a booster effect to people without NAbs
        if nab_pars.nab_init is None:
            errormsg = f'Attempt to administer a vaccine without an initial NAb distribution to {len(no_prior_nab_inds)} unvaccinated people failed.'
            raise ValueError(errormsg)

        # Now draw the initial NAb levels
        init_nab = cvu.sample(**nab_pars.nab_init, size=len(no_prior_nab_inds))
        no_prior_nab = (2 ** init_nab)

        # Next, these initial NAb levels are normalized to be equivalent to "vaccine NAbs".
        # This is done so that when we check immunity, we can calculate immune protection
        # using a single curve and account for multiple sources of immunity (vaccine and natural).
        if symp is not None:
            # Setting up for symptom scaling
            prior_symp = np.full(pars['pop_size'], np.nan)
            prior_symp[symp['asymp']] = nab_pars.rel_imm_symp['asymp']
            prior_symp[symp['mild']] =  nab_pars.rel_imm_symp['mild']
            prior_symp[symp['sev']] =  nab_pars.rel_imm_symp['severe']
            prior_symp[prior_nab_inds] = np.nan
            prior_symp = prior_symp[~np.isnan(prior_symp)]
            # Applying symptom scaling and a normalization factor to the NAbs
            norm_factor = 1 + nab_pars.nab_eff['alpha_inf_diff']
            no_prior_nab = no_prior_nab * prior_symp * norm_factor

        # Update people's peak NAbs
        people.peak_nab[pathogen, no_prior_nab_inds] = no_prior_nab

    # Update time of nab event
    people.t_nab_event[pathogen,inds] = people.t

    return


def update_nab(people, inds, pathogen):
    '''
    Step NAb levels forward in time
    '''
    t_since_boost = people.t - people.t_nab_event[pathogen, inds]
    people.nab[pathogen,inds] += people.pars['pathogens'][pathogen].nab_kin[t_since_boost]*people.peak_nab[pathogen,inds]
    people.nab[pathogen,inds] = np.where(people.nab[pathogen,inds]<0, 0, people.nab[pathogen,inds]) # Make sure nabs don't drop below 0
    people.nab[pathogen,inds] = np.where([people.nab[pathogen,inds] > people.peak_nab[pathogen,inds]], people.peak_nab[pathogen,inds], people.nab[pathogen,inds]) # Make sure nabs don't exceed peak_nab
    return
 
def update_imm(people, inds, pathogen, min_imm, max_imm, days_to_min, days_to_max):
    '''
    Step imm levels forward in time (for pathogens using the generalized immunity system)
    '''  
    people.imm_level[pathogen],people['curr_min'][pathogen] = update_imm_nb(people['p_dead'][pathogen], people['date_p_recovered'][pathogen], people['date_p_dead'][pathogen], people.t_peak_imm[pathogen], people.t_min_imm[pathogen], people['p_exposed'][pathogen], people.t, people['date_p_exposed'][pathogen],  people['curr_min'][pathogen],people['p_recovered'][pathogen],  people.imm_level[pathogen],inds,min_imm, max_imm,days_to_min,days_to_max)
    return

 
@nb.njit()
def update_imm_nb(people_dead, people_date_rec, people_date_dead, people_t_peak, people_t_min, people_exposed, t, people_date_exposed, people_curr_min, people_recovered, people_imm_level, inds, min_imm, max_imm, days_to_min, days_to_max):

    for i in inds:    
        if people_dead[i]:
            continue

        if not np.isnan(people_date_rec[i]):
             rec_or_dead_date =  people_date_rec[i]  
        else:
             rec_or_dead_date =  people_date_dead[i]  


        if (people_t_peak[i] - t) > 0 and people_exposed[i] == True:
            x = t - people_date_exposed[i]
             
            people_imm_level[i] = immunity_growth_function(x, people_curr_min[i], max_imm, days_to_max) #if previously infected, immunity starts at the min value 



        elif (people_t_peak[i] - t) <= 0 and people_exposed[i] == True and rec_or_dead_date > t: #if between peak and recovered
            people_curr_min[i] = min_imm
            continue

        elif (people_t_min[i] - t) > 0 and (people_t_peak[i] - t) <= 0  and rec_or_dead_date <= t and people_recovered[i] == True:
            
            x = t - rec_or_dead_date
            people_imm_level[i] = immunity_decay_function(x, min_imm, max_imm, days_to_min)
            people_curr_min[i] = people_imm_level[i]  

    return people_imm_level, people_curr_min


def validate_imm(i, pathogen, people, mini, maxi, rec_dead_date): 
    if people['p_dead'][pathogen, i]:
        return

    if people.t_peak_imm[pathogen, i] == people.t:
        assert abs(people.imm_level[pathogen, i] - maxi) < 0.01
    if people.t_min_imm[pathogen, i] == people.t: 
         
        assert abs(people.imm_level[pathogen,i] - mini) < 0.01 
    
    if people.t < rec_dead_date and people.t >= people.t_peak_imm[pathogen, i]: 
        assert abs(people.imm_level[pathogen, i] - maxi) < 0.01
         

def update_IgG(people, indices, ratio, b, prev_IgG, t):
    '''
    Step IgG levels forward in time 
    ''' 

    mask = (prev_IgG == 0) & (people['nab'][0, indices] > 0)
    num_to_initialize = np.count_nonzero(mask)

    if num_to_initialize > 0:
        # Initialize IgG using uniform distribution for selected indices
        initialized_indices = indices[mask]
        #people['IgG_level'][initialized_indices] = np.random.uniform(40.63, 162.5, num_to_initialize)
        people['IgG_level'][initialized_indices] = np.random.uniform(162, 162.5, num_to_initialize)
        #average_initialized_IgG = np.mean(people['IgG_level'][initialized_indices])
        #print("Average initialized IgG:", average_initialized_IgG)

    
    # Update IgG using the formula for all other indices
    updated_indices = indices[~mask]
    if len(updated_indices) > 0:
        ratio_mask = np.isin(indices, updated_indices)
        updated_ratio = ratio[ratio_mask]  # Select ratio values for updated indices
        curr_IgG = people['IgG_level'][updated_indices]
        people['IgG_level'][updated_indices] = (updated_ratio * (curr_IgG - b)) + b
        #average_updated_IgG = np.mean(people['IgG_level'][updated_indices])
        #print("Average updated IgG:", average_updated_IgG)

        #print(people['IgG_level'][updated_indices])
    

def calc_VE(nab, ax, pars):
    '''
        Convert NAb levels to immunity protection factors, using the functional form
        given in this paper: https://doi.org/10.1101/2021.03.09.21252641

        Args:
            nab  (arr)  : an array of effective NAb levels (i.e. actual NAb levels, scaled by cross-immunity)
            ax   (str)  : axis of protection; can be 'sus', 'symp' or 'sev', corresponding to the efficacy of protection against infection, symptoms, and severe disease respectively
            pars (dict) : dictionary of parameters for the vaccine efficacy

        Returns:
            an array the same size as NAb, containing the immunity protection factors for the specified axis
         '''

    choices = ['sus', 'symp', 'sev']
    if ax not in choices:
        errormsg = f'Choice {ax} not in list of choices: {sc.strjoin(choices)}'
        raise ValueError(errormsg)

    if ax == 'sus':
        alpha = pars['alpha_inf']
        beta = pars['beta_inf']
    elif ax == 'symp':
        alpha = pars['alpha_symp_inf']
        beta = pars['beta_symp_inf']
    else:
        alpha = pars['alpha_sev_symp']
        beta = pars['beta_sev_symp']

    exp_lo = np.exp(alpha) * nab**beta
    output = exp_lo/(1+exp_lo) # Inverse logit function
    return output


def calc_VE_symp(nab, pars):
    '''
    Converts NAbs to marginal VE against symptomatic disease
    '''

    exp_lo_inf = np.exp(pars['alpha_inf']) * nab**pars['beta_inf']
    inv_lo_inf = exp_lo_inf / (1 + exp_lo_inf)

    exp_lo_symp_inf = np.exp(pars['alpha_symp_inf']) * nab**pars['beta_symp_inf']
    inv_lo_symp_inf = exp_lo_symp_inf / (1 + exp_lo_symp_inf)

    VE_symp = 1 - ((1 - inv_lo_inf)*(1 - inv_lo_symp_inf))
    return VE_symp
 


# %% Immunity methods

def init_immunity(sim, create=False):
    ''' Initialize immunity matrices with all variants that will eventually be in the sim'''

    # Don't use this function if immunity is turned off
    if not sim['use_waning']:
        return
     
    for p in sim.pathogens:
        # Next, precompute the NAb kinetics and store these for access during the sim
        if p.use_nab_framework:
            p.nab_kin = precompute_waning(length=sim.npts, pars=p.nab_decay) 

    return


def check_immunity(people, variant, pathogen):
    '''
    Calculate people's immunity on this timestep from prior infections + vaccination. Calculates effective NAbs by
    weighting individuals NAbs by source and then calculating efficacy.

    There are two fundamental sources of immunity:

           (1) prior exposure: degree of protection depends on variant, prior symptoms, and time since recovery
           (2) vaccination: degree of protection depends on variant, vaccine, and time since vaccination

    '''

    # Handle parameters and indices
    pars = people.pars['pathogens'][pathogen]
    v_cross_imm = pars.immunity[variant,:] # cross-immunity/own-immunity scalars to be applied to NAb level before computing efficacy
    v_cross_imm_multiplier = np.ones(len(people))
     
    if pars.use_nab_framework:
        nab_eff = pars.nab_eff
        current_nabs = sc.dcp(people.nab[pathogen]) 
    else:
        current_imm = sc.dcp(people.imm_level[pathogen])

    date_rec = people.date_p_recovered[pathogen]  # Date recovered
    is_vacc = cvu.true(people.vaccinated)  # Vaccinated
    vacc_source = people.vaccine_source[is_vacc]

    was_inf = cvu.true((people.t >= people.date_p_recovered[pathogen])& (people.dead == False))  # Had a previous exposure, now recovered
    was_inf_same = cvu.true((people.p_recovered_variant[pathogen] == variant) & (people.t >= date_rec) & (people.dead== False))  # Had a previous exposure to the same variant, now recovered

    was_inf_diff = np.setdiff1d(was_inf, was_inf_same)  # Had a previous exposure to a different variant, now recovered
    variant_was_inf_diff = people.p_recovered_variant[pathogen, was_inf_diff]
       
    #print(f'{people.date_p_recovered[pathogen,2489]} with variant {people.p_recovered_variant[pathogen, 2489]} with pathogen {pathogen}') 
    variant_was_inf_diff = variant_was_inf_diff.astype(cvd.default_int)
      
    v_cross_imm_multiplier[was_inf_same] = v_cross_imm[variant]
    v_cross_imm_multiplier[was_inf_diff] = [v_cross_imm[i] for i in variant_was_inf_diff]

    pars = people.pars


    if pars['pathogens'][pathogen].use_nab_framework:

        if len(is_vacc) and len(pars['vaccine_pars']): # if using simple_vaccine, do not apply
            vx_pars = pars['vaccine_pars']
            vx_map = pars['vaccine_map']
            var_key = people.pars['pathogens'][pathogen].get_variants_labels()[variant]
            imm_arr = np.zeros(max(vx_map.keys())+1)
             
            for num,key in vx_map.items():
                if vx_pars[key]['type'] == 'gen': #ignore if vax changes generic immunity system
                    continue
                imm_arr[num] = vx_pars[key][var_key]

            v_cross_imm_multiplier[is_vacc] = imm_arr[vacc_source]

        current_nabs *= v_cross_imm_multiplier
        people.sus_imm[pathogen,variant,:]  = calc_VE(current_nabs, 'sus',  nab_eff)
        people.symp_imm[pathogen,variant,:] = calc_VE(current_nabs, 'symp', nab_eff)
        people.sev_imm[pathogen, variant,:]  = calc_VE(current_nabs, 'sev',  nab_eff)
    else:
         
        current_imm *= v_cross_imm_multiplier 
        clamped_current_imm = cvu.clamp_np_arr(current_imm, 0, pars['pathogens'][pathogen].imm_peak)
         
        people.sus_imm[pathogen,variant,:]  = clamped_current_imm
        people.symp_imm[pathogen,variant,:] = clamped_current_imm
        people.sev_imm[pathogen, variant,:] = clamped_current_imm
        

    return



#%% Methods for computing waning

def precompute_waning(length, pars=None):
    '''
    Process functional form and parameters into values:

        - 'nab_growth_decay' : based on Khoury et al. (https://www.nature.com/articles/s41591-021-01377-8)
        - 'nab_decay'   : specific decay function taken from https://doi.org/10.1101/2021.03.09.21252641
        - 'exp_decay'   : exponential decay. Parameters should be init_val and half_life (half_life can be None/nan)
        - 'linear_decay': linear decay

    A custom function can also be supplied.

    Args:
        length (float): length of array to return, i.e., for how long waning is calculated
        pars (dict): passed to individual immunity functions

    Returns:
        array of length 'length' of values
    '''

    pars = sc.dcp(pars)
    form = pars.pop('form')
    choices = [
        'nab_growth_decay', # Default if no form is provided
        'nab_decay',
        'exp_decay',
    ]

    # Process inputs
    if form is None or form == 'nab_growth_decay':
        output = nab_growth_decay(length, **pars)

    elif form == 'nab_decay':
        output = nab_decay(length, **pars)

    elif form == 'exp_decay':
        if pars['half_life'] is None: pars['half_life'] = np.nan
        output = exp_decay(length, **pars)

    elif callable(form):
        output = form(length, **pars)

    else:
        errormsg = f'The selected functional form "{form}" is not implemented; choices are: {sc.strjoin(choices)}'
        raise NotImplementedError(errormsg)

    return output


def nab_growth_decay(length, growth_time, decay_rate1, decay_time1, decay_rate2, decay_time2):
    '''
    Returns an array of length 'length' containing the evaluated function nab growth/decay
    function at each point.

    Uses linear growth + exponential decay, with the rate of exponential decay also set to
    decay linearly until it reaches a 10-year half life.

    Args:
        length (int): number of points
        growth_time (int): length of time NAbs grow (used to determine slope)
        decay_rate1 (float): initial rate of exponential decay
        decay_time1 (float): time of the first exponential decay
        decay_rate2 (float): the rate of exponential decay in late period
        decay_time2 (float): total time until late decay period (must be greater than decay_time1)
    '''


    def f1(t, growth_time):
        '''Simple linear growth'''
        return (1 / growth_time) * t

    def f2(t, decay_time1, decay_time2, decay_rate1, decay_rate2):
        decayRate = np.full(len(t), fill_value=decay_rate1)
        decayRate[cvu.true(t>decay_time2)] = decay_rate2
        slowing = (1 / (decay_time2 - decay_time1)) * (decay_rate1 - decay_rate2)
        decayRate[cvu.true((t>decay_time1)*(t<=decay_time2))] = decay_rate1 - slowing * (np.arange(len(cvu.true((t>decay_time1)*(t<=decay_time2))), dtype=cvd.default_int))
        titre = np.zeros(len(t))
        for i in range(1, len(t)):
            titre[i] = titre[i-1]+decayRate[i]
        return np.exp(-titre)

    if decay_time2 < decay_time1:
        errormsg = f'Decay time 2 must be larger than decay time 1, but you supplied {decay_time2} which is smaller than {decay_time1}.'
        raise ValueError(errormsg)

    length = length + 1
    t1 = np.arange(growth_time, dtype=cvd.default_int)
    t2 = np.arange(length - growth_time, dtype=cvd.default_int)
    y1 = f1(t1, growth_time)
    y2 = f2(t2, decay_time1, decay_time2, decay_rate1, decay_rate2)
    y  = np.concatenate([y1,y2])
    y = np.diff(y)[0:length]

    return y


def nab_decay(length, decay_rate1, decay_time1, decay_rate2):
    '''
    Returns an array of length 'length' containing the evaluated function nab decay
    function at each point.

    Uses exponential decay, with the rate of exponential decay also set to exponentially
    decay (!) after 250 days.

    Args:
        length (int): number of points
        decay_rate1 (float): initial rate of exponential decay
        decay_time1 (float): time on the first exponential decay
        decay_rate2 (float): the rate at which the decay decays
    '''
    def f1(t, decay_rate1):
        ''' Simple exponential decay '''
        return np.exp(-t*decay_rate1)

    def f2(t, decay_rate1, decay_time1, decay_rate2):
        ''' Complex exponential decay '''
        return np.exp(-t*(decay_rate1*np.exp(-(t-decay_time1)*decay_rate2)))

    t  = np.arange(length, dtype=cvd.default_int)
    y1 = f1(cvu.true(t<=decay_time1), decay_rate1)
    y2 = f2(cvu.true(t>decay_time1), decay_rate1, decay_time1, decay_rate2)
    y  = np.concatenate([[-np.inf],y1,y2])
    y = np.diff(y)[0:length]
    y[0] = 1
    return y


def exp_decay(length, init_val, half_life, delay=None):
    '''
    Returns an array of length t with values for the immunity at each time step after recovery
    '''
    length = length+1
    decay_rate = np.log(2) / half_life if ~np.isnan(half_life) else 0.
    if delay is not None:
        t = np.arange(length-delay, dtype=cvd.default_int)
        growth = linear_growth(delay, init_val/delay)
        decay = init_val * np.exp(-decay_rate * t)
        result = np.concatenate([growth, decay], axis=None)
    else:
        t = np.arange(length, dtype=cvd.default_int)
        result = init_val * np.exp(-decay_rate * t)
    return np.diff(result)


def linear_decay(length, init_val, slope):
    ''' Calculate linear decay '''
    result = -slope*np.ones(length)
    result[0] = init_val
    return result


def linear_growth(length, slope):
    ''' Calculate linear growth '''
    return slope*np.ones(length)

@nb.njit()
def immunity_growth_function(x, min, max, peak_t):
    '''Calculate immunity growth function for generalized immunity system'''
    y = (max-min) / np.exp(peak_t) * (np.exp(x)-1)  + min
    return y

@nb.njit()
def immunity_decay_function(x, min, max, min_t): # https://en.wikipedia.org/wiki/Non-analytic_smooth_function#Smooth_transition_functions
    '''Calculate immunity decay function for generalized immunity system'''
    a = exponeover(x/min_t)
    y = (max-min) * (1 - a / ( a + exponeover(1-(x/min_t)))) + min
    return y

@nb.njit()
def exponeover(x): 
    y = 0
    if (x> 0.00000001):
        y = np.exp(-1/x)
    return y


