import numpy as np
import pathosim.utils as cvu
from math import ceil

# Very simple module of what background ILI looks like

def infect_ILI(sim, pathogen = 0): 
    '''
    Void function that infects some proportion of the population with non-COVID ILI, refreshed each week. Meant to be called at each time step. We assume 
    that people cannot be infected with both COVID and ILI at the same time. Unsure if this is realistic. 
    
    Assumptions: 
        - Everyone infected is symptomatic
        - Everyone infected will remain infected for a week
        - Infections occur spontaneously, and does not involve transmission

    Args: 
        sim  (simulation) : Initialized simulation object 
        p_infect  (float) : Proportion of population to infect. For now, assume 8% of population has ILI on any given day
    '''
    if isinstance(sim.pathogens[pathogen].bkg_ILI, (float, np.floating)):
        ILI_float(sim, pathogen)
    else: 
        ILI_array(sim, pathogen)
        return


def ILI_float(sim, pathogen):
    p_infect = sim.pathogens[pathogen].bkg_ILI
    cur = sim.t
    pop_size = sim.people.pars['pop_size']
    n_infect = int(pop_size * p_infect)

    if cur == 0:  # Initialize tracker of whether people have symptomatic ILI or not
        sim.people.unlock()  # Redundancy 
        sim.people['symptomatic_ILI'] = np.full(pop_size, False)  # Whether each individual has symptomatic ILI 
        sim.people['date_symptomatic_ILI'] = np.full(pop_size, np.nan)  # Date that a person last became symptomatic with ILI 
        sim.people['date_recovered_ILI'] = np.full(pop_size, np.nan)  # Date that a person last recovered from ILI

        # Initialize infections in individuals who do not have COVID
        no_covid = cvu.true(sim.people.exposed == False)
        if len(no_covid) >= n_infect:
            infect_inds = np.random.choice(no_covid, n_infect, replace=False).astype(int)
        elif len(no_covid) == 0:
            infect_inds = np.array([], dtype=int)
        else:
            infect_inds = np.random.choice(no_covid, n_infect, replace=True).astype(int)  # everyone without covid gets selected

        # Determine how long until they recover
        recov_times = np.random.choice(np.arange(7, 11), len(infect_inds)).astype(int)  # Recover individuals in 7-10 days
        sim.people['date_symptomatic_ILI'][infect_inds] = cur
        sim.people['date_recovered_ILI'][infect_inds] = cur + recov_times
        sim.people['symptomatic_ILI'][infect_inds] = True

    else: 
        # Find candidates for new infections - not infected with COVID or ILI
        no_covid_or_ILI = cvu.true(np.logical_and(sim.people.exposed == False, sim.people['symptomatic_ILI'] == False))  

        # Check who is supposed to recover today, and recover them. This happens after previous step to prevent back-to-back infections.
        recov_today = cvu.true(sim.people['date_recovered_ILI'] == cur)
        sim.people['symptomatic_ILI'][recov_today] = False

        # Infect new people, reproduction number is not 1 if n_infect quota was not met in the last infection cycle
        if len(no_covid_or_ILI) >= n_infect:
            infect_inds = np.random.choice(no_covid_or_ILI, n_infect, replace=False).astype(int)
        elif len(no_covid_or_ILI) == 0:
            infect_inds = np.array([], dtype=int)
        else:
            infect_inds = np.random.choice(no_covid_or_ILI, n_infect, replace=True).astype(int)  # everyone without covid or ILI gets selected

        # Determine how long until they recover
        recov_times = np.random.choice(np.arange(7, 11), len(infect_inds))  # Recover individuals in 7-10 days
        sim.people['date_symptomatic_ILI'][infect_inds] = cur
        sim.people['date_recovered_ILI'][infect_inds] = cur + recov_times
        sim.people['symptomatic_ILI'][infect_inds] = True


def ILI_array(sim, pathogen):
    p_infect = sim.pathogens[pathogen].bkg_ILI
    cur = sim.t
    cur_week = ceil( (cur+1) / 7)
    pop_size = sim.people.pars['pop_size']
    n_days = sim.pars['n_days']
    recov_time = 7

    if len(p_infect) < ceil( (n_days+1) / 7):
        raise RuntimeError('Array length of weekly background ILI rates is not enough to cover the full simulation.')
    else:
        if cur == 0: 
            sim.people.unlock()
            sim.people['symptomatic_ILI'] = np.full(pop_size, False)  # Whether each individual has symptomatic ILI 
            sim.people['date_symptomatic_ILI'] = np.full(pop_size, np.nan)  # Date that a person last became symptomatic with ILI 
            sim.people['date_recovered_ILI'] = np.full(pop_size, np.nan)  # Date that a person last recovered from ILI

            n_infect = int(pop_size * p_infect[0])  # Use rate for first week
            no_covid = cvu.true(sim.people.exposed == False)
            infect_inds = np.random.choice(no_covid, n_infect, replace=False)

            sim.people['date_symptomatic_ILI'][infect_inds] = cur 
            sim.people['date_recovered_ILI'][infect_inds] = cur + recov_time
            sim.people['symptomatic_ILI'][infect_inds] = True
        elif cur % recov_time == 0:  # First day of a new week
            no_covid_or_ILI = cvu.true(np.logical_and(sim.people.exposed == False, sim.people['symptomatic_ILI'] == False))  

            sim.people['symptomatic_ILI'][:] = False  # Everyone currently infected recovers

            # Number of people to infect based on ILI prevalence of current week
            n_infect = int(pop_size * p_infect[cur_week-1])
            infect_inds = np.random.choice(no_covid_or_ILI, n_infect, replace=False)  # This could break if too many people have covid

            # Cleanup and set recovery dates
            sim.people['date_symptomatic_ILI'][infect_inds] = cur
            sim.people['date_recovered_ILI'][infect_inds] = cur + recov_time
            sim.people['symptomatic_ILI'][infect_inds] = True


def infect_ILI_old(sim, p_infect=0.04): 
    '''
    Void function that infects some proportion of the population with non-COVID ILI, refreshed each week. Meant to be called at each time step. We assume 
    that people cannot be infected with both COVID and ILI at the same time. Unsure if this is realistic. 
    
    Assumptions: 
        - Everyone infected is symptomatic
        - Everyone infected will remain infected for a week
        - Infections occur spontaneously, and does not involve transmission

    Args: 
        sim  (simulation) : Initialized simulation object 
        p_infect  (float) : Proportion of population to infect. For now, assume 8% of population has ILI on any given day
    '''
    cur = sim.t
    pop_size = sim.people.pars['pop_size']

    if cur == 0:  # Initialize tracker of whether people have symptomatic ILI or not
        sim.people.unlock()  # Redundancy 
        sim.people['symptomatic_ILI'] = np.full(pop_size, False)
        sim.people['date_symptomatic_ILI'] = np.full(pop_size, np.nan)
        sim.people['date_recovered_ILI'] = np.full(pop_size, np.nan)

    # Recover people and infect new people every 7 days
    if cur % 7 == 0: 
        # Check who is symptomatic 
        symp_ILI = cvu.true(sim.people['symptomatic_ILI'])

        # Recover people: Record current day as date of recovery, and reset symptomatic status
        sim.people['date_recovered_ILI'][symp_ILI] = cur
        sim.people['symptomatic_ILI'][symp_ILI] = False

        # Infect new people
        n_infect = int(pop_size * p_infect)
        no_cov_sym = cvu.true(sim.people.exposed == False)  # Indices of people without COVID

        if n_infect <= len(no_cov_sym): 
            ILI_inds = np.random.choice(no_cov_sym, n_infect, replace=False)  # Infect subset of people who do not have COVID
        else:
            ILI_inds = no_cov_sym  # The case where less than p=0.04 of population is free of COVID, so we infect all of them with non-COVID ILI

        sim.people['date_symptomatic_ILI'][ILI_inds] = cur
        sim.people['symptomatic_ILI'][ILI_inds] = True

    return 
