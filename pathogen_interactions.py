import numpy as np
import math

def mod_rel_trans(current_pathogen, p_exposed, n_pathogens, rel_trans, Mtrans):  
    '''
    Function that updates relative transmissibility rel_trans of the people based on current co-infection.

    Parameters:
        current_pathogen(int) = index of the pathogen co-infecting
        p_exposed(bool 2D array) = Array of whether or not a person at index j has been exposed to the pathogen with index i.
        n_pathogen(int) = number of pathogens in the simulation
        rel_trans(float array) = current relative transmissibilities of the people (float array of size # people) against the pathogen with index current_pathogen.
        Mtrans(float matrix) = matrix of relative transmissibility change due to co-infection


    '''

    if n_pathogens == 1:
        return rel_trans

    for i in range(n_pathogens):
        if i == current_pathogen:
            continue
        indices_with_i = p_exposed[i].nonzero()[0]
        rel_trans[indices_with_i] *= Mtrans[current_pathogen, i]

    return rel_trans
 

def mod_rel_sus(current_pathogen, rel_sus, p_exposed, Miimm, Mcimm, people_sus_imm, n_pathogens, pop_size):
    '''
    Function that updates relative susceptibility rel_sus of the people based on current co-infection.

    Parameters:
        current_pathogen(int) = index of the pathogen co-infecting
        rel_sus(float array) = current relative susceptibilities of the people (float array of size # people) against the pathogen with index current_pathogen.
        p_exposed(bool 2D array) = Array of whether or not a person at index j has been exposed to the pathogen with index i.
        Miimm(float matrix) = matrix of change in relative susceptibility due to innate immunity
        Mcimm(float matrix) = matrix of change in relative susceptibility due to cross-immunity
        people_sus_imm(float matrix) = reduction in susceptibility due to current immunity matrix as a scale from 0-1.
        n_pathogen(int) = number of pathogens in the simulation
        pop_size(int) = population size
    '''

    if n_pathogens == 1:
        return rel_sus

    #calculate sus_imm (innate immunity provided by ongoing infection with another pathogen (other than current_pathogen))
    sus_iimm = np.full(pop_size, 1, dtype = float)

    for i in range(n_pathogens):
        if i == current_pathogen:
            continue
        indices_with_i = p_exposed[i].nonzero()[0]
        sus_iimm[indices_with_i] = Miimm[current_pathogen, i]
           
    #calculate sus_cimm (cross immunity provided by prior infection)
    sus_cimm = np.full(pop_size, 1, dtype = float)

    for i in range(n_pathogens):
        if i == current_pathogen:
            continue
        indices_not_currently_inf = np.logical_not(p_exposed[i]).nonzero()[0]
        sus_cimm[indices_not_currently_inf] *= (1- Mcimm[current_pathogen, i] * people_sus_imm[i,0][indices_not_currently_inf])

    #Calculate updated relative susceptibility to Pi:
    rel_sus *= sus_iimm * sus_cimm
    return rel_sus

def get_disease_traj_alpha(pi, pj, Msev): 
    '''
    Function that computes the alpha factor used in modification of disease trajectory due to co-infection. 
    Detailed explanation of the use of alpha and this system is outlined in the documentation.
    '''
    p_s_pi = pi
    p_s_pj = pj
    Msev_ij = Msev

    p_s_pi_pj = 1 - ((1-p_s_pi)*(1-p_s_pj)) 

    beta = Msev_ij*p_s_pi_pj
    beta = beta if beta <= 1 else 1

    alpha1 = ((p_s_pi+p_s_pj) + np.sqrt( (p_s_pi+p_s_pj)**2 - 4*p_s_pi*p_s_pj*beta ))/(2*p_s_pi*p_s_pj) 
    alpha2 = ((p_s_pi+p_s_pj) - np.sqrt( (p_s_pi+p_s_pj)**2 - 4*p_s_pi*p_s_pj*beta ))/(2*p_s_pi*p_s_pj) 
     
    res =  min(alpha1, alpha2) 
    return res


 