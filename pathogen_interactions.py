import numpy as np


def mod_rel_trans(current_pathogen, p_exposed, n_pathogens, rel_trans, Mtrans):  
    '''
    Updates rel_trans based on current co-infection: updated_rel_trans(Pi)= rel_trans(Pi)*Mtrans[i,j]
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
    Upon exposure to Pi we recalculate rel_sus(Pi) based on ongoing & previous infection
    to account for innate immune response & cross immunity
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