# Classes are used to conduct COVID-19 testing on the population. 
# The classes are defined according to the three purposes of testing: Diagnostic, Screening, and Surveillance. 
from infection.test_pipeline import eligible as eg
from infection.test_pipeline import seek as sk
import numpy as np
import sys
import networkx as nx

class Diagnostic:
    # Used to conduct diagnostic testing
    # Takes as input:
    # - Type of test to be used (PCR or RAT)    (mandatory)
    #     - Test-specific details: SE/SP, 
    # - Eligibility criteria                    (optional)
    # - Test seeking probabilities              (optional)

    def __init__(self, criteria, seekprobs): 
        self.system = 'Diagnostic'
        self.criteria = criteria
        self.seekprobs = seekprobs

    def initialize(self): 
        return

    def apply(self, testobj, sim): 
        # Apply test criteria
        if not(testobj.crit_kwargs is None):
            all_eligible = eg.check_criteria(self, testobj, sim, **testobj.crit_kwargs)
        else:
            all_eligible = eg.check_criteria(self, testobj, sim)

        # Apply test seeking probabilities
        all_seekers = sk.get_all_seekers(all_eligible, self.seekprobs)

        return all_eligible, all_seekers


class Screening: 

    def __init__(self, criteria, seekprobs): 
        self.system = 'Screening'
        self.criteria = criteria
        self.seekprobs = seekprobs

    def initialize(self, pop_type, people): 
        if not(people is None) and (pop_type == 'behaviour_module'):
            self.workplaces = people.workplaces
            self.n_workplaces = people.n_workplaces
        return

    def apply(self, testobj, sim): 
        # Apply test criteria
        if not(testobj.crit_kwargs is None):
            all_eligible = eg.check_criteria(self, testobj, sim, **testobj.crit_kwargs)
        else:
            all_eligible = eg.check_criteria(self, testobj, sim)

        # Apply test seeking probabilities
        all_seekers = sk.get_all_seekers(all_eligible, self.seekprobs)

        return all_eligible, all_seekers


class Surveillance: 

    def __init__(self): 
        self.system = 'Surveillance'
        return
    
    def initialize(self): 
        return 

    def apply(self, testobj, sim): 
        return


# Helper functions
def find_workplaces(p1, p2):
    '''
    Infer all workplaces with at least 2 people. Same idea as find_workplaces_custom below, just implemented with networkx.
    '''
    G = nx.Graph()
    G.add_nodes_from(p1)
    G.add_nodes_from(p2)
    G.add_edges_from([(i,j) for i,j in zip(p1,p2)])
    components = [c for c in nx.connected_components(G)]  # Since nodes are added from p1 and p2, components all have at least two nodes

    workplaces = []
    wpid = 0
    for c in components: 
        wp = dict()
        wp['member_uids'] = np.array(set(list(c)))
        wp['wpid'] = wpid
        workplaces.append(wp)
        wpid += 1

    return workplaces


def find_workplaces_custom(p1, p2):
    '''
    Infer all workplaces given a contact network represented by two lists, by determining all connected components. Intended for use in hybrid populations, where workplaces are not specified.
    Notes:
    - Currently only works for small hybrid populations due to recursion limit constraints
    - Assume that distinct connected components correspond to distinct workplaces (i.e., assume disjoint contact networks do not exist within each workplace).
    - Synthetic populations may include workplaces consisting of a single individual. That type of workplace can't be inferred here, since the individual 
    would have no workplace contacts, and hence would not appear in the contact network. 
    - Hybrid populations have workplaces where individuals contact 20 people on average according to a Poisson distribution. It is very rare for single-individual workplaces 
    to be constructed in make_random_contacts (population.py) so the inability to detect single-individual workplaces should not matter much.
    '''
    sys.setrecursionlimit(10**6)  # This may need to be increased for larger populations
    visited = set()
    workplaces = []
    workplace_id = 0

    while len(set(p1).difference(visited)) > 0:
        old = visited.copy()
        vertex = (set(p1).difference(visited)).pop()  # Choose an unvisited vertex to begin a DFT
        DFT(vertex, visited, p1, p2)
        new = visited.difference(old)  # Newly visited vertices are exactly the members of a newly inferred workplace
        workplaces.append({'member_uids': np.array(list(new)), 'wpid': workplace_id})
        workplace_id += 1
        print(workplace_id)
        
    return workplaces

def DFT(vertex, visited, p1, p2):
    '''
    Recursive depth first traversal of a connected component in a graph. 
    '''
    if vertex in visited:
        return 
    else:
        visited.add(vertex)
        connections = set(list(p2[(p1==vertex).nonzero()[0]])).union(set(list(p1[(p2==vertex).nonzero()[0]])))  # Need to determine connections bidirectionally
        for c in connections:
            DFT(c, visited, p1, p2)

