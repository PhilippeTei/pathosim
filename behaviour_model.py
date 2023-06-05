from . import utils as cvu
from . import people_categories as cvpc

def send_to_quarantine(sim, inds, l_Q):
    """
    Full name: State transition - quarantine
    Args: 
        length of quarantine
    """

    sim.people.schedule_quarantine(inds, start_date=sim.t, period=l_Q)

#### Code some mock-up functions to simplify development. ####
def random_tests(sim):
    # Randomly select 5% of people to seek tests; assume instant results. 
    test_inds = cvu.choose(sim.n, 0.05)

    # Filter them by who's exposed. 100% SP. 
    test_inds = test_inds[sim.people.exposed[test_inds]]

    # Set their test date to today. 100% SE. Assume 0 day delay. 
    sim.people.date_pos_test[test_inds] = sim.t
    sim.people.date_diagnosed[test_inds] = sim.t

class BehaviourUpdater:
    def __init__(self, **kwargs):
        self.make_pars()
        self.update_pars(**kwargs)
    
    def make_pars(self):
        self.quarantine_probs = {
            'cs_sx': 0.3,
            'pos_test': 0.8,
            'sw_alert': 0.05
        }
        self.l_Q = 7 # quarantine length. 

    def update_behaviour(self, sim):
        """
        Simple model: each data point gives people a certain probability of quarantine. 
        """
        ppl_cats = cvpc.PeopleCategories(sim, criteria=self.quarantine_probs.keys())
        all_quar_inds = set()

        for cat in ppl_cats.cat_dict:
            # Draw quarantines from the category. Coin flip. 
            quar_draws = cvu.n_binomial(self.quarantine_probs[cat], len(ppl_cats.cat_dict[cat]))

            # Get the indicies corresponding to quar_draws being true.
            quar_inds = ppl_cats.cat_dict[cat][quar_draws.nonzero()[0]]
            all_quar_inds = all_quar_inds.update(quar_inds)

        # Schedule quarantines.
        send_to_quarantine(sim, all_quar_inds, self.l_Q)

    def update_pars(self, **kwargs):
        """
        Args: 
            **kwargs: keyword arguments to update parameters.
        Body:
            Update fields with parameters. New fields cannot be initialized in this way. 
        """
        for k, v in kwargs.items():
            # Check if k is already a field. 
            if hasattr(self, k):
                setattr(self, k, v)
            else:
                raise ValueError('Invalid parameter: %s' % k)