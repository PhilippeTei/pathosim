import sciris as sc
import numpy as np
from . import utils as cvu
import scipy.stats as stats

class RVSampler:
    def __init__(self, dist_name, dist_pars):

        # Check if it's a valid distribution name, from scipy. 
        if dist_name not in stats.distributions.__all__:
            raise ValueError('Invalid distribution name: %s' % dist_name)
        self.dist_name = dist_name
        self.dist_pars = dist_pars

    def sample(self, n=1):
        return stats.distributions.__dict__[self.dist_name].rvs(size=n, **self.dist_pars)/100 # It outputs percentages. Convert to fractions.


class Watches:
    def __init__(self, pars):
        self.make_pars(pars)

    def make_pars(self, pars=None):

        # Configure default parameter values
        watch_pars = sc.objdict()
        watch_pars['mean_fpr']          = 0.08345302977643946 # Mean P(non-inf alert) from Alavi et al
        watch_pars['use_variable_fpr']  = True

        watch_pars['day_i']             = np.arange(-21, 22, 1) # Dates the the above pdf is fitted wrt to.

        watch_pars['loc']               = 3.25 # Day of max probability of alert, relative to the day of symptom onset.
        watch_pars['alpha']             = 1 # Scales the probability of receiving an alert
        watch_pars['usage_rate']        = 1 # Out of people who have smartwatches, the amount who use download the alerting app and stick with it.

        # Override the default values if any others are passed in
        if pars is not None:
            for par_name, par_value in pars.items(): 
                if par_name not in watch_pars:
                    raise KeyError(f"Parameter {par_name} doesn't exist, check watches.py:make_pars()!")
                watch_pars[par_name] = par_value

        # Check if parameter values are valid
        # e.g., check that alpha >= 0

        # Calculate probability distributions. 
        # p_alert isn't a distribution, so we rescale the pdf.
        watch_pars['p_alert']           = (stats.norm(watch_pars['loc'], 6).pdf(watch_pars['day_i'])*1314 + 20)*watch_pars['mean_fpr']/20
        watch_pars['p_inf_alert']       = 1 - (1-watch_pars['p_alert'])/(1-watch_pars['mean_fpr'])

        # Scale alerting probability, if a scaling factor is provided
        if watch_pars['alpha'] != 1:
            p_alert, p_noninf_alert, p_inf_alert = scale_p_alert(watch_pars['alpha'], watch_pars['mean_fpr'], watch_pars['p_alert'])

            # Update P(non-inf alert)
            watch_pars['mean_fpr'] = p_noninf_alert

            # Update P(inf alert | day i)
            watch_pars['p_inf_alert'] = p_inf_alert

        # Assign people false positive rates
        watch_pars['fpr_sampler'] = RVSampler('expon', dict(loc=0, scale=watch_pars['mean_fpr']*100))

        self.pars = watch_pars
        return


    def send_alerts_baseline(self, sim):
        """
        Models the baseline level of alerts, i.e., alerts stemming from the false positive rate of the watch. 
        This comes from noisy non-modelled events, like agents getting non-COVID illnesses, getting stressed, over-drinking, etc. 
        """
        
        # Get indices of people with watches
        watch_inds = sim.people.has_watch.nonzero()[0]

        # Alert these people based on their false positive rate (i.e., determine if factors unrelated to infection will prompt an alert)
        alert_inds = cvu.binomial_arr(sim.people.watch_fpr[watch_inds])
        alert_inds = alert_inds.nonzero()[0] # The indicies of valid people that will get alerts
        to_alert = watch_inds[alert_inds]

        # Mark alerted individuals as alerted
        sim.people.alerted[to_alert] = True
        sim.people.date_alerted[to_alert] = sim.t
        return


    def send_alerts_infected(self, sim):
        """
        Called by sim.py:step. For infected agents, adds additional probability of being 
        alerted. Separate treatment for symptomatic and asymptomatic cases. 
        """

        # Identify exposed individuals that also have watches
        inf_watches = sim.people.has_watch * sim.people.exposed
        inds_inf_watches = np.nonzero(inf_watches)[0]
        
        # Get each individual's day relative to symptom onset
        days_rel_symp = self.get_days_rel_symp(sim, inds_inf_watches)

        # Get each individual's alert probabilities
        alert_probs = self.get_p_inf_alert_on_days(self.pars['p_inf_alert'], days_rel_symp)

        # print("DEBUG: SMARTWATCH ALERT PROBS 2")
        # # Reset daily. 
        # sim.people.all_peoples_infection_days = np.empty(20000) # hardcoded
        # sim.people.all_peoples_infection_days[:] = np.nan
        # sim.people.all_peoples_infection_days[inf_watches] = days_rel_symp
        # #### END DEBUG

        # Distribute alertsbased on the probability
        bin_mask = cvu.binomial_arr(alert_probs)
        alert_inds = inds_inf_watches[bin_mask > 0]
        sim.people.alerted[alert_inds] = True
        sim.people.date_alerted[alert_inds] = sim.t # Results in logical OR when used in conjunction with send_alerts_baseline.
        return


    def get_days_rel_symp(self, sim, inds):
        """
        Return agent-wise number of days relative to symptom onset, for both eventually symptomatic and non-symptomatic agents. 
        Get a Mx1 array of the number of days into the infection valid people are. 
        M is the # valid people: those with watches and infections.
        """
        days_rel_symp = np.zeros(len(inds))

        # isnan returns Mx1 bool array. 
        mask_asymp = np.isnan(sim.people.date_symptomatic[inds])
        mask_symp = ~mask_asymp
        
        # Returns list of indicies of inds
        mask_asymp = np.nonzero(mask_asymp)[0]
        mask_symp = np.nonzero(mask_symp)[0]

        # The inds that are symp or asymp, rel to Nx1 uids. 
        inds_asymp = inds[mask_asymp]
        inds_symp = inds[mask_symp]

        # Find dates for symp people
        days_rel_symp[mask_symp] = sim.t - sim.people.date_symptomatic[inds_symp]

        # Find dates for asymp people.
        # First, sample a hypothetical symptom onset date. 
        num_asymp = len(inds_asymp)
        dur_inf2sym = cvu.sample(**sim.people.pars['dur']['inf2sym'], size=num_asymp)

        t_fake_symp = sim.people.date_infectious[inds_asymp] + dur_inf2sym
        days_rel_symp[mask_asymp] = sim.t - t_fake_symp
        
        return days_rel_symp


    def get_p_inf_alert_on_days(self, p_inf_alert, days_rel_symp):
        '''
        Retrieves the probability of an alert driven by infection for a list of days relative to symptom onset.
        '''

        # Convert days to indices; account for people who are infected for very long or have a very long incubation period
        idx = days_rel_symp + (len(p_inf_alert)-1)//2
        idx[idx < 0] = 0
        idx[idx >= len(p_inf_alert)] = len(p_inf_alert)-1

        # Index the probabilities
        return p_inf_alert[idx.astype(int)]


    def update_alert_histories(self, sim):
        """
        Update the state of the last M days of alerts; agentwise. 
        This history is used for modelling behaviour. 
        """
        sim.people.alert_histories = np.roll(sim.people.alert_histories, shift=-1, axis=1)
        sim.people.alert_histories[:, sim.people.LEN_ALERT_HIST - 1] = sim.people.alerted
        return


# class SwBehaviour:
#     def __init__(self, pars):
#         self.make_pars(pars)

#     def make_pars(self, pars=None):
#         # Set defaults.
#         params = sc.objdict()
#         params['behaviour_policy'] = 'SQ' # Possibilities: simple quarantine (SQ), ...
#         # todo: params.behaviour_policy = sc.objdict(simple_pq=0.1, ...)
#         params['simple_pq'] = 1 # p quarantine for simple behaviour  model. 
#         params['simple_lq'] = 7 # length of quarantine for simple behaviour model.

#         # Override defaults.
#         if pars is not None:
#             for par_name, par_value in pars.items():
#                 if par_name not in params:
#                     raise KeyError(f"Parameter {par_name} doesn't exist, check SWBehaviour:make_pars()!")
#                 params[par_name] = par_value
        
#         self.pars = params
    
#     def consume_alerts(self, sim, policy="SQ"):
#         """
#         The person is going to look at their alert and do something. 
#         This is the master function for the behaviour modelling.
#         For now, we just have people quarantine with 10% chance.

#         Policies:
#         SQ: Simple quarantine, with a 10% chance. 
#         """
#         policy = self.pars['behaviour_policy']
#         # Get the people who have alerts. 
#         alerted = sim.people.true('sw_alarmed') # TODO: THIS IS OUTDATED; WE NO LONGER RESET IT HERE BUT IN CHECK_ALERTED.
#         # Set those alert flags to False.
#         sim.people.sw_alarmed[alerted] = False

#         if policy == "SQ":
#             # Quarantine. 
#             self.simple_quarantine(sim, alerted)
#         else:
#             raise ValueError(f"Policy {policy} does not exist!")

#     def simple_quarantine(self, sim, alerted_inds):
#         p_Q = self.pars.simple_pq
#         l_Q = self.pars.simple_lq

#         # Set p_Q% of input people to quarantine.
#         quar_draws = cvu.n_binomial(p_Q, len(alerted_inds))
#         quar_inds = quar_draws.nonzero()[0] # The indicies that will quarantine

#         to_quar = alerted_inds[quar_inds]

#         sim.people.schedule_quarantine(to_quar, start_date=sim.t, period=l_Q)  # Schedule quarantine for the notified people to start on the date they will be notified


# Function that takes in a scaling factor, false positive rate, and p_alert and returns the three profiles
def scale_p_alert(factor, fpr, p_alert_base):
    """
    Analogous to sweeping sensitivity/specificity tradeoff via adjusting classifier logit threshold tau. 
        factor: the scaling alpha. 
        fpr: mean_p_noninf_alert
        p_alert_base: the total p_alert combining infected and noninfected causes. 
    """

    # Calculate new P(alert) by multiplying by a scaling factor
    p_alert = np.minimum(factor*p_alert_base, 1)

    # Calculate corresponding P(inf alert)
    if np.minimum(factor*fpr, 1) == 1:
        p_inf_alert = np.ones(len(p_alert_base))
    else:
        p_inf_alert = 1 - (1-p_alert)/(1-(factor*fpr))
    
    return p_alert, np.minimum(factor*fpr, 1), p_inf_alert