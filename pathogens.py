import numpy as np
import sciris as sc
from . import utils as utls
from . import defaults as df
from . import parameters as ptpar
from . import interventions as itvs

__all__ = ['Pathogen']

class Pathogen(sc.prettyobj):
    
    class Variant(sc.prettyobj):
        '''
        Add a new variant to the sim

        Args:
            variant (str/dict): name of variant, or dictionary of parameters specifying information about the variant
            days   (int/list): day(s) on which new variant is introduced
            label       (str): if variant is supplied as a dict, the name of the variant
            n_imports   (int): the number of imports of the variant to be added
            rescale    (bool): whether the number of imports should be rescaled with the population

        **Example**::

            alpha    = cv.variant('alpha', days=10) # Make the alpha variant B117 active from day 10
            p1      = cv.variant('p1', days=15) # Make variant P1 active from day 15
            my_var  = cv.variant(variant={'rel_beta': 2.5}, label='My variant', days=20)
            sim     = cv.Sim(variants=[alpha, p1, my_var]).run() # Add them all to the sim
            sim2    = cv.Sim(variants=cv.variant('alpha', days=0, n_imports=20), pop_infected=0).run() # Replace default variant with alpha
        '''

        def __init__(self, variant, days, label=None, n_imports=1, rescale=True):
            self.days = days # Handle inputs
            self.n_imports = int(n_imports)
            self.rescale   = rescale
            self.index     = None # Index of the variant in the sim; set later
            self.label     = None # Variant label (used as a dict key)
            self.p         = None # This is where the parameters will be stored
            self.parse(variant=variant, label=label) # Variants can be defined in different ways: process these here
            self.initialized = False
            return


        def parse(self, variant=None, label=None):
            ''' Unpack variant information, which may be given as either a string or a dict '''

            # Option 1: variants can be chosen from a list of pre-defined variants
            if isinstance(variant, str):

                choices, mapping = ptpar.get_variant_choices()
                known_variant_pars = ptpar.get_variant_pars()

                label = variant.lower()
                for txt in ['.', ' ', 'variant', 'variant', 'voc']:
                    label = label.replace(txt, '')

                if label in mapping:
                    label = mapping[label]
                    variant_pars = known_variant_pars[label]
                else:
                    errormsg = f'The selected variant "{variant}" is not implemented; choices are:\n{sc.pp(choices, doprint=False)}'
                    raise NotImplementedError(errormsg)

            # Option 2: variants can be specified as a dict of pars
            elif isinstance(variant, dict):

                default_variant_pars = ptpar.get_variant_pars(default=True)
                default_keys = list(default_variant_pars.keys())

                # Parse label
                variant_pars = variant
                label = variant_pars.pop('label', label) # Allow including the label in the parameters
                if label is None:
                    label = 'custom'

                # Check that valid keys have been supplied...
                invalid = []
                for key in variant_pars.keys():
                    if key not in default_keys:
                        invalid.append(key)
                if len(invalid):
                    errormsg = f'Could not parse variant keys "{sc.strjoin(invalid)}"; valid keys are: "{sc.strjoin(df.variant_pars)}"'
                    raise sc.KeyNotFoundError(errormsg)

                # ...and populate any that are missing
                for key in default_keys:
                    if key not in variant_pars:
                        variant_pars[key] = default_variant_pars[key]

            else:
                errormsg = f'Could not understand {type(variant)}, please specify as a dict or a predefined variant:\n{sc.pp(choices, doprint=False)}'
                raise ValueError(errormsg)

            # Set label and parameters
            self.label = label
            self.p = variant_pars

            return


        def initialize(self, sim):
            ''' Update variant info in sim '''
            self.days = itvs.process_days(sim, self.days) # Convert days into correct format
            sim['variant_pars'][self.label] = self.p  # Store the parameters
            self.index = list(sim['variant_pars'].keys()).index(self.label) # Find where we are in the list
            sim['variant_map'][self.index]  = self.label # Use that to populate the reverse mapping
            self.initialized = True
            return


        def apply(self, sim):
            ''' Introduce new infections with this variant '''
            for ind in itvs.find_day(self.days, sim.t, interv=self, sim=sim): # Time to introduce variant
                susceptible_inds = utls.true(sim.people.susceptible)
                rescale_factor = sim.rescale_vec[sim.t] if self.rescale else 1.0
                scaled_imports = self.n_imports/rescale_factor
                n_imports = sc.randround(scaled_imports) # Round stochastically to the nearest number of imports
                if self.n_imports > 0 and n_imports == 0 and sim['verbose']:
                    msg = f'Warning: {self.n_imports:n} imported infections of {self.label} were specified on day {sim.t}, but given the rescale factor of {rescale_factor:n}, no agents were infected. Increase the number of imports or use more agents.'
                    print(msg)
                importation_inds = np.random.choice(susceptible_inds, n_imports, replace=False) # Can't use utls.choice() since sampling from indices
                sim.people.infect(inds=importation_inds, layer='importation', variant=self.index)
                sim.results['n_imports'][sim.t] += n_imports
            return

