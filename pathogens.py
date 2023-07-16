import numpy as np
import sciris as sc
from . import utils as utls
from . import defaults as df
from . import parameters as ptpar
from . import interventions as itvs

__all__ = ['Pathogen', 'SARS_COV_2']

class Pathogen(sc.prettyobj):    #PATHOGEN BASE CLASS

    def __init__(self, pop_infected, label = "New_Pathogen"): 
        
        self.label = label
        self.pop_infected = pop_infected #Number of initial infections
        self.variants = [] #variants excluding wild type 


        self.make_default_pathogen_pars()


        return

    def initialize(self, sim):
        if len(self.variants) > 0:
            for variant in self.variants:  
                variant.initialize(sim)
                
        self.update_runtime_pars(sim)
        self.make_variant_cross_immunity_matrix() 
        self.convert_prognoses()
         

        return



    #------------ USER CONFIGURE PATHOGEN FUNCTIONS: ADD VARIANTS, CONFIGURE PARAMETERS OF EXISTING PATHOGEN OR FULLY CUSTOM PATHOGEN ------------#
    #region 
    
    def add_custom_variant(self, label = "custom", days = 0, n_imports = 10, rel_beta = 1.0, rel_symp_prob = 1.0, rel_severe_prob = 1.0, rel_crit_prob = 1.0, rel_death_prob = 1.0):
        '''
        Add a custom variant to the pathogen and configure its parameters
        '''
        params = {
            "rel_beta": rel_beta,
            "rel_symp_prob":rel_symp_prob,
            "rel_severe_prob":rel_severe_prob,
            "rel_crit_prob":rel_crit_prob,
            "rel_death_prob":rel_death_prob
            }
        variant = self.Variant(label = label, variant = params, days=days, n_imports=n_imports, base_pathogen = self)
        self.variants.append(variant)  
        return
    
    def add_existing_variant(self, name = "wild", days = 10, n_imports = 10): 
        '''
        Add a variant to the pathogen that is already pre-configured
        '''
        variant  = self.Variant(name, days=days, n_imports=n_imports, base_pathogen = self)
        self.variants.append(variant) 
        return

    def configure_variant_cross_immunity(self, cross_immunity_dict): 
        self.make_variant_cross_immunity_matrix(False) #init and fill self.default_variant_cross_immunity with 1s
 
        for key in cross_immunity_dict.keys():
            for k in cross_immunity_dict[key].keys():
                self.default_variant_cross_immunity[key][k] = cross_immunity_dict[key][k] 

    def configure_transmission(self, beta_dist = None, viral_dist= None, beta= None, asymp_factor= None, viral_levels= None):
        if beta_dist is not None: 
            self.beta_dist    = beta_dist
        if viral_dist is not None:
            self.viral_dist   = viral_dist
        if beta is not None: 
            self.beta         = beta
        if asymp_factor is not None: 
            self.asymp_factor = asymp_factor
        if viral_levels is not None: 
            self.viral_levels = viral_levels

    def configure_disease_state_durations(self, exp2inf= None, inf2sym= None, sym2sev= None, sev2crit= None, asym2rec= None, mild2rec= None, sev2rec= None, crit2rec= None, crit2die= None):
 
        if exp2inf is not None: 
            self.dur['exp2inf']  = exp2inf 
        if inf2sym is not None:            
            self.dur['inf2sym']  = inf2sym 
        if sym2sev is not None:            
            self.dur['sym2sev']  = sym2sev 
        if sev2crit is not None:           
            self.dur['sev2crit'] = sev2crit
        
        if  asym2rec is not None:          
            self.dur['asym2rec'] = asym2rec
        if mild2rec is not None:           
            self.dur['mild2rec'] =mild2rec 
        if sev2rec is not None:            
            self.dur['sev2rec']  = sev2rec 
        if crit2rec is not None:           
            self.dur['crit2rec'] = crit2rec
        if crit2die is not None:           
            self.dur['crit2die'] = crit2die

    def configure_prognoses(self, prognoses, prog_by_age = True): 
        self.prog_by_age      = prog_by_age 
        for key in prognoses.keys():
            self.prognoses[key] = prognoses[key]

    def configure_generalized_immunity(self, min_immunity = 0, max_immunity = 1, duration_to_min_immunity = 180, duration_to_max_immunity = 14):
        self.use_nab_framework = False
        self.imm_peak = max_immunity
        self.imm_final_value = min_immunity
        self.imm_days_to_peak = duration_to_max_immunity
        self.imm_days_to_final_value = duration_to_min_immunity

    def configure_imports(self, imports = 0):
        self.n_imports = imports
    #endregion

    #------------------ DEFAULT PATHOGEN PARAMETRS, ALL ARE TO POSSIBLY OVERRIDE IN DERIVED CLASS ----------------#
    #region
    #TO OVERRIDE
    def make_default_pathogen_pars(self):
        
        self.n_imports  = 0 # Average daily number of imported cases (actual number is drawn from Poisson distribution)
        self.n_variants = 1 # The number of variants of this pathogen, NOT SET BY USER

        #Basic disease transmission parameters
        self.beta_dist    = dict(dist='neg_binomial', par1=1.0, par2=0.45, step=0.01)   # Distribution to draw individual level transmissibility
        self.viral_dist   = dict(frac_time=0.3, load_ratio=2, high_cap=4)               #  The time varying viral load (transmissibility)
        self.beta         = 0.016                                                       # Beta per symptomatic contact; absolute value
        self.asymp_factor = 1.0                                                         # Multiply beta by this factor for asymptomatic cases
        self.viral_levels = dict(min_vl=0.75, max_vl=2)                                 # Specifies the range within which viral load should be scaled so it can contribute to relative transmissibility

         # Variant-specific disease transmission parameters. By default, these are set up for a single variant, but can all be modified for multiple variants
        self.rel_beta        = 1.0 # Relative transmissibility varies by variant
 
          
        # Parameters used to calculate immunity 
        self.use_nab_framework = True

        self.nab_init     = dict(dist='normal', par1=0, par2=2)  # Parameters for the distribution of the initial level of log2(nab) following natural infection
        self.nab_decay    = dict(form='nab_growth_decay', growth_time=21, decay_rate1=np.log(2) / 50, decay_time1=150, decay_rate2=np.log(2) / 250, decay_time2=365)
        self.nab_kin    = None # Constructed during sim initialization using the nab_decay parameters
        self.nab_boost    = 1.5 # Multiplicative factor applied to a person's nab levels if they get reinfected
        self.nab_eff      = dict(alpha_inf=1.08, alpha_inf_diff=1.812, beta_inf=0.967, alpha_symp_inf=-0.739, beta_symp_inf=0.038, alpha_sev_symp=-0.014, beta_sev_symp=0.079) # Parameters to map nabs to efficacy
        self.rel_imm_symp = dict(asymp=0.85, mild=1, severe=1.5) # Relative immunity from natural infection varies by symptoms. Assumption.
        self.trans_redux  = 0.59  # Reduction in transmission for breakthrough infections

        #immunity
        #growth rate after infection (exponential), max immunity,  decay rate
        self.imm_peak = 0.95  
        self.imm_final_value = 0.1
        self.imm_days_to_peak = 14  
        self.imm_days_to_final_value = 180 


        #Cross-variant immunity
        self.immunity     = None  # Matrix of immunity and cross-immunity factors, set later

        # Duration parameters: time for disease progression
        self.dur = {}
        self.dur['exp2inf']  = dict(dist='lognormal_int', par1=4.5, par2=1.5) # Duration from exposed to infectious;
        self.dur['inf2sym']  = dict(dist='lognormal_int', par1=1., par2=1) # Duration from infectious to symptomatic; 
        self.dur['sym2sev']  = dict(dist='lognormal_int', par1=6, par2=4.5) # Duration from symptomatic to severe symptoms;
        self.dur['sev2crit'] = dict(dist='lognormal_int', par1=1, par2=2.0) # Duration from severe symptoms to requiring ICU;
        
        # Duration parameters: time for disease recovery
        self.dur['asym2rec'] = dict(dist='lognormal_int', par1=8.0,  par2=2.0) # Duration for asymptomatic people to recover;  
        self.dur['mild2rec'] = dict(dist='lognormal_int', par1=8.0,  par2=2.0) # Duration for people with mild symptoms to recover;  
        self.dur['sev2rec']  = dict(dist='lognormal_int', par1=20, par2=6) # Duration for people with severe symptoms to recover, 24.7 days total; 
        self.dur['crit2rec'] = dict(dist='lognormal_int', par1=20, par2=6) # Duration for people with critical symptoms to recover; 
        self.dur['crit2die'] = dict(dist='lognormal_int', par1=11, par2=5) # Duration from critical symptoms to death, 18.8 days total;  

         # Severity parameters: probabilities of symptom progression
        self.rel_symp_prob    = 1.0  # Scale factor for proportion of symptomatic cases
        self.rel_severe_prob  = 1.0  # Scale factor for proportion of symptomatic cases that become severe
        self.rel_crit_prob    = 1.0  # Scale factor for proportion of severe cases that become critical
        self.rel_death_prob   = 1.0  # Scale factor for proportion of critical cases that result in death
        
        # ILI: Proportion of population with symptomatic ILI on any given day. Can be a list, in which case length needs to be at least ceil((n_days + 1) / 7).
        self.bkg_ILI = 0.04
         
        self.default_variant_cross_immunity = self.get_existing_variant_cross_immunity()
        
        self.prog_by_age      = True # Whether to set disease progression based on the person's age
        self.prognoses        = self.get_default_prognoses(self.prog_by_age) # prognoses
        
 
    #TO OVERRIDE
    def get_default_prognoses(self, by_age=True):
        '''
        Return the default parameter values for prognoses 
        The prognosis probabilities are conditional given the previous disease state. 
        ''' 

        #Generic prognoses
        prognoses = dict(
            age_cutoffs   = np.array([0]),
            sus_ORs       = np.array([1.00]),
            trans_ORs     = np.array([1.00]),
            symp_probs    = np.array([1.00]),
            comorbidities = np.array([1.00]),
            severe_probs  = np.array([0.10]),
            crit_probs    = np.array([0.05]),
            death_probs   = np.array([0.01]),
        )                                                                                              
        return prognoses
       
    #TO OVERRIDE
    def get_existing_variants(self):
        ''' Define the parameters of pre-existing variants, made to be overriden by classes that inherit from class Pathogen ''' 
        return dict()

    #TO OVERRIDE
    def get_existing_variant_cross_immunity(self):
        '''Return existing variant cross_immunity, made to be overriden '''
        pars = dict(wild = dict(wild  = 1.0))
        return pars

    #endregion 


    #------------------------------ OTHER PATHOGEN AND VARIANT FUNCTIONS ---------------------------------------------#
    #region 
         
    def get_variants_labels(self):
        '''
        Returns a list of the labels ("names") of all the variants including the wild type
        '''

        labels = ['wild']
        for v in self.variants:
            labels.append(v.label)

        return labels

    def make_variant_cross_immunity_matrix(self, init_arrays = True):
        '''Adds any missing cross-variant immunity values in the dictionariy self.variant_cross_immunity, 
        assuming complete cross-immunity between variants as a default (default value of 1)'''

        variant_labels = ['wild']
        for v in self.variants:
            variant_labels.append(v.label)

        #complete the dict of variant cross-immunity with missing values for custom variants
        for variant in variant_labels:   
            if not variant in self.default_variant_cross_immunity:
                self.default_variant_cross_immunity[variant] = {'wild':1}
                for v in variant_labels: 
                        self.default_variant_cross_immunity[variant][v] = 1 #Assume complete cross-variant immunity 
            else:
                for v in variant_labels:
                    if v not in self.default_variant_cross_immunity[variant]:
                        self.default_variant_cross_immunity[variant][v] = 1 #Assume complete cross-variant immunity
         
        if not init_arrays:
            return

        # Convert the dict into a np 2D array
        nv = self.n_variants
        self.immunity = np.ones((nv, nv), dtype=df.default_float)  # Fill with defaults
 
        default_cross_immunity =self.default_variant_cross_immunity
        for i in range(nv):
            label_i = self.get_variants_labels()[i]
            for j in range(nv):
                label_j = self.get_variants_labels()[j]
                if label_i in default_cross_immunity and label_j in default_cross_immunity:
                    self.immunity[j][i] = default_cross_immunity[label_j][label_i]
           
    def update_runtime_pars(self, sim):
        '''updates any parameters that should be set automatically in the beginning of a simulation'''
        self.n_variants = len(self.variants)+1
        self.pathogen_index = self.get_pathogen_index(sim)

        if self.imm_final_value <= 0.001:
            self.imm_final_value = 0.005

    def convert_prognoses(self):  
        out = sc.dcp(self.prognoses)
        out['death_probs']  /= out['crit_probs']   # Conditional probability of dying, given critical symptoms
        out['crit_probs']   /= out['severe_probs'] # Conditional probability of symptoms becoming critical, given severe
        out['severe_probs'] /= out['symp_probs']   # Conditional probability of symptoms becoming severe, given symptomatic

        self.prognoses = out

    def get_pathogen_index(self, sim):
        for p in range(len(sim.pathogens)):
            if sim.pathogens[p] == self:
                return p
                

    #endregion
        
    class Variant(sc.prettyobj):
        '''
        Add a new variant to the sim, this should not be created directly by the user!

        Args:
            variant (str/dict): name of variant, or dictionary of parameters specifying information about the variant
            days   (int/list): day(s) on which new variant is introduced
            label       (str): if variant is supplied as a dict, the name of the variant
            n_imports   (int): the number of imports of the variant to be added 
            base_pathogen (Pathogen): the pathogen from which the variant is from
        '''

        def __init__(self, variant, days, label=None, n_imports=1, base_pathogen = None):
            self.days = days # Handle inputs
            self.n_imports = int(n_imports)
            self.index     = None # Index of the variant in the sim; set later
            self.label     = None # Variant label (used as a dict key)
            self.p         = None # This is where the parameters will be stored
            self.pathogen_index = 0
            self.initialized = False
            self.base_pathogen = base_pathogen

            
            self.parse(variant=variant, label=label) # Variants can be defined in different ways: process these here
            return


        def parse(self, variant=None, label=None):
            ''' Unpack variant information, given by the pathogen class''' 
            if isinstance(variant, str):
                  
                known_variant_pars = self.base_pathogen.get_existing_variants()  

                label = variant.lower()

                if label in known_variant_pars.keys(): 
                    variant_pars = known_variant_pars[label]
                else: 
                    raise NotImplementedError(f'The selected variant "{variant}" is not implemented; choices are:\n{sc.pp(known_variant_pars.keys(), doprint=False)}')

            elif isinstance(variant, dict):
                 
                # Parse label
                variant_pars = variant
                label = variant_pars.pop('label', label) # Allow including the label in the parameters
                if label is None:
                    label = 'custom' 
            else:
                errormsg = f'Could not understand {type(variant)}, please specify as a dict or a predefined variant:\n{sc.pp(known_variant_pars.keys(), doprint=False)}'
                raise ValueError(errormsg)
              
            self.label = label
            self.p = variant_pars 
            return


        def initialize(self, sim):
            ''' Update variant info in sim '''
            self.days = itvs.process_days(sim, self.days) # Convert days into correct format
              
            #Set variant index in variants list in Pathogen(0 = wild type)
            for v in range(len(self.base_pathogen.variants)):
                if self.base_pathogen.variants[v].label == self.label:
                    self.index = v+1 
                    break 

            #Set variant's pathogen index, in the list of pathogens in the Sim class
            for p in range(len(sim.pathogens)):
                if sim.pathogens[p] == self.base_pathogen:
                    self.pathogen_index = p
                    break

            self.initialized = True
            return


        def apply(self, sim):
            ''' Introduce new infections with this variant '''
            for ind in itvs.find_day(self.days, sim.t, interv=self, sim=sim): # Time to introduce variant
                susceptible_inds = utls.true(sim.people.p_susceptible[self.pathogen_index])
                scaled_imports = self.n_imports
                n_imports = sc.randround(scaled_imports) # Round stochastically to the nearest number of imports
                if self.n_imports > 0 and n_imports == 0 and sim['verbose']:
                    msg = f'Warning: {self.n_imports:n} imported infections of {self.label} were specified on day {sim.t}. Increase the number of imports or use more agents.'
                    print(msg)
                importation_inds = np.random.choice(susceptible_inds, n_imports, replace=False) # Can't use utls.choice() since sampling from indices
                sim.people.infect(inds=importation_inds, layer='importation', variant=self.index, pathogen_index = self.pathogen_index)
                sim.results[self.pathogen_index]['n_imports'][sim.t] += n_imports
            return



#-------------------------------------------  PRE-DEFINED PATHOGENS ----------------------------------------------# 

class SARS_COV_2(Pathogen):
    
    def make_default_pathogen_pars(self):
          
        super().make_default_pathogen_pars() #Only need to override the parameters you want to

        if self.label == 'New_Pathogen':
            self.label = 'SARS-CoV-2'

        #Basic disease transmission parameters
        self.beta_dist    = dict(dist='neg_binomial', par1=1.0, par2=0.45, step=0.01) # Distribution to draw individual level transmissibility; dispersion from https://www.researchsquare.com/article/rs-29548/v1
        self.viral_dist   = dict(frac_time=0.3, load_ratio=2, high_cap=4) # The time varying viral load (transmissibility); estimated from Lescure 2020, Lancet, https://doi.org/10.1016/S1473-3099(20)30200-0
        self.beta         = 0.016  # Beta per symptomatic contact; absolute value, calibrated
        self.asymp_factor = 1.0  # Multiply beta by this factor for asymptomatic cases; no statistically significant difference in transmissibility: https://www.sciencedirect.com/science/article/pii/S1201971220302502
        self.viral_levels = dict(min_vl=0.75, max_vl=2) # Specifies the range within which viral load should be scaled so it can contribute to relative transmissibility
         
        # Parameters used to calculate immunity 
        self.use_nab_framework = True
        self.nab_init     = dict(dist='normal', par1=0, par2=2)  # Parameters for the distribution of the initial level of log2(nab) following natural infection, taken from fig1b of https://doi.org/10.1101/2021.03.09.21252641
        self.nab_decay    = dict(form='nab_growth_decay', growth_time=21, decay_rate1=np.log(2) / 50, decay_time1=150, decay_rate2=np.log(2) / 250, decay_time2=365)
        self.nab_kin    = None # Constructed during sim initialization using the nab_decay parameters
        self.nab_boost    = 1.5 # Multiplicative factor applied to a person's nab levels if they get reinfected. No data on this, assumption.
        self.nab_eff      = dict(alpha_inf=1.08, alpha_inf_diff=1.812, beta_inf=0.967, alpha_symp_inf=-0.739, beta_symp_inf=0.038, alpha_sev_symp=-0.014, beta_sev_symp=0.079) # Parameters to map nabs to efficacy
        self.rel_imm_symp = dict(asymp=0.85, mild=1, severe=1.5) # Relative immunity from natural infection varies by symptoms. Assumption.
        self.immunity     = None   
        self.trans_redux  = 0.59  # Reduction in transmission for breakthrough infections, https://www.medrxiv.org/content/10.1101/2021.07.13.21260393v

           
        self.dur = {}
        self.dur['exp2inf']  = dict(dist='lognormal_int', par1=4.5, par2=1.5) # Duration from exposed to infectious; see Lauer et al., https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7081172/, appendix table S2, subtracting inf2sym duration
        self.dur['inf2sym']  = dict(dist='lognormal_int', par1=1.1, par2=0.9) # Duration from infectious to symptomatic; see Linton et al., https://doi.org/10.3390/jcm9020538, from Table 2, 5.6 day incubation period - 4.5 day exp2inf from Lauer et al.
        self.dur['sym2sev']  = dict(dist='lognormal_int', par1=6.6, par2=4.9) # Duration from symptomatic to severe symptoms; see Linton et al., https://doi.org/10.3390/jcm9020538, from Table 2, 6.6 day onset to hospital admission (deceased); see also Wang et al., https://jamanetwork.com/journals/jama/fullarticle/2761044, 7 days (Table 1)
        self.dur['sev2crit'] = dict(dist='lognormal_int', par1=1.5, par2=2.0) # Duration from severe symptoms to requiring ICU; average of 1.9 and 1.0; see Chen et al., https://www.sciencedirect.com/science/article/pii/S0163445320301195, 8.5 days total - 6.6 days sym2sev = 1.9 days; see also Wang et al., https://jamanetwork.com/journals/jama/fullarticle/2761044, Table 3, 1 day, IQR 0-3 days; std=2.0 is an estimate
        self.dur['asym2rec'] = dict(dist='lognormal_int', par1=8.0,  par2=2.0) # Duration for asymptomatic people to recover; see Wölfel et al., https://www.nature.com/articles/s41586-020-2196-x
        self.dur['mild2rec'] = dict(dist='lognormal_int', par1=8.0,  par2=2.0) # Duration for people with mild symptoms to recover; see Wölfel et al., https://www.nature.com/articles/s41586-020-2196-x
        self.dur['sev2rec']  = dict(dist='lognormal_int', par1=18.1, par2=6.3) # Duration for people with severe symptoms to recover, 24.7 days total; see Verity et al., https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext; 18.1 days = 24.7 onset-to-recovery - 6.6 sym2sev; 6.3 = 0.35 coefficient of variation * 18.1; see also https://doi.org/10.1017/S0950268820001259 (22 days) and https://doi.org/10.3390/ijerph17207560 (3-10 days)
        self.dur['crit2rec'] = dict(dist='lognormal_int', par1=18.1, par2=6.3) # Duration for people with critical symptoms to recover; as above (Verity et al.)
        self.dur['crit2die'] = dict(dist='lognormal_int', par1=10.7, par2=4.8) # Duration from critical symptoms to death, 18.8 days total; see Verity et al., https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext; 10.7 = 18.8 onset-to-death - 6.6 sym2sev - 1.5 sev2crit; 4.8 = 0.45 coefficient of variation * 10.7
         
        self.rel_symp_prob    = 1.0  # Scale factor for proportion of symptomatic cases
        self.rel_severe_prob  = 1.0  # Scale factor for proportion of symptomatic cases that become severe
        self.rel_crit_prob    = 1.0  # Scale factor for proportion of severe cases that become critical
        self.rel_death_prob   = 1.0  # Scale factor for proportion of critical cases that result in death
         
        self.bkg_ILI = 0.04 

          

    def get_default_prognoses(self, by_age=True):
        '''
        Return the default parameter values for prognoses

        The prognosis probabilities are conditional given the previous disease state.

        Args:
            by_age (bool): whether to use age-specific values (default true)

        Returns:
            prog_pars (dict): the dictionary of prognosis probabilities
        '''

        if not by_age: # All rough estimates -- almost always, prognoses by age (below) are used instead
            prognoses = dict(
                age_cutoffs   = np.array([0]),
                sus_ORs       = np.array([1.00]),
                trans_ORs     = np.array([1.00]),
                symp_probs    = np.array([0.75]),
                comorbidities = np.array([1.00]),
                severe_probs  = np.array([0.10]),
                crit_probs    = np.array([0.04]),
                death_probs   = np.array([0.01]),
            )
        else:
            prognoses = dict(
                age_cutoffs   = np.array([0,       10,      20,      30,      40,      50,      60,      70,      80,      90,]),     # Age cutoffs (lower limits)
                sus_ORs       = np.array([0.34,    0.67,    1.00,    1.00,    1.00,    1.00,    1.24,    1.47,    1.47,    1.47]),    # Odds ratios for relative susceptibility -- from Zhang et al., https://science.sciencemag.org/content/early/2020/05/04/science.abb8001; 10-20 and 60-70 bins are the average across the ORs
                trans_ORs     = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Odds ratios for relative transmissibility -- no evidence of differences
                comorbidities = np.array([1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00,    1.00]),    # Comorbidities by age -- set to 1 by default since already included in disease progression rates
                symp_probs    = np.array([0.50,    0.55,    0.60,    0.65,    0.70,    0.75,    0.80,    0.85,    0.90,    0.90]),    # Overall probability of developing symptoms (based on https://www.medrxiv.org/content/10.1101/2020.03.24.20043018v1.full.pdf, scaled for overall symptomaticity)
                severe_probs  = np.array([0.00050, 0.00165, 0.00720, 0.02080, 0.03430, 0.07650, 0.13280, 0.20655, 0.24570, 0.24570]), # Overall probability of developing severe symptoms (derived from Table 1 of https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-16-COVID19-Report-9.pdf)
                crit_probs    = np.array([0.00003, 0.00008, 0.00036, 0.00104, 0.00216, 0.00933, 0.03639, 0.08923, 0.17420, 0.17420]), # Overall probability of developing critical symptoms (derived from Table 1 of https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-16-COVID19-Report-9.pdf)
                death_probs   = np.array([0.00002, 0.00002, 0.00010, 0.00032, 0.00098, 0.00265, 0.00766, 0.02439, 0.08292, 0.16190]), # Overall probability of dying -- from O'Driscoll et al., https://www.nature.com/articles/s41586-020-2918-0; last data point from Brazeau et al., https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-34-ifr/
            ) 
        return prognoses
        
    def get_existing_variants(self):
        '''
        Define the parameters of pre-existing variants, made to be overriden by classes that inherit from class Pathogen
        '''
        pars = dict(
 
            alpha = dict(
                rel_beta        = 1.67, # Midpoint of the range reported in https://science.sciencemag.org/content/372/6538/eabg3055
                rel_symp_prob   = 1.0,  # Inconclusive evidence on the likelihood of symptom development. See https://www.thelancet.com/journals/lanpub/article/PIIS2468-2667(21)00055-4/fulltext
                rel_severe_prob = 1.64, # From https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3792894, and consistent with https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2021.26.16.2100348 and https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/961042/S1095_NERVTAG_update_note_on_B.1.1.7_severity_20210211.pdf
                rel_crit_prob   = 1.0,  # Various studies have found increased mortality for B117 (summary here: https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(21)00201-2/fulltext#tbl1), but not necessarily when conditioned on having developed severe disease
                rel_death_prob  = 1.0,  # See comment above
            ),

            beta = dict(
                rel_beta        = 1.0, # No increase in transmissibility; B1351's fitness advantage comes from the reduction in neutralisation
                rel_symp_prob   = 1.0,
                rel_severe_prob = 3.6, # From https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2021.26.16.2100348
                rel_crit_prob   = 1.0,
                rel_death_prob  = 1.0,
            ),

            gamma = dict(
                rel_beta        = 2.05, # Estimated to be 1.7–2.4-fold more transmissible than wild-type: https://science.sciencemag.org/content/early/2021/04/13/science.abh2644
                rel_symp_prob   = 1.0,
                rel_severe_prob = 2.6, # From https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2021.26.16.2100348
                rel_crit_prob   = 1.0,
                rel_death_prob  = 1.0,
            ),

            delta = dict(
                rel_beta        = 2.2, # Estimated to be 1.25-1.6-fold more transmissible than B117: https://www.researchsquare.com/article/rs-637724/v1
                rel_symp_prob   = 1.0,
                rel_severe_prob = 3.2, # 2x more transmissible than alpha from https://mobile.twitter.com/dgurdasani1/status/1403293582279294983
                rel_crit_prob   = 1.0,
                rel_death_prob  = 1.0,
            ),

            omicron = dict(
                rel_beta        = 2.5, # Made-up parameters. (Learned in calibration)
                rel_symp_prob   = 1.0,
                rel_severe_prob = 0.5,
                rel_crit_prob   = 1.0,
                rel_death_prob  = 1.0,
            )
        ) 
        return pars
     
    def get_existing_variant_cross_immunity(self):
        '''Return existing variant cross_immunity, made to be overriden'''
        pars = dict( 
            wild = dict(
                wild  = 1.0, # Default for own-immunity
                alpha = 0.5, # Assumption
                beta  = 0.5, # https://www.nature.com/articles/s41586-021-03471-w
                gamma = 0.34, # Assumption
                delta = 0.374, # Assumption
                omicron = 15/315 # https://www.medrxiv.org/content/10.1101/2021.12.24.21268317v6
            ),

            alpha = dict(
                wild  = 0.5, # Assumption
                alpha = 1.0, # Default for own-immunity
                beta  = 0.8, # Assumption
                gamma = 0.8, # Assumption
                delta = 0.689,  # Assumption
                omicron = 0.5*6/301 # https://www.medrxiv.org/content/10.1101/2021.12.24.21268317v6
            ),

            beta = dict(
                wild  = 0.066, # https://www.nature.com/articles/s41586-021-03471-w
                alpha = 0.5,   # Assumption
                beta  = 1.0,   # Default for own-immunity
                gamma = 0.5,   # Assumption
                delta = 0.086,    # Assumption
                omicron = 0.066*8/91 # https://www.medrxiv.org/content/10.1101/2021.12.24.21268317v6
            ),

            gamma = dict(
                wild  = 0.34, # Previous (non-P.1) infection provides 54–79% of the protection against infection with P.1 that it provides against non-P.1 lineages: https://science.sciencemag.org/content/early/2021/04/13/science.abh2644
                alpha = 0.4,  # Assumption based on the above
                beta  = 0.4,  # Assumption based on the above
                gamma = 1.0,  # Default for own-immunity
                delta = 0.088,   # Assumption
                omicron = 0.05 # ASSUMPTION
            ),

            delta = dict( # Parameters from https://www.cell.com/cell/fulltext/S0092-8674(21)00755-8
                wild  = 0.374,
                alpha = 0.689,
                beta  = 0.086,
                gamma = 0.088,
                delta = 1.0, # Default for own-immunity
                omicron = 0.374*42/464 # https://www.medrxiv.org/content/10.1101/2021.12.24.21268317v6
            ),

            omicron = dict(
                wild  = 0.2, # ASSUMPTION
                alpha = 0.2, # ASSUMPTION
                beta  = 0.2, # ASSUMPTION
                gamma = 0.2, # ASSUMPTION
                delta = 28/46, # https://www.nature.com/articles/s41586-022-04830-x
                omicron = 1.0 # Default for own-immunity
            ),
        )
        return pars
