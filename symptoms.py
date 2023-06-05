import numpy as np 
import infection.utils as cvu


class SymptomGenerator_ArrBased: 

    def __init__(self, prev_COVID=None, prev_ILI=None): 
        '''
        Attributes:
            cs_symp          (list) : Types of COVID specific symptoms
            ncs_symp         (list) : Types of non-COVID specific symptoms (not in use)
            prev_c_symp     (dict) : Prevalence of different symptoms (specific or non-specific) within symptomatic COVID-19 
            prev_ILI_symp   (dict) : Prevalence of different symptoms (specific or non-specific) within symptomatic ILI
        '''
        self.cs_symp = ['fever', 'cough', 'sore_throat']

        if prev_COVID is None: 
            self.prev_COVID = {  # Order of symptoms needs to be the same as cs_symp 
                'fever': 0.8,
                'cough': 0.6,
                'sore_throat': 0.1
            }
           
        else: 
            self.prev_COVID = prev_COVID

        if prev_ILI is None: 
            self.prev_ILI = {
                'fever': 0.3,
                'cough': 0.8,
                'sore_throat': 0.7
            }
        else: 
            self.prev_ILI = prev_ILI 

        return


    def assign_symptoms(self, sim): 
        '''
        Main function call of this class. Calculate the specific type of symptoms each infected person in the population has, whether they have COVID or ILI. 
        '''
        if sim.t == 0: 
            # Initialize symptoms tracker 
            setattr(sim.people, 'symptom_types', np.full((len(self.cs_symp), sim.pars['pop_size']), False))  # symptom_types: p x q shaped boolean array, where p is the number of COVID-specific symptoms, and q is the population size.

        # Reset symptoms of people who recovered today, or died
        self.reset_symptoms(sim, 'COVID')
        self.reset_symptoms(sim, 'ILI')

        # Distribute COVID and ILI symptoms 
        self.calculate_symptom_prevalence(sim, 'prev_COVID')
        self.calculate_symptom_prevalence(sim, 'prev_ILI')

        return 


    def calculate_symptom_prevalence(self, sim, identifier): 
        '''
        On the current day, for those that became symptomatic, calculate the symptoms that they have. Assumed that for any batch of people that become symptomatic, prevalence of symptom varieties remains constant.

        Args: 
            sim         (sim): Simulation object
            identifier  (str): Name of symptom prevalence dictionary, corresponding to COVID or ILI
        '''
        cur = sim.t

        if identifier == 'prev_COVID': 
            prev_dict = self.prev_COVID
            sym_cur = cvu.true(sim.people['date_symptomatic'] == cur)  # Those who became symptomatic with COVID today 
        elif identifier == 'prev_ILI': 
            prev_dict = self.prev_ILI   
            sym_cur = cvu.true(sim.people['date_symptomatic_ILI'] == cur)  # Those who became symptomatic with ILI today
        else: 
            raise RuntimeError('Identifier must be either COVID or ILI ')

        for i in range(len(prev_dict.keys())):
            n_symptom = int(len(sym_cur) * prev_dict[self.cs_symp[i]])  # Get symptom prevalence among those who became symptomatic today
            has_symptom = np.random.choice(sym_cur, n_symptom, replace=False)  # Determine who gets a symptom
            
            sim.people.symptom_types[i][has_symptom] = True
                
    
    def reset_symptoms(self, sim, identifier): 
        '''
        On the current day, identify who has recovered or died from COVID/ILI, and reset their symptoms

        Args: 
            sim          (sim): Simulation object 
            identifier   (str): Whether to reset COVID or ILI symptoms
        '''
        if identifier == 'COVID': 
            date_recovered = sim.people.date_recovered
            date_dead = sim.people.date_dead
        elif identifier == 'ILI': 
            date_recovered = sim.people.date_recovered_ILI 
            date_dead = None  # We are not currently modelling ILI deaths
        else: 
            raise RuntimeError('Identifier must be either COVID or ILI ')

        # Remove symptoms for recovered people
        recovered_today = cvu.true(date_recovered == sim.t)
        sim.people.symptom_types[:, recovered_today] = False

        # Remove symptoms for dead people (only applicable for COVID)
        if not(date_dead is None):
            died_today = cvu.true(date_dead == sim.t)
            sim.people.symptom_types[:, died_today] = False


class Symptom_generator: 

    def __init__(self, prev_COVID=None, prev_ILI=None): 
        '''
        Attributes:
            cs_symp          (list) : Types of COVID specific symptoms
            ncs_symp         (list) : Types of non-COVID specific symptoms (not in use)
            prev_c_symp     (dict) : Prevalence of different symptoms (specific or non-specific) within symptomatic COVID-19 
            prev_ILI_symp   (dict) : Prevalence of different symptoms (specific or non-specific) within symptomatic ILI
        '''
        self.cs_symp = ['fever', 'cough', 'sore_throat']
        # self.ncs_symp = ['tight_chest', 'fatigue', 'nausea', 'runny nose']

        if prev_COVID is None: 
            self.prev_COVID = {
                'fever': 0.8,
                'cough': 0.6,
                'sore throat': 0.1
            }
           
        else: 
            self.prev_COVID = prev_COVID

        if prev_ILI is None: 
            self.prev_ILI = {
                'fever': 0.3,
                'cough': 0.8,
                'sore throat': 0.7
            }
        else: 
            self.prev_ILI = prev_ILI 

        return


    def assign_symptoms(self, sim): 
        '''
        Main function call of this class. Calculate the specific type of symptoms each infected person in the population has, whether they have COVID or ILI. 
        '''
        cur = sim.t

        if cur == 0: 
            # Initialize symptoms tracker 
            pop_size = sim.pars['pop_size']
            setattr(sim.people, 'symptom_types', np.array([set() for _ in range(pop_size)]))  # symptom_types: Array of sets, where each set describes the symptoms an individual has

        # Reset symptoms of people who recovered today, or died
        self.reset_symptoms(sim, 'COVID')
        self.reset_symptoms(sim, 'ILI')

        # Distribute COVID and ILI symptoms 
        self.calculate_symptom_prevalence(sim, 'prev_COVID')
        self.calculate_symptom_prevalence(sim, 'prev_ILI')

        # counter = 0
        # for symptom_set in sim.people.symptom_types: 
        #     if len(symptom_set) > 0: 
        #         counter += 1 
        # print('Number of people with symptoms in population', counter)
        # print("All symptoms in population: ", sim.people.symptom_types)

        return 


    def calculate_symptom_prevalence(self, sim, identifier): 
        '''
        On the current day, for those that became symptomatic, calculate the symptoms that they have. Assumed that for any batch of people that become symptomatic, prevalence of symptom varieties remains constant.

        Args: 
            sim         (sim): Simulation object
            identifier  (str): Name of symptom prevalence dictionary, corresponding to COVID or ILI
        '''
        cur = sim.t

        if identifier == 'prev_COVID': 
            prev_dict = self.prev_COVID
            sym_cur = cvu.true(sim.people['date_symptomatic'] == cur)  # Those who became symptomatic with COVID today 
            # print("Became symptomatic due to COVID: ", len(sym_cur))
        elif identifier == 'prev_ILI': 
            prev_dict = self.prev_ILI   
            sym_cur = cvu.true(sim.people['date_symptomatic_ILI'] == cur)  # Those who became symptomatic with ILI today
            # print('Became symptomatic due to ILI: ', len(sym_cur))
        else: 
            raise RuntimeError('Identifier must be either COVID or ILI ')

        for symptom in prev_dict.keys(): 
            n_symptom = int(len(sym_cur) * prev_dict[symptom])  # Get symptom prevalence among those who became symptomatic today
            has_symptom = np.random.choice(len(sym_cur), n_symptom, replace=False)  # Determine who gets a symptom
            
            for i in has_symptom: 
                sim.people.symptom_types[sym_cur[i]].add(symptom)  # Add this to the agent's set of symptoms
                
    
    def reset_symptoms(self, sim, identifier): 
        '''
        On the current day, identify who has recovered or died from COVID/ILI, and reset their symptoms

        Args: 
            sim          (sim): Simulation object 
            identifier   (str): Whether to reset COVID or ILI symptoms
        '''
        if identifier == 'COVID': 
            date_recovered = sim.people.date_recovered
            date_dead = sim.people.date_dead
        elif identifier == 'ILI': 
            date_recovered = sim.people.date_recovered_ILI 
            date_dead = None  # We are not currently modelling ILI deaths
        else: 
            raise RuntimeError('Identifier must be either COVID or ILI ')

        # Reset symptoms for recovered people
        recovered_today = cvu.true(date_recovered == sim.t)
        for recov_ind in recovered_today: 
            sim.people.symptom_types[recov_ind].clear()

        # Reset symptoms for dead people (only applicable for COVID)
        if not(date_dead is None):
            died_today = cvu.true(date_dead == sim.t)
            for died_ind in died_today: 
                sim.people.symptom_types[died_ind].clear()



# NOTE: This entire class below is deprecated

# TODO: Eventually, we can get more granularity on test seeking probabilities given specific symptoms, using something like this (although, hopefully constructed more simply)
# For now, it's sufficient to just identify specific vs nonspecific symptoms
# class Symptom_seeking: 
#     '''
#     Used to determine who seeks a test in two cases: 
#         1) The agent has symptoms
#         2) The agent is in contact with someone who has symptoms
#     '''
#     def __init__(self): 
#         '''
#         Attributes: 
#             cs_symp          (list) : Types of COVID specific symptoms
#             ncs_symp         (list) : Types of non-COVID specific symptoms
#             cs_symp_seek     (dict) : Types of COVID specific symptoms (associated with higher probability of test seeking), and the posterior probability of test seeking
#             ncs_symp_seek    (dict) : Types of non-COVID specific symptoms (associated with lower probability of test seeking), and posterior probability of test seeking
#             prev_c_symp     (dict) : Prevalence of different symptoms (specific or non-specific) within symptomatic COVID-19 
#             prev_ILI_symp   (dict) : Prevalence of different symptoms (specific or non-specific) within symptomatic ILI
#         '''
#         self.cs_symp = set(['fever', 'cough', 'sore_throat', 'short_breath', 'loss_taste_smell', 'muscle_joint_pain', 'chills_shaking'])
#         self.ncs_symp = set(['tight_chest', 'fatigue', 'nausea', 'runny nose'])
#         self.cs_symp_seek = {
#             'fever': 0.8,
#             'cough': 0.8, 
#             'sore_throat': 0.8, 
#             'short_breath': 0.8,
#             'loss_taste_smell': 0.8,
#             'muscle_joint_pain': 0.8, 
#             'chills_shaking': 0.8
#         }

#         self.ncs_symp_seek = {  # NOTE: Some of these symptoms actually are associated with seeking a COVID PCR test
#             'tight_chest': 0.4,
#             'fatigue': 0.4,
#             'nausea': 0.4, 
#             'runny nose': 0.4,
#         }

    
#     def stratify_sx(self, sim, sym_uid): 
#         '''
#         Stratify those who have symptoms by whether their symptoms are COVID or non-COVID specific.

#         Args: 
#             sim      (sim)      : Simulation object
#             sym_uid  (np.array) : UID of symptomatic people
#         '''
#         # Stratify by COVID or non-COVID specific symptoms
#         cs_uid = set() 
#         ncs_uid = set() 

#         for uid in sym_uid: 
#             if len(sim.people.symptom_types[uid].intersection(self.cs_symp)) > 0:  # Has COVID specific symptoms 
#                 cs_uid.add(uid)
#             elif len(sim.people.symptom_types[uid].intersection(self.cs_symp)) == 0:  # Does not have COVID specific symptoms
#                 ncs_uid.add(uid)

#         return cs_uid, ncs_uid


#     def stratify_contacts_sx(self, sim, sym_uid): 
#         '''
#         Get contacts of people who are symptomatic with either COVID-specific or non COVID-specific symptoms.

#         Args: 
#             sim         (sim) : Simulation object
#             sym_uid     (set) : Set of people with symptoms (either due to COVID or ILI)

#         Returns: 
#             contacts    (dict): Dictionary where keys reference sets of people who either contact people with COVID-specific or non COVID-specific symptoms
#         '''
#         cs_uid, ncs_uid = self.stratify_sx(sim, sym_uid)

#         contacts = dict()
#         cs_contacts = set()
#         ncs_contacts = set()

#         for l in sim.people.contacts.keys(): 
#             cs_contacts.update(sim.people.contacts[l].find_contacts(np.array(list(cs_uid))))
#             ncs_contacts.update(sim.people.contacts[l].find_contacts(np.array(list(ncs_uid))))

#         contacts['cs_contacts'] = cs_contacts
#         contacts['ncs_contacts'] = ncs_contacts

#         return contacts

#     def get_cont_sx_seekers(self, sim, sym_uid, probs):
#         '''
#         Among those who are contacting symptomatic individuals, decide who seeks a test, depending on the nature of contact. 
#         Args: 
#             sim         (sim)  : Simulation object
#             sym_uid     (set)  : Set of people with symptoms (either due to COVID or ILI) 
#             probs       (dict) : Probability of seeking a test given they contact someone with COVID-specific or non COVID-specific symptoms 
        
#         Returns: 
#         '''
#         contacts = self.stratify_contacts_sx(sim, sym_uid)
#         seekers = set()

#         for cont_type in contacts.keys(): 

#             prob = probs[cont_type]

#             for c in contacts[cont_type]: 
#                 if np.random.binomial(1, prob): 
#                     seekers.add(c)

#         return seekers 
    
#     def get_sx_seekers(self, sim, sym_uid, probs): 
#         '''
#         Among those who are symptomatic with COVID or ILI, determine who seeks a test, depending on the nature of symptoms. 

#         Args: 
#             sim         (sim)  : Simulation object
#             sym_uid     (set)  : Set of people with symptoms (either due to COVID or ILI) 
#             probs       (dict) : Probability of seeking a test given they have COVID-specific or non COVID-specific symptoms
#         '''
#         cs_uid, ncs_uid = self.stratify_sx(sim, sym_uid)

#         seekers = set()

#         for c in cs_uid: 
#             if np.random.binomial(1, probs['cs_symp']):
#                 seekers.add(c)

#         for nc in ncs_uid: 
#             if np.random.binomial(1, probs['ncs_symp']): 
#                 seekers.add(nc)

#         return seekers
        


# Archived symptom prevalences
 # self.prev_COVID = {  # NOTE: Prevalence numbers are entirely unreliable, only here for the architecture
            #     'fever': 0.8,
            #     'cough': 0.6, 
            #     'sore_throat': 0.1, 
            #     'short_breath': 0.2,
            #     'loss_taste_smell': 0.3,
            #     'muscle_joint_pain': 0.2, 
            #     'chills_shaking': 0.2,
            #     'tight_chest': 0.1,
            #     'fatigue': 0.3,
            #     'nausea': 0.1, 
            #     'runny nose': 0.1,
            # }

# self.prev_ILI = {
            #     'fever': 0.7,
            #     'cough': 0.9, 
            #     'sore_throat': 0.8, 
            #     'short_breath': 0.3,
            #     'loss_taste_smell': 0.2,
            #     'muscle_joint_pain': 0.9, 
            #     'chills_shaking': 0.9,
            #     'tight_chest': 0.1,
            #     'fatigue': 0.3,
            #     'nausea': 0.1, 
            #     'runny nose': 0.1,
            # }