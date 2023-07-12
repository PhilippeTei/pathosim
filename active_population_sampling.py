import numpy as np
from population_sampling_pars import sample_frame_mapping

class active_population_sampling_program:
    def __init__(self, people, tests, study_design_pars, sample_frame, sample_frame_pars): 
        self.people = people 
        self.tests = tests 
        self.Sampler(people, study_design_pars, sample_frame, sample_frame_pars)
        pass


class Sampler: 
    def __init__(self, people, study_design_pars, sample_frame, sample_frame_pars = None):
        #Initializes sampling group by assigning donations to everyone in the population based on design 
        self.people = people #sim.people object
        self.study_design_pars = study_design_pars
       
        design_type = study_design_pars.get('design_type') #make sure a design type is inputted
        if design_type is None:
            raise ValueError("'design_type' parameter is missing in 'study_design_pars'.")
        
        if study_design_pars['design_type'] == 'longitudinal': #need sampling interval and dropout rate
            self.sampling_interval = study_design_pars['sampling_interval']
            self.dropout_rate = study_design_pars['dropout_rate']

            if self.sampling_interval is None:
                raise ValueError("'sampling_interval' parameter is missing for 'longitudinal' design.")
            if self.dropout_rate is None:
                raise ValueError("'dropout_rate' parameter is missing for 'longitudinal' design.")
            
        elif study_design_pars['design_type'] == 'cross_sectional': #need sampling periods and num people per period
            self.sampling_periods = study_design_pars['sampling_periods']
            self.num_people_captured = study_design_pars['num_people_captured']
            
            if self.sampling_periods is None:
                raise ValueError("'sampling_periods' parameter is missing for 'cross-sectional' design.")
            if self.num_people_captured is None:
                raise ValueError("'num_people_captured' parameter is missing for 'cross-sectional' design.")
            if len(self.sampling_periods) != len(self.num_people_captured):
                raise ValueError("Length of 'sampling_periods' should match the length of 'num_people_captured'.")

        self.sample_frame = sample_frame 
        if sample_frame in sample_frame_mapping: #assign the construction parameters based on the sample frame
            self.sample_frame_pars = sample_frame_mapping[sample_frame]

        elif sample_frame_pars != None: 
            self.sample_frame_pars = sample_frame_pars
     
        else:
           raise ValueError("Invalid sample_frame value.") 
        
        self.assign_donation_probs() #assign provides_sample_prob
        
    def define_cohort(): 
        pass
        

    def assign_donation_probs(self): 
        if self.study_design_pars['design_type'] == 'longitudinal':
            #call define cohort function 
            pass

        elif self.study_design_pars['design_type'] == 'cross_sectional':
            pop_size = self.people.pars['pop_size'] #used in later probability calculations
            #find people in each age range of the sample frame 
            age_array = np.round(self.people['age'])
            interval_arrays = [[] for _ in self.sample_frame_pars['age_intervals']]
            for i, interval in enumerate(self.sample_frame_pars['age_intervals']):
                lower_bound, upper_bound = interval
                mask = np.logical_and(lower_bound <= age_array, age_array <= upper_bound)
                indices_of_people_in_interval = mask.nonzero()[0]
                interval_arrays[i].extend(indices_of_people_in_interval)
                print(f"Indices in interval {self.sample_frame_pars['age_intervals'][i]}: {indices_of_people_in_interval}")

            
            '''
            #seperate indices of males and femlaes in population
            not sure if we even need this since our P(sex) will be based on the breakdown 
            sex_array = self.people['sex']
            mask = (sex_array == 0)
            indices_of_female_agents = mask.nonzero()[0]
            indices_of_male_agents = (~mask).nonzero()[0]
            '''

            for i, interval in enumerate(self.sample_frame_pars['age_intervals']):
                daily_prob_scaling = ((self.sample_frame_pars['donor_breakdown_per_interval'][i]*
                                          (self.sample_frame_pars['num_donors_per_day']/pop_size))/ pop_size)
                print(daily_prob_scaling)
                #everyone in the same age group will have the same scaling factor, and the same age probability 
                age_prob = daily_prob_scaling * (len(interval_arrays[i])/pop_size)
                print(age_prob)
                
                for ID in interval_arrays[i]:
                    if self.people['sex'][ID] == 0: 
                        sex_prob = daily_prob_scaling * self.sample_frame_pars['sex_breakdown']['female'] 
                    else: 
                        sex_prob = daily_prob_scaling * self.sample_frame_pars['sex_breakdown']['male']
                    
                    if age_prob == 0 or sex_prob == 0:
                        self.people['provides_sample_prob'][ID] = 0.0
                    else: 
                        self.people['provides_sample_prob'][ID] = ((age_prob + sex_prob)/2) #should double check how rounding affects 

               

    #np.where(self.people['age'] == self.people['age'][ID])             

    def apply(self, num_samples=None):
        #num_samples should come from study_design_pars 
        #not working yet 
        #Function that is called every day which returns a list of samples to move onto testing
        probability_list = self.people['provides_sample_prob']
        indices_to_sample = np.arange(len(self.people)) 
        selected_indices = np.random.choice(indices_to_sample, size=num_samples, p=probability_list, replace=True)
        selected_to_sample = self.people[selected_indices]
        return selected_to_sample



                





                




        #Function to assign everyona a probability of donating on any given day
        #If some distribution has been passed in then use this or otherwise 
        #analogous to parameters.py in a serperate file that has predefined arrays or dictionary 
         


class Immunoassay: 
    def __init__(self, error, LOD, binary_result=False, pos_threshold=None):
        self.error = error
        self.LOD = LOD

    def apply(self, people, people_indices):
        #Function to test people based on their IgG levels
        
        pass

class Results(): 
    #Results class to store test results 

    def __init__(self) -> None:
        pass

