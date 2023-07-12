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
        
    def define_cohort(): 
        pass
        
    def is_within_age_range(number, age_range):
        lower_bound, upper_bound = age_range
        return lower_bound <= number <= upper_bound


    def assign_donation_probs(self, sample_frame_pars): #think about this with people
        if self.study_design_pars['design_type'] == 'longitudinal':
            #call define cohort function 
            pass

        elif self.study_design_pars['design_type'] == 'cross_sectional':
            #find people in each age range of the sample frame 
            for ID in range(self.people.pars['pop_size']):
                #get associated age interval from the age of the agent
                index = 0 #by default where to look for the associated prob
                for i, interval in enumerate(self.sample_frame_pars['age_interval']): 
                    if (is_within_age_range(self.people['age'][ID], interval)):
                        index = i
                        break
                daily_prob_scaling = ((self.sample_frame_pars['donor_breakdown_per_interval'][index] * (self.sample_frame_pars['num_donors_per_day'] 
                                                                                                        / self.people.pars['pop_size'])) / (self.people['age'][self.people['age'][ID]] 
                                                                                                                                            / self.people.pars['pop_size']))
                #fix sim pars part, also make sure that we can get all people of a certain age in this way, how do we get all people in a certain age range?
                age_prob = daily_prob_scaling * (self.people['age'][ID]/self.people['age']) #all people of this age)
                sex_prob = daily_prob_scaling * (self.people['sex'][ID]/self.people['sex']) #all people of this sex
                if age_prob == 0 or sex_prob == 0:
                    self.people['provides_sample_prob'][ID] = 0.0
                else: 
                    self.people['provides_sample_prob'][ID] = round(((age_prob + sex_prob)/2), 2) #should double check how rounding affects 


               

    #np.where(self.people['age'] == self.people['age'][ID])             

                

                





                




        #Function to assign everyona a probability of donating on any given day
        #If some distribution has been passed in then use this or otherwise 
        #analogous to parameters.py in a serperate file that has predefined arrays or dictionary 
         

    def apply(probability_list, num_samples=None):
        #Function that is called every day which returns a list of samples to move onto testing 
        pass

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

