import numpy as np
import matplotlib.pyplot as plt
from population_sampling_pars import sample_frame_mapping


class active_population_sampling_program:
    '''
    Class not fully functioning yet. Will combine the Sampler, Test and Results class. 
    '''
    def __init__(self, people, test, study_design_pars, sample_frame, sample_frame_pars = None): 
        self.people = people 
        self.test = test
        self.sampler = Sampler(people, study_design_pars, sample_frame, sample_frame_pars)
        self.results = Results()

        if study_design_pars['design_type'] == 'longitudinal':
            self.time = self.sampler.sampling_interval
            #at end of every self.sampling interval (day #) in the simulation: 
                #apply a test
                #send results somewhere as (day #, results of all people tested)
                #could directly call plot test results, make modifications to only look at interval tested 
        
        elif study_design_pars['design_type'] == 'cross_sectional':
            self.time = self.sampler.sampling_periods
        

                #should i set these people who have now been sampled to have a probability of 0? or use their resample factor somehow


        #for loop for every test period, with each number of samples at each iteration 
     

        #based on the created active surveillance program, the simulation object will sample a certain number of people per day
        #or will continue to take "samples" from the same peopl, to be implemented later 

    def apply(self, num_people=None):
        #at end of every self.sampling interval (day #) in the simulation: 
                #apply a test
                #send results somewhere as (day #, results of all people tested)
                #could directly call plot test results, make modifications to only look at interval tested 
        if self.sampler.study_design_pars['design_type'] == 'longitudinal':
            people_to_test = self.sampler.apply(self.sampler.num_participants)
        elif self.sampler.study_design_pars['design_type'] == 'cross_sectional':
            people_to_test = self.sampler.apply(num_people)
        
        test_results = self.test.apply(self.people, people_to_test)
        return people_to_test, test_results
    
    #should i set these people who have now been sampled to have a probability of 0? or use their resample factor somehow
class Sampler: 
    '''
    A class which uses a inputted sample frame and design to assign a probability to each person in the population of being sampled. 

    Args: 
        people (People): the population to draw samples from
        study_design_pars (dict): whether it is a longitudinal or cross-sectional study and the associated parameters for these designs
        sample_frame (string): to map to an existing sample frame with age and sex breakdown parameters, or, a name for a new sample frame
        sample_frame_pars (dict): a dictionary with sample frame age and sex breakdown pars

    Example: 
    sim = inf.Sim(pop_size=100, pop_type ='random', pathogens=[pathogen])
    sim.run()
    sampler_test = Sampler(sim.people, {'design_type' : 'cross_sectional', 'sampling_periods' : [(0, 14), (15, 30)], 
                                        'num_people_captured': [50, 50]}, 'canadian_blood_donors')
    '''
    def __init__(self, people, study_design_pars, sample_frame, sample_frame_pars = None):
        #Initializes sampling group by assigning donations to everyone in the population based on design 
        self.people = people #sim.people object
        self.study_design_pars = study_design_pars
       
        design_type = study_design_pars.get('design_type') #make sure a design type is inputted
        if design_type is None:
            raise ValueError("'design_type' parameter is missing in 'study_design_pars'.")
        
        if study_design_pars['design_type'] == 'longitudinal': #need sampling interval and dropout rate
            self.sampling_interval = study_design_pars['sampling_interval'] #shoud this be kept in a dictionary rather than two seperate attributes 
            self.dropout_rate = study_design_pars['dropout_rate']
            self.num_participants = study_design_pars['num_participants'] #initial number of participants to follow 

            if self.sampling_interval is None:
                raise ValueError("'sampling_interval' parameter is missing for 'longitudinal' design.")
            if self.dropout_rate is None:
                raise ValueError("'dropout_rate' parameter is missing for 'longitudinal' design.")
            if self.num_participants <= 0: 
                raise ValueError("'num_participants' parameter is invalid for 'longitudinal' design")

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
        
    def define_cohort(self, people): 
        '''
        Function to define a cohort as it's own sample frame to follow longitudinally. When called in assign_donation_probs it will 
        assign all agents not in the cohort to have a probability of 0 and those in the cohort to a probability fluctuating with a drop 
        out rate.

        Could potentially take the same approach of using a sample frame (CanPath), or could choose a random cohort 
        '''
        #find people in each age range of the sample frame 
        age_array = np.round(self.people['age'])
        interval_arrays = [[] for _ in self.sample_frame_pars['age_intervals']]
        for i, interval in enumerate(self.sample_frame_pars['age_intervals']):
            lower_bound, upper_bound = interval
            mask = np.logical_and(lower_bound <= age_array, age_array <= upper_bound)
            indices_of_people_in_interval = mask.nonzero()[0]
            sex_array = self.people[indices_of_people_in_interval]['sex']
            mask = (sex_array == 0)
            indices_of_female_agents = mask.nonzero()[0]
            indices_of_male_agents = (~mask).nonzero()[0]
            #Find how many people you need to fill up in that age and sex bucket 
            


            interval_arrays[i].extend(indices_of_people_in_interval)
        
            #FOR TESTING: 
            #print(f"Indices in interval {self.sample_frame_pars['age_intervals'][i]}: {indices_of_people_in_interval}")



        
    

    def assign_donation_probs(self): 
        '''
        Assigns a probability to each agent in the simulaiton indicating how likely they are to provide a sample given the sample frame's
        age and sex breakdown. 

        '''
        if self.study_design_pars['design_type'] == 'longitudinal':
            #call define cohort function 
            pass

        elif self.study_design_pars['design_type'] == 'cross_sectional':
            pop_size = self.people.pars['pop_size'] 

            #seperate people by sex to use later in calculating p_sex
            sex_array = self.people['sex']
            mask = (sex_array == 0)
            indices_of_female_agents = mask.nonzero()[0]
            indices_of_male_agents = (~mask).nonzero()[0]

            #find people in each age range of the sample frame 
            age_array = np.round(self.people['age'])
            interval_arrays = [[] for _ in self.sample_frame_pars['age_intervals']]
            for i, interval in enumerate(self.sample_frame_pars['age_intervals']):
                lower_bound, upper_bound = interval
                mask = np.logical_and(lower_bound <= age_array, age_array <= upper_bound)
                indices_of_people_in_interval = mask.nonzero()[0]
                interval_arrays[i].extend(indices_of_people_in_interval)
                #FOR TESTING: 
                #print(f"Indices in interval {self.sample_frame_pars['age_intervals'][i]}: {indices_of_people_in_interval}")
            
            
            for i, interval in enumerate(self.sample_frame_pars['age_intervals']):
                
                p_donation = self.sample_frame_pars['num_donors_per_day']/pop_size
                #everyone in the same age group will have the same scaling factor, and the same age probability 
                p_age = len(interval_arrays[i])/pop_size
                
                p_age_given_donation = self.sample_frame_pars['donor_breakdown_per_interval'][i]
                
                #final age factor for provides_donation_prob
                p_donation_given_age = p_age_given_donation * p_donation / p_age
                
                #age_prob = daily_prob_scaling * (len(interval_arrays[i])/pop_size)
                
                #calculate sex probability, average sex and age probability to assign 'provides_sample_prob'
                for ID in interval_arrays[i]:

                    if self.people['sex'][ID] == 0: 
                        p_sex = len(indices_of_female_agents)/pop_size
                        p_sex_given_donation = self.sample_frame_pars['sex_breakdown']['female'] 
                        
                        p_donation_given_sex = p_sex_given_donation * p_donation / p_sex

                        #final probability assignment for females in this age group:
                        if p_donation_given_age == 0 or p_donation_given_sex == 0:
                            self.people['provides_sample_prob'][ID] = 0.0
                        else: 
                            self.people['provides_sample_prob'][ID] = ((p_donation_given_sex + p_donation_given_sex)/2) 

                        
                    else: 
                        p_sex = len(indices_of_male_agents)/pop_size
                        p_sex_given_donation = self.sample_frame_pars['sex_breakdown']['male'] 
                        
                        p_donation_given_sex = p_sex_given_donation * p_donation / p_sex

                        #final probability assignment for males in this age group:
                        if p_donation_given_age == 0 or p_donation_given_sex == 0:
                            self.people['provides_sample_prob'][ID] = 0.0
                        else: 
                            self.people['provides_sample_prob'][ID] = ((p_donation_given_sex + p_donation_given_age)/2) 

               
           

    def apply(self, num_samples):
        '''
        Function that is called every day which returns a list of samples to move onto testing. Returns num_samples of people indices
        in the population based on provides_sample_prob attribute, bernoulli trials, and random selection. 

        Args: 
            num_samples (int): The number of people to select from the population to move onto testing. 

        Returns: 
            An array of people indices of length num_samples. 
        '''
        
        probability_list = self.people['provides_sample_prob']
        #normalized_probs = probability_list / np.sum(probability_list)
        #print(probability_list)
        can_sample = bernoulli_trials(probability_list)
        #print(can_sample)
        # Get the indices of True values in the can_sample list
        true_indices = np.where(can_sample)[0]

        # If the number of True values is less than num_samples, return all of them
        if len(true_indices) <= num_samples:
            return true_indices

        # Randomly select num_samples True values from the list of indices
        selected_indices = np.random.choice(true_indices, num_samples, replace=False)

        return selected_indices
       


class Immunoassay: 
    '''
    A class to create immunoassay testing objects for IgG levels of people. This class can simulate multiple different 
    types of tests such as quantitative and qualitative and allows for different error margins. 

    Args: 
        LOD (float): limit of detection, lower boundary in order to be able to obtain results from a test. 
        error (tuple): range for error, converted to standard deviation used in giving test results. 
        binary_result (bool): indicating whether to run a quantitative or binary test
        pos_threshold (float): value used for binary tests to give a positive or negative result

    **Examples**
    binary_immunoassay_test_false = Immunoassay(3, binary_result=True, pos_threshold=7)
    quantitative_immunoassay_test = Immunoassay(2)

    '''
    def __init__(self, LOD, error=(0.01, 0.4), binary_result=False, pos_threshold=None):
        self.error = error #immunoassays have been reported to have an error within this range https://journals.sagepub.com/doi/full/10.1258/acb.2011.011073 
        self.LOD = LOD
        self.binary_result = binary_result
        self.pos_threshold = pos_threshold


    def apply(self, people, people_indices):
        '''
        Initial form of apply function, need to combine with apply_direct. Applies a test to the people indicated by people_indices. 

        Args: 
            people (People): the population where the people indices come from. 
            people_indices (array): An array of indices in the population to apply the test onto. 

        Returns: An array of test results the same length as people_indices. 
        '''
        IgG_levels = people['IgG_level'][people_indices] #from the Sampler, which people will be tested
        std_deviation = np.random.uniform(self.error[0], self.error[1])
        #assign measured IgG values as a normal distribution with a std_deviation selected from the error range
        measured_IgG = np.random.normal(IgG_levels, std_deviation) #an array of test measurements 
        #print(measured_IgG) FOR TESTING 
        mask = measured_IgG > self.LOD
        if self.binary_result: 
            test_result = np.logical_and(mask, measured_IgG > self.pos_threshold)
        else: 
            test_result = np.where(mask, measured_IgG, False)
        return test_result

    def apply_direct(self, IgG_levels, people_indices):
        '''
        Second version of the apply test function. Tests people at the people indices given and returns the result. 
        Useful for plotting the test results of the people selected over the course of the simulation. 
        
        Args: 
            IgG_levels (array): IgG levels of everyone in the population (from a given day in the simulation) 
            people_indices (array): An array of the people to test. 

        Returns: An array of test results the same length as people_indices. 

        '''
        test_input = IgG_levels[people_indices]
        std_deviation = np.random.uniform(self.error[0], self.error[1])
        #assign measured IgG values as a normal distribution with a std_deviation selected from the error range
        measured_IgG = np.random.normal(IgG_levels, std_deviation) #an array of test measurements 
        #print(measured_IgG) FOR TESTING 
        mask = measured_IgG > self.LOD
        if self.binary_result: 
            test_result = np.logical_and(mask, measured_IgG > self.pos_threshold)
        else: 
            test_result = np.where(mask, measured_IgG, False)
        return test_result


        

class Results(): 
    #Results class to store test results 
    def __init__(self):
        self.test_results = {}

    def store_results(self, day, people_indices, test_results):
        """
        Store the test results for multiple people on a given day.

        Args:
            day (int): The day of the simulation.
            people_indices (list): A list of person indices.
            test_results (list): A list of test results corresponding to the people indices.
        """
        for person_index, test_result in zip(people_indices, test_results):
            if person_index not in self.test_results:
                self.test_results[person_index] = []
            self.test_results[person_index].append((day, test_result))

    def get_results(self, person_index):
        """
        Get the stored test results for a specific person.

        Args:
            person_index (int): The index of the person.

        Returns:
            list: A list of tuples representing test results for the person, where each tuple is (day, test_result).
        """
        return self.test_results.get(person_index, [])
    
    def plot_individual_result(self, person_index):
        """
        Plot the test results for a specific person.

        Args:
            person_index (int): The index of the person.
        """
        results = self.get_results(person_index)
        
        if not results:
            print("No results found for the specified person.")
            return

        days = [day for day, _ in results]
        test_results = [test_result for _, test_result in results]
        
        plt.figure(figsize=(10, 6))
        plt.plot(days, test_results, marker='o')
        plt.title(f"Test Results for Person {person_index}")
        plt.xlabel("Day")
        plt.ylabel("Test Result")
        plt.grid(True)
        plt.show()

    def plot_average_results_each_day(self):
        """
        Plot the average test results for each day.
        """
        if not self.test_results:
            print("No test results found.")
            return
        
        # Initialize a dictionary to store daily average results
        daily_average_results = {}
        
        for person_results in self.test_results.values():
            for day, test_result in person_results:
                if day not in daily_average_results:
                    daily_average_results[day] = []
                daily_average_results[day].append(test_result)
        
        days = []
        average_results = []
        
        for day, results_list in sorted(daily_average_results.items()):
            days.append(day)
            average_result = np.mean(results_list)
            average_results.append(average_result)
        
        plt.figure(figsize=(10, 6))
        plt.plot(days, average_results, marker='o')
        plt.title("Average Test Results Each Day")
        plt.xlabel("Day")
        plt.ylabel("Average Test Result")
        plt.grid(True)
        plt.show()

    def print_all_results(self):
        """
        Print all stored test results.
        """
        for person_index, results in self.test_results.items():
            print(f"Person {person_index} results:")
            for day, test_result in results:
                print(f"Day {day}: {test_result}")
            print()


def plot_true_IgG_levels(sim, people, people_indices):
    '''
    Plots the true IgG levels of people at specified indices over time. 

    Args: 
        sim (Sim): simulation object 
        people (People): people object (can be from sim object)
        people_indices (array): array of IDs to specify which people to plot 

    '''
    days = sim.results['t']
    for current_pathogen in range(len(sim.pathogens)): 
        IgG_levels = sim.results[current_pathogen]['IgG_level']
        for person_index in people_indices:
            person_id = people['uid'][person_index]  
            plt.plot(days, IgG_levels[:, person_index], label=f'Person {person_id}')

    plt.xlabel('Days')
    plt.ylabel('IgG Levels')
    plt.title('IgG Levels Over Time for Specific People')
    plt.legend()
    plt.grid(True)
    plt.show()

def plot_test_IgG_levels(sim, people_indices, test):
    '''
    Plots the tested IgG levels of people at specified indices. Currently for quantitative results, basically follows a 'cohort'
    (the people indices) over all the days in the simulation. 

    Args: 
        sim (Sim): simulation object 
        people_indices (array): array of IDs to specify which people to plot 
        test (Immunoassay): The test to apply on the specified people indices. 

    '''
    days = sim.results['t']
    num_days = len(days)
    num_people = len(people_indices)

    # Initialize an array to store test results for each day and each person
    test_results_over_days = [[] for _ in range(num_people)]

    for i, day in enumerate(days):
        for current_pathogen in range(len(sim.pathogens)):
            IgG_levels = sim.results[current_pathogen]['IgG_level'][i]
            test_levels = test.apply_direct(IgG_levels, people_indices)

            # Append the test results for each person to the corresponding array
            for j, person_index in enumerate(people_indices):
                test_results_over_days[j].append(test_levels[j])

    # Plot the test results over the course of all days
    plt.figure()
    for j, person_index in enumerate(people_indices):
        plt.plot(days, test_results_over_days[j], label=f'Person {person_index}')
    plt.xlabel('Day')
    plt.ylabel('Test Result')
    plt.title('Changing Test Results Over Days')
    plt.legend()
    plt.show()

def bernoulli_trials(prob_arr):
    '''
    Copy of binomial_arr function from utils

    Binomial (Bernoulli) trials each with different probabilities.

    Args:
        prob_arr (array): array of probabilities

    Returns:
         Boolean array of which trials on the input array succeeded

    **Example**::

        outcomes = bernoulli_trials([0.1, 0.1, 0.2, 0.2, 0.8, 0.8]) # Perform 6 trials with different probabilities
    '''
    return np.random.random(len(prob_arr)) < prob_arr
