#A FILE WHICH CONTAINS PARS FOR EACH OF THE PREDEFINED SAMPLING FRAMES
'''
Set the parameters for the population sampling by mapping the name of the sample 
frame to the demographic breakdown. 
'''

sample_frame_mapping = {
    'canadian_blood_donors': {
        'age_intervals' : [(17, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 69), (70, 79), (80, 100)], #check max age of agents
        'donor_breakdown_per_interval' : [0.02, 0.17, 0.18, 0.16, 0.21, 0.19, 0.06, 0.01], #based on 2020 data from Statistics Canada
        'sex_breakdown' : {'male' : 0.53, 'female' : 0.47}, 
        'donors_in_population' : 0.04, #dont know if this is even necessary 
        'num_donors_per_day' : 765809//365, #based on 2020 data from Statistics Canada
        'redonation_prob' : 0.9, #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9017805/
    }, 
#  Statistics Canada age breakdown of pregnancies (assuming all pregnancies will provide an antenatal sample which could be included in the sampling frame).
#  Taken from most recent available data (2021) (ref). 

    'antenatal': {
        'age_intervals' : [(12, 15), (15, 19), (20, 24), (25, 29), (30, 34), (35, 39), (40, 44), (45, 49)], 
        'donor_breakdown_per_interval' : [0.000160, 0.0135, 0.0863, 0.265, 0.383, 0.208, 0.0419, 0.00264], 
        'sex_breakdown' : {'male' : 0, 'female' : 1},
        'cohort_size' : 370155,
        'actual_pop_size' : 8566164

    },
    
    'households' : {}, 

    'CanPath_preexisting_cohort' : {
        #BC generations project 
        'age_intervals' : [(35, 44), (45, 54), (55, 64), (64, 100)],
        'donor_breakdown_per_interval' : [0.133, 0.271, 0.405, 0.191],
        'sex_breakdown' : {'male' : 0.31, 'female' : 0.69}, 
        'cohort_size' : 30000,
        'actual_pop_size' : 5000000 #BC population size 
        #insert drop out rate 

    }


}