#A FILE WHICH CONTAINS PARS FOR EACH OF THE PREDEFINED SAMPLING FRAMES
'''
Set the parameters for the population sampling by mapping the name of the sample 
frame to the demographic breakdown. 
'''

sample_frame_mapping = {
    'canadian_blood_donors': {
        'age_intervals' : [(17, 19), (20, 29), (30, 39), (40, 49), (50, 59), (60, 69), (80, 100)], #check max age of agents
        'donor_breakdown_per_interval' : [0.02, 0.17, 0.18, 0.16, 0.21, 0.19, 0.06, 0.01], #based on 2020 data from Statistics Canada
        'sex_breakdown' : {'male' : 0.53, 'female' : 0.47}, 
        'donors_in_population' : 0.04, 
        'num_donors_per_year' : 765809, #based on 2020 data from Statistics Canada
        'redonation_prob' : 0.5
    }, 

    'antenatal': {},
    
    'households' : {}, 

    'preexisting_cohort' : {}


}
