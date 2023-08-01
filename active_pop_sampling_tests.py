import active_population_sampling as aps
import behaviour as bh
import sciris as sc
import matplotlib.pyplot as plt
import pathosim as inf
import numpy as np

pathogen = inf.SARS_COV_2(20)
pathogen.add_existing_variant('alpha', days=5, n_imports=10)
sim = inf.Sim(pop_size=100, pop_type ='random', pathogens=[pathogen])
sim.run()
sampler_test = aps.Sampler(sim.people, {'design_type' : 'cross_sectional', 'sampling_periods' : [(0, 14), (15, 30)], 
                                        'num_people_captured': [50, 50]}, 'canadian_blood_donors')

#different tests to run 
b_immunoassay_test_true = aps.Immunoassay(3, binary_result=True, pos_threshold=4)
b_immunoassay_test_false = aps.Immunoassay(3, binary_result=True, pos_threshold=7)
q_immunoassay_test = aps.Immunoassay(3)
LOD_under_test = aps.Immunoassay(6)
LOD_under_b_test = aps.Immunoassay(6, binary_result=True, pos_threshold=5)

#use sampler parameters to select 10 people
people_to_test = sampler_test.apply(10)
#test these 10 people 
#q_immunoassay_test.apply(sampler_test.people, people_to_test)

aps.plot_true_IgG_levels(sim, sim.people, people_to_test)
aps.plot_test_IgG_levels(sim, people_to_test, q_immunoassay_test)

#test_results = LOD_under_b_test.apply(sampler_test.people, people_to_test)

#print(test_results)


#age_array = np.round(sampler_test.people['age'])
#people_over_30 = sampler_test.people[np.where(age_array > 30)]
#print(people_over_30)


#age_array = np.round(sampler_test.people['age'])
#people_over_30 = np.where(age_array > 30)[0]  # Get the indices of individuals over 30
#filtered_people = sampler_test.people[people_over_30]


#sampler_test.assign_donation_probs()
##print(sampler_test.people['age'])
#print(sampler_test.people['uid'])
#people_to_test = sampler_test.apply(10)
#print(people_to_test)
#print(sampler_test.people['IgG_level'][people_to_test])
#print(sampler_test.people['nab'])
#print(sampler_test.people['IgG_level'])
