import active_population_sampling as aps
import behaviour as bh
import sciris as sc
import pathosim as inf
import numpy as np

pathogen = inf.SARS_COV_2(20)
pathogen.add_existing_variant('alpha', days=5, n_imports=10)
sim = inf.Sim(pop_size=100, pop_type ='random', pathogens=[pathogen])
sim.run()


sampler_test = aps.Sampler(sim.people, {'design_type' : 'cross_sectional', 'sampling_periods' : [(0, 14), (15, 30)], 
                                        'num_people_captured': [50, 50]}, 'canadian_blood_donors')


for ID in range(sampler_test.people.pars['pop_size']):
    print(round(sampler_test.people['age'][ID]))
    people_in_interval = np.where((round(sampler_test.people['age'][ID])) > 30)


print(people_in_interval)
