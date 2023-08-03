import unittest
import behaviour as bh
import pathosim as inf 
import numpy as np

# class test_contact(unittest.TestCase):

#     def test_viralloads(self):
        
#         pars = dict(
#             pop_size=2000,
#             pop_type='behaviour_module',
#             rand_seed = 0,
#             start_day = '2020-01-03',
#             end_day = '2020-02-13',
#         )
       
#         MyNewPathogen = inf.SARS_COV_2(500)   
#         sim = inf.Sim(pars=pars, pathogens = [MyNewPathogen], verbose = 0)
#         simtest = sim.run_as_generator()
#         for t in simtest:
#             self.assertEqual(sim.people['viral_load'], sim.viral_loads)

# if __name__ == '__main__':
#     unittest.main()

pars = dict(
    pop_size=10,
    pop_type='behaviour_module',
    rand_seed = 0,
    start_day = '2020-01-03',
    end_day = '2020-02-13',
)
       
MyNewPathogen = inf.SARS_COV_2(5)   
sim = inf.Sim(pars=pars, pathogens = [MyNewPathogen], verbose = 0)
simtest = sim.run_as_generator()
for t in simtest:
    print(sim.contact_parameters)