import unittest
import pathosim as inf 
import numpy as np


result_stocks = {
    'susceptible': 'Number susceptible',
    'exposed':     'Number exposed',
    'infectious':  'Number infectious',
    'symptomatic': 'Number symptomatic',
    'severe':      'Number of severe cases',
    'critical':    'Number of critical cases',
    'recovered':   'Number recovered',
    'dead':        'Number dead',
    'diagnosed':   'Number of confirmed cases',
    'known_dead':  'Number of confirmed deaths',
    'quarantined': 'Number in quarantine',
    'vaccinated':  'Number of people vaccinated',
}

result_cum_stocks = {
    'infections': 'Number infections',
    'reinfections':     'Number reinfections',
    'infectious':  'Number infectious',
    'symptomatic': 'Number symptomatic',
    'severe':      'Number of severe cases',
    'critical':    'Number of critical cases',
    'recoveries':   'Number recovered',
    'deaths':        'Number dead'
}

variant_results = ['cum_infections_by_variant',
                   'cum_symptomatic_by_variant', 'cum_severe_by_variant', 'cum_infectious_by_variant', 'cum_diagnoses_by_variant', 
                   'new_infections_by_variant', 'new_symptomatic_by_variant', 'new_severe_by_variant', 'new_infectious_by_variant', 
                   'new_diagnoses_by_variant', 'n_exposed_by_variant', 'n_infectious_by_variant']

  
class test_stratification(unittest.TestCase):
   
    def test_stratify(self):

        n_days = 150
        pop_size = 5000 
        covid = inf.SARS_COV_2(50)   
        covid.add_existing_variant("alpha", 50, 10)
        sim = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid], rand_seed = 0, verbose = 0)  
        sim.run()

         
        covid = inf.SARS_COV_2(50)   
        covid.add_existing_variant("alpha", 50, 10)
        constraint = { 
                        'age': [('>', 25), ('<=', 50)]
        }  
        sim1 = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid], enable_stratifications = True, stratification_pars = constraint, rand_seed = 0, verbose = 0)  
        sim1.run() 
        
        covid = inf.SARS_COV_2(50)   
        covid.add_existing_variant("alpha", 50, 10)
        constraint = { 
                        'age': [('<=', 25)]
        }  
        sim2 = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid], enable_stratifications = True, stratification_pars = constraint, rand_seed = 0, verbose = 0)  
        sim2.run()

        covid = inf.SARS_COV_2(50)   
        covid.add_existing_variant("alpha", 50, 10)
        constraint = { 
                        'age': [('>', 50)]
        }  
        sim3 = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid], enable_stratifications = True, stratification_pars = constraint, rand_seed = 0, verbose = 0)  
        sim3.run()
         
        assert pop_size == len(sim1.stratification_indices) + len(sim2.stratification_indices)+ len(sim3.stratification_indices)
          
        assert len(np.intersect1d(sim1.stratification_indices, sim2.stratification_indices)) == 0
        assert len(np.intersect1d(sim2.stratification_indices, sim3.stratification_indices)) == 0
        assert len(np.intersect1d(sim1.stratification_indices, sim3.stratification_indices)) == 0

        for t in range(n_days):  
            assert abs(sim.results[0]['pop_nabs'][t] - (sim1.results[0]['pop_nabs'][t]*len(sim1.stratification_indices) + sim2.results[0]['pop_nabs'][t]*len(sim2.stratification_indices)+ sim3.results[0]['pop_nabs'][t]*len(sim3.stratification_indices))/pop_size) < 0.01
            assert abs(sim.results[0]['pop_protection'][t] - (sim1.results[0]['pop_protection'][t]*len(sim1.stratification_indices) + sim2.results[0]['pop_protection'][t]*len(sim2.stratification_indices)+ sim3.results[0]['pop_protection'][t]*len(sim3.stratification_indices))/pop_size) < 0.00001
            assert abs(sim.results[0]['pop_symp_protection'][t] - (sim1.results[0]['pop_symp_protection'][t]*len(sim1.stratification_indices) + sim2.results[0]['pop_symp_protection'][t]*len(sim2.stratification_indices)+ sim3.results[0]['pop_symp_protection'][t]*len(sim3.stratification_indices))/pop_size) < 0.00001
              
            for key in result_stocks: 
                assert sim.results[0][f'n_{key}'][t] == (sim1.results[0][f'n_{key}'][t]+sim2.results[0][f'n_{key}'][t]+sim3.results[0][f'n_{key}'][t])
            for key in result_cum_stocks: 
                assert sim.results[0][f'new_{key}'][t] == (sim1.results[0][f'new_{key}'][t]+sim2.results[0][f'new_{key}'][t]+sim3.results[0][f'new_{key}'][t])
                assert sim.results[0][f'cum_{key}'][t] == (sim1.results[0][f'cum_{key}'][t]+sim2.results[0][f'cum_{key}'][t]+sim3.results[0][f'cum_{key}'][t])
            for key in variant_results: 
                for i in range(2): 
                    assert sim.results[0]['variant'][key][i][t] == (sim1.results[0]['variant'][key][i][t]+sim2.results[0]['variant'][key][i][t]+sim3.results[0]['variant'][key][i][t])

 
if __name__ == '__main__':
    unittest.main()
