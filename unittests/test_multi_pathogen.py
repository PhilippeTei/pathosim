
import unittest
import pathosim as inf 
import numpy as np

class test_multi_pathogen_sim(unittest.TestCase):
   
    def test_multi_cov(self):
        """
        Check simulations with multiple COVID-19, different immunity systems, different combinations. No interaction or cross immunity enabled!!
        """
      
       

        n_days = 100 
        pop_size = 10000
         
        #Test with 3 COVID-19 using the nab system
        pathogen1 = inf.SARS_COV_2(50, "Covid 1") 
        pathogen2 = inf.SARS_COV_2(50, "Covid 2") 
        pathogen3 = inf.SARS_COV_2(50, "Covid 3") 
           
        sim1 = inf.Sim(pop_size=pop_size, n_days=n_days, verbose = 0, pathogens = [pathogen1, pathogen2,pathogen3])  
        sim1.run()
          
         
        #Test 3 COVID-19 with a mix of nab and generalized immunity systems
        pathogen1 = inf.SARS_COV_2(50, "Covid 1") 
        pathogen1.configure_generalized_immunity(0.1, 0.95, 180, 14)
        pathogen2 = inf.SARS_COV_2(50, "Covid 2") #uses nabs
        pathogen3 = inf.SARS_COV_2(50, "Covid 3") 
        pathogen3.configure_generalized_immunity(0.1, 0.95, 180, 14)
        

        sim2 = inf.Sim(pop_size=pop_size, n_days=n_days, verbose = 0, pathogens = [pathogen1, pathogen2, pathogen3])  
        sim2.run()

        
        #Test 3 COVID-19 with a only the generalized immunity system
        pathogen1 = inf.SARS_COV_2(50, "Covid 1") 
        pathogen1.configure_generalized_immunity(0.1, 0.95, 180, 14)
        pathogen2 = inf.SARS_COV_2(50, "Covid 2")  
        pathogen2.configure_generalized_immunity(0.1, 0.95, 180, 14)
        pathogen3 = inf.SARS_COV_2(50, "Covid 3") 
        pathogen3.configure_generalized_immunity(0.1, 0.95, 180, 14)
         
        sim3 = inf.Sim(pop_size=pop_size, n_days=n_days, verbose = 0, pathogens = [pathogen1, pathogen2, pathogen3])  
        sim3.run()

      
        assert sim1.results[0].keys() == sim1.results[1].keys()
        assert sim1.summary[0].keys() == sim1.summary[1].keys()
        assert sim1.results[1].keys() == sim1.results[2].keys()
        assert sim1.summary[1].keys() == sim1.summary[2].keys() 
           
        #If one of these fail retry the simulation, if it fails again then theres likely to be a bug

        for i in range(3):
            assert approx(sim1.summary[i]['cum_infections'], 11150, 2000) == True
            assert approx(sim1.summary[i]['cum_reinfections'], 2250, 500) == True
            assert approx(sim1.summary[i]['cum_infectious'], 11200, 2000) == True
            assert approx(sim1.summary[i]['cum_symptomatic'], 7150, 2000) == True
            assert approx(sim1.summary[i]['cum_recoveries'], 11000, 2000) == True
            assert approx(sim1.summary[i]['cum_severe'], 500, 250) == True
            assert approx(sim1.summary[i]['cum_critical'], 70, 50) == True
            
            assert approx(sim2.summary[i]['cum_infections'], 11150, 2000) == True
            assert approx(sim2.summary[i]['cum_reinfections'], 2250, 500) == True
            assert approx(sim2.summary[i]['cum_infectious'], 11200, 2000) == True
            assert approx(sim2.summary[i]['cum_symptomatic'], 7150, 2000) == True
            assert approx(sim2.summary[i]['cum_recoveries'], 11000, 2000) == True
            assert approx(sim2.summary[i]['cum_severe'], 500, 250) == True
            assert approx(sim2.summary[i]['cum_critical'], 70, 50) == True

            assert approx(sim3.summary[i]['cum_infections'], 11150, 2000) == True
            assert approx(sim3.summary[i]['cum_reinfections'], 2250, 500) == True
            assert approx(sim3.summary[i]['cum_infectious'], 11200, 2000) == True
            assert approx(sim3.summary[i]['cum_symptomatic'], 7150, 2000) == True
            assert approx(sim3.summary[i]['cum_recoveries'], 11000, 2000) == True
            assert approx(sim3.summary[i]['cum_severe'], 500, 250) == True
            assert approx(sim3.summary[i]['cum_critical'], 70, 50) == True
       

          
    
    def test_multi_pathos(self):
        #not asserting anything, used for printing for now
        pop_size = 20000
        initial_infections = 100
         
        MyNewPathogen = inf.Pathogen(initial_infections, "DeadlyNewPathogen") #You can use Pathogen as the starting base pathogen, or any pre-programmed pathogen such as SARS-CoV-2, and then modify its parameters 
        beta_dist    = dict(dist='neg_binomial', par1=1.0, par2=0.2, step=0.05) 
        beta         = 0.0320      
        prognoses = dict(
        age_cutoffs   = np.array([0]),
        sus_ORs       = np.array([1.00]),
        trans_ORs     = np.array([1.00]),
        symp_probs    = np.array([0.75]),
        comorbidities = np.array([1.00]),
        severe_probs  = np.array([0.50]),
        crit_probs    = np.array([0.50]),
        death_probs   = np.array([1.0]))
         
        MyNewPathogen.configure_transmission(beta_dist = beta_dist, beta = beta) 
        MyNewPathogen.configure_prognoses(prognoses, False) 
         
        MyNewPathogen.add_custom_variant(label = "custom_alpha", days = 5, n_imports = 25, rel_beta = 1, rel_death_prob = 2.0)
        MyNewPathogen.add_custom_variant(label = "custom_beta", days = 15, n_imports = 40, rel_beta = 2, rel_symp_prob=1.5, rel_crit_prob=1.5, rel_death_prob=0.5)
        MyNewPathogen.add_custom_variant(label = "custom_gamma", days = 25, n_imports = 80, rel_beta = 1.5, rel_crit_prob=2.0)
          
        custom_cross_immunity = dict( 

                custom_alpha = dict(
                    custom_beta  = 0.9,
                    custom_gamma = 0.1),

                custom_gamma = dict(
                    custom_beta  = 0.2,
                    custom_alpha = 0.1), 

                custom_beta = dict(custom_gamma  = 0.1))

        MyNewPathogen.configure_variant_cross_immunity(custom_cross_immunity)
 
        MyNewPathogen.configure_generalized_immunity(0.5, 1, 300, 14)


        COVID = inf.SARS_COV_2(100)
        sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen, COVID], n_days = 150)
        sim.run()


  
def approx(value, desired, delta):
    return True if (abs(value-desired)<delta) else False

if __name__ == '__main__':
    unittest.main()
