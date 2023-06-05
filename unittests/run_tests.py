import unittest
from test_bhPop import * 
from test_dataLoaders import * 
from test_mortality import *
from test_infPop import *
from test_covSimInterventions import *
from test_covMultiSim import *
from test_covVariantsSim import *

#To run ALL unittests:
#python run_tests.py

#Run individual unittests modules:
#python unittest_file.py

disable_warnings = True

def all_tests():
     
    if disable_warnings:
        import warnings
        warnings.filterwarnings("ignore")
          
    suite = unittest.TestSuite()

    #test generation of behaviour pop
    suite.addTest(test_bhPop('test_simplePop'))

    #test loading of country households data
    suite.addTest(test_DataLoaders('test_country_households'))
    
    #test initialization of infection module population from behaviour population
    suite.addTest(test_infPop('test_different_pop_types'))

    #test simulations with very high or very low disease probability
    suite.addTest(test_diseaseMortalityTests_COVID('test_default_death_prob_one')) 
    suite.addTest(test_diseaseMortalityTests_COVID('test_default_death_prob_zero'))

    #test simulations with different variants and custom variants with differing parameters
    suite.addTest(test_covVariantsSim('test_oneVariant'))
    suite.addTest(test_covVariantsSim('test_multiVariants'))
    
    #test simulations with basic interventions (contact tracing, testing, vaccination, change_beta)
    suite.addTest(test_covSimWithInterventions('test_interventions'))

    #test multi simulations
    suite.addTest(test_covMultiSim('test_multiSim'))
     
    return suite

if __name__ == '__main__': 
    runner = unittest.TextTestRunner()
    runner.run(all_tests())
