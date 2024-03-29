import unittest
import test_bhPop as t1
import test_dataLoaders as t2 
import test_mortality as t4
import test_covSimInterventions as t5 
import test_covVariantsSim as t6
import test_mergeStates as t7
import test_multi_pathogen as t8
import test_multi_region as t9
import test_stratification as t10

#To run ALL unittests:
#python run_tests.py

#Run individual unittests modules:
#python unittest_file.py

disable_warnings = True
generate_baselines = False

def all_tests():
     
    if disable_warnings:
        import warnings
        warnings.filterwarnings("ignore")
          
    suite = unittest.TestSuite()

    #test generation of behaviour pop
    suite.addTest(t1.test_bhPop('test_simplePop'))

    #test loading of country households data
    suite.addTest(t2.test_DataLoaders('test_country_households'))
     
    #test simulations with very high or very low disease probability
    suite.addTest(t4.test_diseaseMortalityTests_COVID('test_default_death_prob_one')) 
    suite.addTest(t4.test_diseaseMortalityTests_COVID('test_default_death_prob_zero'))

    #test simulations with different variants and custom variants with differing parameters
    suite.addTest(t6.test_covVariantsSim('test_oneVariant'))
    suite.addTest(t6.test_covVariantsSim('test_multiVariants'))
    
    #test simulations with basic interventions (contact tracing, testing, vaccination, change_beta)
    suite.addTest(t5.test_covSimWithInterventions('test_changeBeta'))
    suite.addTest(t5.test_covSimWithInterventions('test_vaccinateNumAndBoosters'))

    #test merge states functions 
    suite.addTest(t7.test_mergeStates("test_mergeStates"))

    #test multi-pathogen simulations 
    suite.addTest(t8.test_multi_pathogen_sim("test_multi_cov"))
     
    #test multi-region module 
    suite.addTest(t9.test_multi_region("test_multi_reg"))
     
    #test stratification module 
    suite.addTest(t10.test_stratification("test_stratify"))
     
     
    return suite

if __name__ == '__main__': 

    if(not generate_baselines):
        runner = unittest.TextTestRunner()
        runner.run(all_tests())
    else:
        print("Generating baseline files for all unittests which require one")
         
        t1.generate_baseline()  
        t5.generate_baseline() 
        t6.generate_baseline()
        t9.generate_baseline()
