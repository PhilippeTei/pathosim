import unittest
import behaviour as bh
import pathosim as inf 
import numpy as np

pars = dict(
    pop_size=2000,
    pop_type='behaviour_module',
    rand_seed = 0,
    start_day = '2020-01-03',
    end_day = '2020-02-13',
    enable_surveillance = True,
    enable_random_testing = True,
    random_test_start_date = '2020-01-05',
    random_test_end_date = '2020-04-05',
    random_test_frequency = 7,
)

MyNewPathogen = inf.SARS_COV_2(500)   
sim = inf.Sim(pars=pars, pathogens = [MyNewPathogen], verbose = 0)
simtest = sim.run_as_generator()


class test_dates(unittest.TestCase):

    def test_current_day(self):
        
        pars = dict(
            pop_size=2000,
            pop_type='behaviour_module',
            rand_seed = 0,
            start_day = '2020-01-03',
            end_day = '2020-02-13',
        )
        MyNewPathogen = inf.SARS_COV_2(500)   
        sim = inf.Sim(pars=pars, pathogens = [MyNewPathogen], verbose = 0)
        simtest = sim.run_as_generator()
        for t in simtest:
            self.assertEqual(sim.current_day, sim.datevec[sim.t - 1])

    def test_surveillanceEfrequency(self):
            
            syndromictestpars = dict(
                pop_size=2000,
                pop_type='behaviour_module',
                rand_seed = 0,
                start_day = '2020-01-03',
                end_day = '2020-02-13',
                enable_surveillance = True,
                enable_syndromic_testing = True,
            )

            contacttestpars = dict(
                pop_size=2000,
                pop_type='behaviour_module',
                rand_seed = 0,
                start_day = '2020-01-03',
                end_day = '2020-02-13',
                enable_surveillance = True,
                enable_contact_testing = True,
                contact_test_start_date = '2020-01-05',
                contact_test_end_date = '2020-02-05',
                contact_test_frequency = 1,
            )

            severitytestpars = dict(
                pop_size=2000,
                pop_type='behaviour_module',
                rand_seed = 0,
                start_day = '2020-01-03',
                end_day = '2020-02-13',
                enable_surveillance = True,
                enable_severity_testing = True,
                severity_test_start_date = '2020-01-05',
                severity_test_end_date = '2020-02-05',
                severity_test_frequency = 1,
            )

            agetestpars = dict(
                pop_size=2000,
                pop_type='behaviour_module',
                rand_seed = 0,
                start_day = '2020-01-03',
                end_day = '2020-02-13',
                enable_surveillance = True,
                enable_age_testing = True,
                age_test_start_date = '2020-01-05',
                age_test_end_date = '2020-02-05',
                age_test_frequency = 1,
            )

            randomtestpars = dict(
                pop_size=2000,
                pop_type='behaviour_module',
                rand_seed = 0,
                start_day = '2020-01-03',
                end_day = '2020-02-13',
                enable_surveillance = True,
                enable_random_testing = True,
                random_test_start_date = '2020-01-05',
                random_test_end_date = '2020-02-05',
                random_test_frequency = 1,
            )

            MyNewPathogen = inf.SARS_COV_2(500)
            syndromicsim = inf.Sim(pars=syndromictestpars, pathogens = [MyNewPathogen], verbose = 0)
            contactsim = inf.Sim(pars=contacttestpars, pathogens = [MyNewPathogen], verbose = 0)
            severitysim = inf.Sim(pars=severitytestpars, pathogens = [MyNewPathogen], verbose = 0)
            agesim = inf.Sim(pars=agetestpars, pathogens = [MyNewPathogen], verbose = 0)
            randomsim = inf.Sim(pars=randomtestpars, pathogens = [MyNewPathogen], verbose = 0)

            syndromicsimtest = syndromicsim.run_as_generator()
            contactsimtest = contactsim.run_as_generator()
            severitysimtest = severitysim.run_as_generator()
            agesimtest = agesim.run_as_generator()
            randomsimtest = randomsim.run_as_generator()  

            for t in contactsimtest:
                if sim.surveillance_done_today == True:
                    self.assertTrue(sim.current_day % sim.syndromic_test_frequency == 0)
                    self.assertTrue(sim.contact_test_start_date <= current_date <= sim.contact_test_end_date)

            for t in severitysimtest:
                if sim.surveillance_done_today == True:
                    self.assertTrue(sim.current_date % sim.severity_test_frequency == 0)
                    self.assertTrue(sim.severity_test_start_date <= current_date <= sim.severity_test_end_date)

            for t in agesimtest:
                if sim.surveillance_done_today == True:
                    self.assertTrue(sim.current_day % sim.age_test_frequency == 0)
                    self.assertTrue(sim.age_test_start_date <= current_date <= sim.age_test_end_date)

            for t in randomsimtest:
                if sim.surveillance_done_today == True:
                    self.assertTrue(sim.current_day % sim.random_test_frequency == 0)
                    self.assertTrue(sim.random_test_start_date <= current_date <= sim.random_test_end_date)

if __name__ == '__main__':
    unittest.main()
