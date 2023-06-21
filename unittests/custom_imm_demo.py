import behaviour as bh
import infection as inf 
import numpy as np
 
pop_size = 1000 
 
#300 days immunity
MyNewPathogen = inf.SARS_COV_2(50)   
MyNewPathogen.configure_generalized_immunity(0, 1, 20, 7) 

sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 100)
sim.run() 
sim.plot()


 
#Short immunity (20 days)

MyNewPathogen = inf.SARS_COV_2(1000)   
MyNewPathogen.configure_generalized_immunity(0, 1, 20, 7) 
MyNewPathogen.configure_imports([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50])#at day 50, infect 50 people
cb = inf.change_beta([0,100], [0,1])

prognoses = dict(
    age_cutoffs   = np.array([0]),
    sus_ORs       = np.array([1.00]),
    trans_ORs     = np.array([1.00]),
    symp_probs    = np.array([0.75]),
    comorbidities = np.array([1.00]),
    severe_probs  = np.array([0.250]),
    crit_probs    = np.array([0.1]),
    death_probs   = np.array([0.001]))

MyNewPathogen.configure_prognoses(prognoses)

simb = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 200, interventions = [cb])
simb.run()   
simb.plot() 


