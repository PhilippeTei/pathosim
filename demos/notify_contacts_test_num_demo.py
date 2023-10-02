import behaviour as bh
import pathosim as inf 
import numpy as np
import sciris as sc


n_days = 365
pop_size = 20000
            
bh_pars = sc.objdict(n = pop_size, 
                        country_location = "canada")
        
pop = bh.BehaviourModel(bh_pars)


covid = inf.SARS_COV_2(100)  
second_path = inf.SARS_COV_2(100)
second_path.label = 'Covid 2' 

#create testing interventions and specify as an argument the pathogen it is for
tn = inf.test_num(daily_tests=500, pathogen = 0)
tn1 = inf.test_num(daily_tests=500, pathogen = 1) 
#tp = inf.test_prob(symp_prob=0.2, asymp_prob=0.001, symp_quar_prob=1.0, asymp_quar_prob=1.0, do_plot=False, pathogen = 0)
#tp2 = inf.test_prob(symp_prob=0.2, asymp_prob=0.001, symp_quar_prob=1.0, asymp_quar_prob=1.0, do_plot=False, pathogen = 1) 

#create the contact tracing intervention and specify which pathogens are affected as an array of pathogens

#ct = inf.contact_tracing(trace_probs=dict(h=1.0, s=0.5, w=0.5, c=0.3), do_plot=False, pathogens = [0,1])   #Quarantining after notified
nc = inf.notify_contacts(trace_probs={'h': 1, 's': 0.2, 'w': 0.2, 'c': 0.03},                               #Just notification
                                  trace_time={'h': 0, 's': 1, 'w': 1, 'c': 2}, pathogens = [0])

sim = inf.Sim(pop_size=pop_size, people = pop, n_days=n_days, interventions = [tn1, tn], pathogens = [covid, second_path])  
sim.run()     
 