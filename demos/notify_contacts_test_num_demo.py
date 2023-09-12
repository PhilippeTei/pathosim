import behaviour as bh
import pathosim as inf 
import numpy as np
import sciris as sc


n_days = 365
pop_size = 40000
            
bh_pars = sc.objdict(n = pop_size, 
                        country_location = "canada")
        
pop = bh.BehaviourModel(bh_pars)


covid = inf.SARS_COV_2(500)   

tp = inf.test_prob(symp_prob=0.2, asymp_prob=0.001, symp_quar_prob=1.0, asymp_quar_prob=1.0, do_plot=False)
ct = inf.contact_tracing(trace_probs=dict(h=1.0, s=0.5, w=0.5, c=0.3), do_plot=False)  
 

sim = inf.Sim(pop_size=pop_size, people = pop, n_days=n_days,interventions = [tp, ct],  pathogens = [covid])  
sim.run()   
     
 