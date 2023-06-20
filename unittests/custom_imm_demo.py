import behaviour as bh
import infection as inf 
import numpy as np
 
pop_size = 1000 
 
#300 days immunity
MyNewPathogen = inf.SARS_COV_2(1000)   
MyNewPathogen.configure_generalized_immunity(0.2, 1, 300, 7) 
MyNewPathogen.configure_imports([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50])#at day 50, infect 50 people
cb = inf.change_beta([0,50], [0,1])

sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 150, interventions = [cb])
sim.run() 
sim.plot()


 
#Short immunity (20 days)

MyNewPathogen = inf.SARS_COV_2(1000)   
MyNewPathogen.configure_generalized_immunity(0.2, 1, 20, 7) 
MyNewPathogen.configure_imports([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,50])#at day 50, infect 50 people
cb = inf.change_beta([0,50], [0,1])

simb = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 150, interventions = [cb])
simb.run()   
simb.plot() 