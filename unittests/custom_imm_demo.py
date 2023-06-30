import behaviour as bh
import infection as inf 
import numpy as np
 
pop_size = 20000 
  
MyNewPathogen = inf.SARS_COV_2(50)   
MyNewPathogen.configure_generalized_immunity(0.35, 0.95, 180, 14)  
sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 500)
sim.run() 
print(sim.people.imm_level[0])
sim.plot()

   
MyNewPathogen = inf.SARS_COV_2(50)     
sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 500)
sim.run() 
sim.plot()

