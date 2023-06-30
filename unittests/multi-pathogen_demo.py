import behaviour as bh
import infection as inf 
import numpy as np
 
pop_size = 20000 
 
#300 days immunity
COVID = inf.SARS_COV_2(50)    



sim = inf.Sim(pop_size = pop_size, pathogens = [COVID], n_days = 200)
sim.run() 
sim.plot()


'''
pop_size = 20000  
#300 days immunity
MyNewPathogen = inf.SARS_COV_2(50)    
sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 200)
sim.run() 
sim.plot()
'''
