import behaviour as bh
import pathosim as inf 
import numpy as np
 
n_days = 300
pop_size = 20000
          
covid = inf.SARS_COV_2(100)  
        
#Define conditions for stratification
#List of conditions for that attribute (age), first element of a tuple is the operator, second is the value to compare the attribute against
constraints = {
                'age': [('>', 50), ('<=', 80)], 
                'sex': [('=', 0)]
}
 
sim = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid], enable_stratifications = True, stratification_pars = constraints)  
sim.run()
 