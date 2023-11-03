import behaviour as bh
import pathosim as inf 
import numpy as np
 
# ---- DEMO SHOWCASING THE STRATIFICATION OF RESULTS SYSTEM ---- #

#Define general simulation parameters
n_days = 300
pop_size = 20000
          
#Create a pathogen, here, SARS-CoV-2
covid = inf.SARS_COV_2(100)  
        
#Define conditions for stratification
#List of conditions for that attribute (age), first element of a tuple is the operator, second is the value to compare the attribute against
#See Pathosim Documentation.
constraints = {
                'age': [('>', 50), ('<=', 80)], 
                'sex': [('=', 0)]
}
 
#Run the simulation and see stratified results (plotting as well)
sim = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid], enable_stratifications = True, stratification_pars = constraints)  
sim.run()
 