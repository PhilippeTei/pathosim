import behaviour as bh
import pathosim as inf 
import numpy as np
 
# ---- DEMO SHOWCASING HOW TO CONFIGURE A GENERIC PATHOGEN'S IMMUNITY (USING THE GENERALIZED IMMUNITY SYSTEM) ---- #

pop_size = 20000 
  
#Create the circulating pathogen
MyNewPathogen = inf.SARS_COV_2(50)   

#Configure the immunity of that pathogen. 4 parameters have to be passed in this function. See Pathosim Documentation.
MyNewPathogen.configure_generalized_immunity(0.35, 0.95, 180, 14)  #Here, we pass parameters such that the spread of MyNewPathogen ressembles COVID-19's spread.

#Initialize and run the simulation
sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 365)
sim.run() 
sim.plot() 
