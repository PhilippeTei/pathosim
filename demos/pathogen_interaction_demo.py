import behaviour as bh
import pathosim as inf 
import numpy as np
  
# ---- DEMO SHOWCASING MULTI-PATHOGEN SPREAD SIMULATION AND INTERACTIONS ---- #

#Define basic simulation parameters
n_days = 100
pop_size = 40000
         
#Define the first pathogen: SARS-COV-2
covid = inf.SARS_COV_2(500)   
            
#Define the second pathogen: A HIGH DEATH RATE PATHOGEN
path2 = inf.Pathogen(400, "DeadlyPathogen") #You can use Pathogen as the starting base pathogen, or any pre-programmed pathogen such as SARS-CoV-2, and then modify its parameters

#Configure the pathogen's beta, prognoses (to make it deadlier) and immunity parameters
path2.beta = 0.0065
prognoses = dict(
            age_cutoffs   = np.array([0]),
            sus_ORs       = np.array([1.00]),
            trans_ORs     = np.array([1.00]),
            symp_probs    = np.array([1.00]),
            comorbidities = np.array([1.00]),
            severe_probs  = np.array([0.2]),
            crit_probs    = np.array([0.5]),
            death_probs   = np.array([0.]),
            )
   
path2.configure_prognoses(prognoses, False)  
path2.configure_generalized_immunity(0.65, 0.95, 350, 14)
     

#Configure the matrices for pathogen-pathogen interactions
#See Pathosim Documentation: [General Parameters | Multi-pathogen simulation parameters] for details.

Mtrans = [[1,.5],
            [2,1]]

Miimm = [[1,0],
            [1,1]] 

Mcimm = [[1,1],
            [0,1]] 


Msev = [[1,2],  
        [2,1]]

Mdur = [[1,.5],
        [.5,1]]
 
#Run the simulation
sim = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid, path2], rand_seed = 1, Msev = Msev, Mdur = Mdur)  
sim.run()  
sim.plot('multi-pathogen') #Plot using specific display settings for multi-pathogen. Shows co-infections.


