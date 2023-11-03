import behaviour as bh
import pathosim as inf 
import numpy as np
  
# ---- DEMO SHOWCASING THE VACCINATION INTERVENTION IN MULTI-PATHOGEN SIMULATION ---- #

#Define general simulation parameters
n_days = 365
pop_size = 40000
         
#Define the 2 circulating pathogens
covid = inf.SARS_COV_2(500)   
            
#PATHOGEN 2: A HIGH DEATH RATE PATHOGEN
path2 = inf.Pathogen(200, "DeadlyPathogen") #You can use Pathogen as the starting base pathogen, or any pre-programmed pathogen such as SARS-CoV-2, and then modify its parameters

#Configure the pathogen's beta, prognoses (to make it deadlier) and immunity parameters
path2.beta = 0.0045
prognoses = dict(
            age_cutoffs   = np.array([0]),
            sus_ORs       = np.array([1.00]),
            trans_ORs     = np.array([1.00]),
            symp_probs    = np.array([1.00]),
            comorbidities = np.array([1.00]),
            severe_probs  = np.array([0.2]),
            crit_probs    = np.array([0.5]),
            death_probs   = np.array([0.9]),
            )
   
path2.configure_prognoses(prognoses, False)  
covid.configure_generalized_immunity(0.1, 0.95, 180, 14) #Make covid use generic system of immunity
     
#Create a vaccine for covid, includes predefined vaccines such as pfizer
#vx = inf.vaccinate_num_cov('pfizer',500, False,None, None, 1) #SPECIFY WHICH PATHOGEN THE VACCINE IS FOR (index in the list of pathogens passed to sim)

#OR Create a custom vaccine for pathogens using the generalized system for modeling immunity
#See Pathosim Documentation
vax_pars = dict(imm_peak = 1,imm_time_to_peak = 5, imm_time_to_final_value = 250, imm_final_value = 0.5) 
vx2 = inf.vaccinate_num_generic(vax_pars, 300, False,None, None, 1)  #here the vaccine is for pathogen with index 1, that uses generic immunity. 

#Run simulation
sim = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [path2, covid], interventions = [vx2],  rand_seed = 1)  
sim.run()    
sim.plot('multi-pathogen') #Plot using specific display settings for multi-pathogen. Shows co-infections.