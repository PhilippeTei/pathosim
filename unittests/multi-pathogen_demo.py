import behaviour as bh
import pathosim as inf 
import numpy as np
 
n_days = 300
pop_size = 20000
         
#PATHOGEN 1: SARS-COV-2
covid = inf.SARS_COV_2(20)  
           

#PATHOGEN 2: A HIGH DEATH RATE PATHOGEN
MyNewPathogen = inf.Pathogen(5, "DeadlyPathogen") #You can use Pathogen as the starting base pathogen, or any pre-programmed pathogen such as SARS-CoV-2, and then modify its parameters
 
beta_dist    = dict(dist='neg_binomial', par1=1.0, par2=0.2, step=0.05) 
beta         = 0.01      
prognoses = dict(
age_cutoffs   = np.array([0]),
sus_ORs       = np.array([1.00]),
trans_ORs     = np.array([1.00]),
symp_probs    = np.array([0.75]),
comorbidities = np.array([1.00]),
severe_probs  = np.array([0.75]),
crit_probs    = np.array([0.75]),
death_probs   = np.array([0.85]))
 
MyNewPathogen.configure_transmission(beta_dist = beta_dist, beta = beta) 
MyNewPathogen.configure_prognoses(prognoses, False)
   
MyNewPathogen.add_custom_variant(label = "custom_variant1", days = 5, n_imports = 10, rel_beta = 1, rel_death_prob = 2.0)
MyNewPathogen.add_custom_variant(label = "custom_variant2", days = 15, n_imports = 10, rel_beta = 2, rel_symp_prob=1.5, rel_crit_prob=1.5, rel_death_prob=0.5) 


#Add some cross immunity between these new variants. No need to put in all values, by default, when immunity between two variants is not explicitly set, it will have by default full cross-immunity
custom_cross_immunity = dict( 

        custom_variant1 = dict(
            custom_variant2  = 0.95),

        custom_variant2 = dict(
            custom_variant1  = 0.92))

MyNewPathogen.configure_variant_cross_immunity(custom_cross_immunity)

sim = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid, MyNewPathogen])  
sim.run()
