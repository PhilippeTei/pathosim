import behaviour as bh
import pathosim as inf 
import numpy as np

# ---- DEMO SHOWCASING HOW TO CREATE AND CONFIGURE AN ARBITRARY PATHOGEN ---- #

#simulation with random population 
#set up basic parameters for the simulation  
pop_size = 10000
initial_infections = 100


#Create a Pathogen object
#You can use Pathogen as the starting base pathogen, or any pre-programmed pathogen such as SARS-CoV-2, and then modify its parameters
MyNewPathogen = inf.Pathogen(initial_infections, "MyNewPathogen")


#Configure the pathogen, change parameters you want to change 
#Lets configure some transmission parameters and some prognoses. Other parameters can be configured, see Pathosim Documentation: [Pathogen-Specific Parameters section]
beta_dist    = dict(dist='neg_binomial', par1=1.0, par2=0.2, step=0.05) 
beta         = 0.015     

prognoses = dict(
age_cutoffs   = np.array([0]),
sus_ORs       = np.array([1.00]),
trans_ORs     = np.array([1.00]),
symp_probs    = np.array([0.75]),
comorbidities = np.array([1.00]),
severe_probs  = np.array([0.250]),
crit_probs    = np.array([0.5]),
death_probs   = np.array([1.0]))


MyNewPathogen.configure_transmission(beta_dist = beta_dist, beta = beta) 
MyNewPathogen.configure_prognoses(prognoses, False)

 

#Creating some variants with different relative transmission and probabilities of disease outcomes.
MyNewPathogen.add_custom_variant(label = "custom_alpha", days = 5, n_imports = 10, rel_beta = 1, rel_death_prob = 2.0)
MyNewPathogen.add_custom_variant(label = "custom_beta", days = 15, n_imports = 10, rel_beta = 2, rel_symp_prob=1.5, rel_crit_prob=1.5, rel_death_prob=0.5)
MyNewPathogen.add_custom_variant(label = "custom_gamma", days = 25, n_imports = 10, rel_beta = 1.5, rel_crit_prob=2.0)


#Configuring some cross-variant immunity between these new variants created above. 
#No need to explicetly write in all values, by default, when immunity between two variants is not defined in the dict, 
#a default value of 1 (full cross-variant immunity) is assumed
custom_cross_immunity = dict( 

        custom_alpha = dict(
            custom_beta  = 0.95,
            custom_gamma = 0.92),

        custom_gamma = dict(
            custom_beta  = 0.96,
            custom_alpha = 0.98), 

        custom_beta = dict(custom_gamma  = 0.95))

MyNewPathogen.configure_variant_cross_immunity(custom_cross_immunity)

#Configure the immunity model of the arbirary pathogen
MyNewPathogen.configure_generalized_immunity(0.15, 0.9, 180, 10)

#Initialize and run the simulation, along with a plot
sim = inf.Sim(pop_size = pop_size, pathogens = [MyNewPathogen], n_days = 150)
sim.run()
sim.plot('variants')