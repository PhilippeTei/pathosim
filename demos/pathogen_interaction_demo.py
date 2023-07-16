import behaviour as bh
import pathosim as inf 
import numpy as np
 
n_days = 300
pop_size = 20000
         
#PATHOGEN 1: SARS-COV-2
covid = inf.SARS_COV_2(20)  
           

#PATHOGEN 2: A HIGH DEATH RATE PATHOGEN
MyNewPathogen = inf.Pathogen(10, "DeadlyPathogen") #You can use Pathogen as the starting base pathogen, or any pre-programmed pathogen such as SARS-CoV-2, and then modify its parameters
MyNewPathogen.beta = 0.0075

prognoses = dict(
            age_cutoffs   = np.array([0]),
            sus_ORs       = np.array([1.00]),
            trans_ORs     = np.array([1.00]),
            symp_probs    = np.array([1.00]),
            comorbidities = np.array([1.00]),
            severe_probs  = np.array([0.2]),
            crit_probs    = np.array([0.28]),
            death_probs   = np.array([0.5]))
   

MyNewPathogen.configure_prognoses(prognoses, False) 
MyNewPathogen.configure_generalized_immunity(0.65, 0.95, 350, 14)

#When co-infected: transmitting covid *= 0.8, transmitting deadly pathogen *= 1.5, to represent increased transmission due to the nature of the symptoms
Mtrans = [[1,.5], #Mtrans[i,j] is row i column j
          [2,1]] #Diag should all be 1s

Miimm = [[1,0], #innate immunity
          [1,1]] 

Mcimm = [[1,1], #cross immunity
          [0,1]] 

sim = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid, MyNewPathogen], Mtrans = Mtrans, Miimm = Miimm, Mcimm = Mcimm)  
sim.run()
