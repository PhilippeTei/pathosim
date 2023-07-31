import behaviour as bh
import pathosim as inf 
import numpy as np
  
#!!! MULTI PATHOGEN INTERACTIONS NOT COMPLETE
sum_co_inf = 0
sum_cum_deaths = 0
sum_cum_rec = 0
sims = 1

cum_inf_deadly = 0
cum_inf_covid = 0
cum_coinf_ded =0

for i in range(sims):
    n_days = 100
    pop_size = 40000
         
    #PATHOGEN 1: SARS-COV-2
    #covid = inf.SARS_COV_2(20)  
    covid = inf.SARS_COV_2(500)  
    covid3 = inf.SARS_COV_2(500)  
           

    #PATHOGEN 2: A HIGH DEATH RATE PATHOGEN
    path2 = inf.Pathogen(400, "DeadlyPathogen") #You can use Pathogen as the starting base pathogen, or any pre-programmed pathogen such as SARS-CoV-2, and then modify its parameters
    path2.beta = 0.005

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
   
    #MyNewPathogen.prognoses = MyNewPathogen.get_default_prognoses(False)
    path2.configure_prognoses(prognoses, False)  
    path2.configure_generalized_immunity(0.65, 0.95, 350, 14)
     
    #When co-infected: transmitting covid *= 0.8, transmitting deadly pathogen *= 1.5, to represent increased transmission due to the nature of the symptoms
    Mtrans = [[1,.5], #Mtrans[i,j] is row i column j
                [2,1]] #Diag should all be 1s

    Miimm = [[1,0], #innate immunity
                [1,1]] 

    Mcimm = [[1,1], #cross immunity
                [0,1]] 


    Msev = [[1,0],  
            [0,1]]

    Mdur = [[1,1],
            [1,1]]
 
    sim = inf.Sim(pop_size=pop_size, n_days=n_days, pathogens = [covid, path2,covid3], rand_seed = i)  
    sim.run()
  
    cum_inf_deadly += (sim.results[1]['cum_infections'][n_days]) 
    cum_inf_covid += (sim.results[0]['cum_infections'][n_days]) 
    sum_cum_deaths += (sim.results['cum_deaths'][n_days]) 
    sum_cum_rec += (sim.results[0]['cum_recoveries'][n_days]) 
    sum_co_inf += (sum(sim.results['co-infections'])) 
    cum_coinf_ded +=(sum(sim.results['co-infected_deaths']))  

print('cumulative co-infections', sum_co_inf/sims)   
print('cumulative deaths', sum_cum_deaths/sims) 
print('cumulative recoveries',sum_cum_rec/sims)    
print('cumulative infections deadly ',cum_inf_deadly/sims)    
print('cumulative infections covid',cum_inf_covid/sims)    
print('cumulative coinfection deaths',cum_coinf_ded/sims)    

#print(sim.results[0]['new_infections']) 
#print(sim.results[1]['new_infections']) 
#print(sim.results.keys())

#WITH RECALCULATION DISEASE TRAJ
''' 
''' 

#WITHOUT RECALCULTATIO DISEASE TRAJ
'''
MSev = 100:
cumulative co-infections 1304.0
cumulative deaths 1652.0
cumulative recoveries 43567.0
cumulative infections deadly  9414.0
cumulative infections covid 44413.0
cumulative coinfection deaths 456.0

'''