import behaviour as bh
import pathosim as inf 
import numpy as np

# Sim Parameters
Sim1Pars = dict(
    pop_size=20,
    pop_type='random',
    n_days = 10,
    rand_seed = 0,
    start_day = '2020-01-01'
)

Sim2Pars = dict(
    pop_size=20,
    pop_type='random',
    n_days = 10,
    rand_seed = 0,
    start_day = '2020-01-01'
)

# Initialize Sims
Sim1 = inf.Sim(pars=Sim1Pars, pathogens = [inf.SARS_COV_2(4)], verbose = 0, label = "Simulation without high-contact layer")
Sim1.initialize()
Sim2 = inf.Sim(pars=Sim2Pars, pathogens = [inf.SARS_COV_2(4)], verbose = 0, label = "Simulation with high-contact layer")
Sim2.initialize()

# new high-contact layer
n_people = len(Sim2.people)
top_pct = 0.5  # Represents x%
num_adjusted_people = int(n_people * top_pct)
additional_contacts = 5  # Represents y

# Manually specify p1 and p2 for the new layer for the top 25% 
p1 = []
p2 = []
for i in range(num_adjusted_people):
    for _ in range(additional_contacts):
        p1.append(i)
        # You can choose contacts randomly from the rest of the population or use other criteria.
        # Here, I am selecting them sequentially for simplicity.
        contact = (i + 1 + _) % n_people
        p2.append(contact)

beta = np.ones(len(p1))
layer = inf.Layer(p1=p1, p2=p2, beta=beta)

# Add  high-contact layer to Sim2
Sim2.people.contacts.add_layer(highcontact=layer)
Sim2.reset_layer_pars()
Sim2.label = f'High-contact layer with manually specified contacts for top 25%'

# #Double the contacts for Sim1
# p1_double = np.concatenate([p1, p1])
# p2_double = np.concatenate([p2, p2])
# beta_double = np.ones(len(p1_double))
# layer1 = inf.Layer(p1=p1_double, p2=p2_double, beta=beta_double)
# Sim1.people.contacts.add_layer(layer_name=layer1)
# Sim1.reset_layer_pars()

# Run the simulations
Sim1.run()
Sim2.run()


#TESTING

#Print the results

# print(
#     (((Sim2.contact_parameters) / (Sim2.pars['n_days']+1))/2) -
#     (Sim1.contact_parameters) / (Sim1.pars['n_days']+1)
# )

#print cumulative infections
print(Sim1.results[0]['cum_infections'].values)
print(Sim2.results[0]['cum_infections'].values)
