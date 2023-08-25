import covasim as cv

pars = dict(
pop_size=20000,
pop_type='random',
n_days = 50,
pop_infected = 3,
rand_seed = 0,
start_day = '2020-01-30',
)

sim = cv.Sim(pars=pars, verbose = 1)
sim.run()