from pyabc import (ABCSMC, RV, Distribution, PNormDistance)
from pyabc.sampler import MulticoreEvalParallelSampler
from functools import partial, update_wrapper
import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
from collections import OrderedDict
from clone_competition_simulation.parameters import Parameters
import pyabc.visualization.credible as credible
import sys

DATA_FILE = "09-03-21 Final clonal counting dataset.xlsx"

# Functions and fixed parameters for simulations
GRID_SHAPE = (500, 500)
DIVISION_RATE = 0.27
CELLS = GRID_SHAPE[0]*GRID_SHAPE[1]

ERROR_OBJECT = {'distance': 99999}
NUM_CLONES = 100

# Each 500x500 grid is 0.25 million cells.
# If 1million cells in the oesophagus, then 4 grids is a full oesophagus.
# Asking for 100 clones in total for each time point.
# Run simulations repeatedly until enough clones exist for each time point.
# Don't want to get stuck in a loop if the induction/fitness is too low and no clones survive
# Limit the number of loops. If 4 grids per oeosophagus, then 50 grids is 12.5 mice worth, more than sampled number.
LOOP_LIMITS = 50


def get_grid(fitness, induction, grid_shape, cells):
    # Make the initial grid with randomly placed induced cells
    initial_grid = np.zeros(grid_shape, dtype=int)
    total_mutants = int(induction*cells)
    mutant_locs = np.random.choice(grid_shape[0]*grid_shape[1], total_mutants, replace=False)
    mutant_locs = [(m // grid_shape[1], m % grid_shape[1]) for m in mutant_locs]  # Convert the random draws into array indices

    count = 0
    for i in range(total_mutants):
        initial_grid[mutant_locs[count]] = i + 1
        count += 1

    fitness_array = [1] + [fitness]*total_mutants
    label_array = [0] + [1]*total_mutants
    return initial_grid, fitness_array, label_array


def distance_ks(target, sim_results):
    return ks_2samp(target, sim_results).statistic


def run_sim(parameters, target_data, return_clone_sizes=False):
    fitness, induction = parameters['fitness'], parameters['induction']
    times = [t for t in target_data]

    try:
        initial_grid, fitness_array, label_array = get_grid(fitness, induction, GRID_SHAPE, CELLS)
        if len(fitness_array) == 1:  # Induction rate too low. No mutants on grid.
            return ERROR_OBJECT
        p = Parameters(algorithm='WF2D', initial_grid=initial_grid, times=times, fitness_array=fitness_array,
                       label_array=label_array,
                       print_warnings=False, division_rate=DIVISION_RATE,
                       cell_in_own_neighbourhood=True)
        clone_counts = np.zeros(len(times), dtype=int)
        sim_distributions = [[] for i in range(len(times))]
        if p:
            for loop in range(LOOP_LIMITS):
                s = p.get_simulator()
                s.run_sim()
                for i, t in enumerate(times):
                    clone_sizes = s.get_clone_sizes_array_for_non_mutation(t=t, index_given=False, exclude_zeros=True,
                                                                       label=1)
                    clone_counts[i] += len(clone_sizes)
                    sim_distributions[i].extend(clone_sizes)

                if np.all(clone_counts >= NUM_CLONES):
                    # Have got enough clones
                    break

            if not np.all(clone_counts >= NUM_CLONES):
                # Reached loop limit before getting required number of clones.
                return ERROR_OBJECT

            res = OrderedDict([
                (t, np.random.choice(sim_distributions[i], size=NUM_CLONES, replace=False)) for i, t in enumerate(times)
            ])
            if return_clone_sizes:
                return res

            total_distance = sum([distance_ks(target_data[t], res[t]) for t in target_data])
            return {'distance': total_distance}
        else:
            return ERROR_OBJECT
    except (Exception, SystemExit) as e:
        print('Error')
        print(e)
        return ERROR_OBJECT


def _read_basal_counts(df):
    return OrderedDict([(t, _read_basal_counts_column(df[t])) for t in df.columns])


def _read_basal_counts_column(col):
    return np.array([v for v in col.values if (not np.isnan(v) and v > 0)], dtype=int)


def load_data(target_data, data_file):
    target_data = target_data.lower()
    if target_data == 'wt':
        cols = "C,E,G,I,K"
        times = [10, 14, 28, 63, 91]
    elif target_data == 'het':
        cols = "M,O,Q,S"
        times = [10, 28, 63, 91]
    elif target_data == 'het_ctl':
        cols = "U,W,Y,AA"
        times = [10, 28, 63, 91]
    elif target_data == 'hom':
        cols = "AC,AE,AG"
        times = [10, 14, 28]
    elif target_data == 'hom_ctl':
        cols = "AI,AK,AM"
        times = [10, 14, 28]
        
    df = pd.read_excel(data_file, sheet_name="Splitted data table", skiprows=8, 
                       usecols=cols, header=None, names=times, engine='openpyxl')
    return _read_basal_counts(df)


if __name__ == "__main__":
    # Set up ABC run
    priors = Distribution(
            fitness=RV("uniform", 0, 50),
            induction=RV("uniform", 0, 0.1)
        )

    distance = PNormDistance()
    sampler = MulticoreEvalParallelSampler(n_procs=10)  # Can change the number of cores used here

    # Pick the data to fit to
    target_data = sys.argv[1]   # Use 'wt', 'het_ctl', 'hom_ctl', 'het' or 'hom'
    target_data_options = ['wt', 'het_ctl', 'hom_ctl', 'het', 'hom']
    if target_data.lower() not in target_data_options:
        print('Run script with one of the following arguments: {}'.format(target_data_options))
        sys.exit(1)
    f = partial(run_sim, target_data=load_data(target_data, DATA_FILE))
    update_wrapper(f, run_sim)  # Adds the __name__ attribute that pyabc uses

    # Run the fitting using PyABC
    abc = ABCSMC(f, priors, distance, population_size=1000, sampler=sampler)
    db_path = ("sqlite:///" + target_data+'_pyabc.db')
    abc.new(db_path, {'distance': 0})
    history = abc.run(minimum_epsilon=0.1, max_nr_populations=15)


    # The database is stored and can be reloaded if needed - see pyabc documentation
    # Here just print the median and 95% credible intervals.
    def get_estimate_and_ci_for_param(param, df, w, confidence=0.95):
        vals = np.array(df[param])
        lb, ub = credible.compute_credible_interval(vals, w, confidence)
        median = credible.compute_quantile(vals, w, 0.5)
        return {'median': median, 'CI_lower_bound': lb, 'CI_upper_bound':ub}


    df, w = history.get_distribution()  # Get the accepted parameters from the last generation

    for p in ['fitness', 'induction']:
        print(p, get_estimate_and_ci_for_param(p, df, w))
