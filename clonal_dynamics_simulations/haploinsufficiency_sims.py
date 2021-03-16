from clone_competition_simulation.parameters import Parameters
from clone_competition_simulation.fitness_classes import *
from clone_competition_simulation.stop_conditions import EndConditionError
import matplotlib.pyplot as plt

# Sim functions
def stop_when_fully_mutant(sim):
    # Stops the simulation when the entire tissue is double Notch1 mutant.
    # Will speed up simulations
    double_mutants = sim.population_array[:, sim.plot_idx-1].toarray()[np.where(~np.isnan(sim.raw_fitness_array[:, 3]))]
    double_mutant_pop = double_mutants.sum()
    if double_mutant_pop == sim.total_pop:
        raise EndConditionError()


def sim_suff(initial_cells, mutation_rate, max_time):
    # Single allele mutations do not increase fitness
    alleles = [Gene('Notch1_allele1', FixedValue(1), synonymous_proportion=0),
               Gene('Notch1_allele2', FixedValue(1), synonymous_proportion=0)]

    # Double allele mutations have the Hom fitness
    epistatics = [('Notch1_hom', 'Notch1_allele1', 'Notch1_allele2', FixedValue(NOTCH1_HOM_FITNESS))]

    mut_gen_epi = MutationGenerator(genes=alleles, epistatics=epistatics,
                                    combine_mutations='replace', combine_array='priority',
                                    multi_gene_array=True)

    p = Parameters(algorithm='WF2D', initial_cells=initial_cells,
                   mutation_generator=mut_gen_epi,
                   mutation_rates=mutation_rate, max_time=max_time,
                   print_warnings=False, division_rate=DIVISION_RATE, samples=NUM_SAMPLES,
                   cell_in_own_neighbourhood=True, end_condition_function=stop_when_fully_mutant)
    s = p.get_simulator()
    s.run_sim()
    return s


def sim_insuff(initial_cells, mutation_rate, max_time):
    alleles = [Gene('Notch1_allele1', FixedValue(NOTCH1_HET_FITNESS), synonymous_proportion=0),
               Gene('Notch1_allele2', FixedValue(NOTCH1_HET_FITNESS), synonymous_proportion=0)]

    epistatics = [('Notch1_hom', 'Notch1_allele1', 'Notch1_allele2', FixedValue(NOTCH1_HOM_FITNESS))]

    mut_gen_epi = MutationGenerator(genes=alleles, epistatics=epistatics,
                                    combine_mutations='replace', combine_array='priority',
                                    multi_gene_array=True)

    p = Parameters(algorithm='WF2D', initial_cells=initial_cells,
                   mutation_generator=mut_gen_epi,
                   mutation_rates=mutation_rate, max_time=max_time,
                   print_warnings=False, division_rate=DIVISION_RATE, samples=NUM_SAMPLES,
                   cell_in_own_neighbourhood=True, end_condition_function=stop_when_fully_mutant)
    s = p.get_simulator()
    s.run_sim()
    return s


def get_double_mutant_pop(sim):
    double_mutants = sim.population_array.toarray()[np.where(~np.isnan(sim.raw_fitness_array[:, 3]))]
    double_mutant_pop = double_mutants.sum(axis=0)
    if double_mutant_pop.max() == sim.total_pop:
        double_mutant_pop[np.argmax(double_mutant_pop):] = sim.total_pop
    return double_mutant_pop


def get_multiple_runs(initial_cells, mutation_rate, sim_function, num_runs, max_time):
    res = np.empty((num_runs, NUM_SAMPLES + 1))
    for n in range(num_runs):
        np.random.seed(n)
        p = get_double_mutant_pop(sim_function(initial_cells, mutation_rate, max_time))
        res[n, :len(p)] = p

    return res[:, :len(p)]


def plot_results(res, colour, times, label=None, alpha=0.1):
    medians = np.median(res, axis=0)
    lower = np.quantile(res, 0.025, axis=0)
    upper = np.quantile(res, 0.975, axis=0)

    # Plot median result
    plt.plot(times, medians / GRID_SIZE, c=colour, label=label, linewidth=2, zorder=3)
    # Plot the 95% interval
    plt.fill_between(times, lower / GRID_SIZE, upper / GRID_SIZE, color=colour, alpha=alpha)


if __name__ == "__main__":
    # Parameters
    NOTCH1_HOM_FITNESS = 7.0
    NOTCH1_HET_FITNESS = 2.3
    DIVISION_RATE = 0.27
    MUTATION_RATE = 0.000005
    NUM_SIMS = 100
    NUM_SAMPLES = 200
    GRID_SIZE = 500 ** 2
    MAX_TIME = 5000

    # Run a quick simulation to get the simulation time points for plotting
    TIMES = sim_insuff(initial_cells=16, mutation_rate=0, max_time=MAX_TIME).times

    # Run haplosufficient simulations
    print('Running haplosufficient simulations')
    res_suff = get_multiple_runs(GRID_SIZE, MUTATION_RATE, sim_suff, NUM_SIMS, max_time=MAX_TIME)

    # Run haploinsufficient simulations
    print('Running haploinsufficient simulations')
    res_insuff = get_multiple_runs(GRID_SIZE, MUTATION_RATE, sim_insuff, NUM_SIMS, max_time=MAX_TIME)

    # Experiment data for plotting against the simulations
    wt_notch_neg_area = np.array([1.4, 3.1, 3.0, 6.4, 12.0]) / 100
    wt_notch_neg_area_sem = np.array([1.2, 1.8, 0.3, 3.2, 4.2]) / 100
    data_times = np.array([91, 183, 274, 365, 456])  #  3, 6, 9, 12, 15 months

    plt.figure(figsize=(6, 4))
    plt.errorbar(data_times, wt_notch_neg_area, yerr=wt_notch_neg_area_sem, zorder=5, color='k', label='Data')
    plot_results(res_insuff, colour='C1', label='Haploinsufficient', times=TIMES, alpha=0.2)
    plot_results(res_suff, colour='C2', label='Haplosufficient', times=TIMES, alpha=0.2)
    plt.ylabel('Proportion of tissue Notch1-/-')
    plt.xlabel('Time (days)')
    plt.xlim([0, 5000])
    plt.ylim([0, 1])
    plt.legend(loc='upper center', ncol=3, bbox_to_anchor=(0.5, 1.12))
    plt.tight_layout()
    plt.savefig("Haplosufficiency.pdf")
