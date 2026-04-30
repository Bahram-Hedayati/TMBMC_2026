import matplotlib.pyplot as plt
import pandas as pd
from publicFns import *

# This file provides useful functions to report or visualize the results obtained by different navigation mechanisms

def print_statistics(single_tumor, mean_effectiveness_hybrid, mean_effectiveness_chemotaxis, mean_effectiveness_random):
    print('=' * 70)
    if single_tumor:
        print("Single Tumor Scenario")
    else:
        print("Multiple Tumors Scenario")

    print("Mean Effectiveness of Hybrid: ", round(mean_effectiveness_hybrid, 2))
    print("Mean Effectiveness of Chemotaxis: ", round(mean_effectiveness_chemotaxis, 2))
    print("Mean Effectiveness of Random: ", round(mean_effectiveness_random, 2))
    print('=' * 70)

def plot_effectiveness(results, tumors_sum, title, color):
    normalized_results = [(x / tumors_sum) * 100 for x in results]

    plt.figure(figsize=(10, 5))
    plt.scatter(list(range(len(normalized_results))), normalized_results, marker='o', s=1, color=color)
    plt.grid()
    plt.title(title)
    plt.show()

def plot_compare_effectiveness(results_hybrid, results_chemotaxis, results_random, tumors_sum, title):
    normalized_hybrid = [(x / tumors_sum) * 100 for x in results_hybrid]
    normalized_chemotaxis = [(x / tumors_sum) * 100 for x in results_chemotaxis]
    normalized_random = [(x / tumors_sum) * 100 for x in results_random]

    plt.figure(figsize=(10, 5))
    plt.scatter(list(range(len(normalized_hybrid))), normalized_hybrid, label="Hybrid", marker='o', s=1, color='purple')
    plt.scatter(list(range(len(normalized_chemotaxis))), normalized_chemotaxis, label="Chemotaxis", marker='o', s=1, color='green')
    plt.scatter(list(range(len(normalized_random))), normalized_random, label="Random", marker='o', s=1, color='blue')
    plt.title(title)

    plt.grid()
    plt.legend()
    plt.show()

def plot_effectiveness_violin(results_hybrid, results_chemotaxis, results_random, title, colors):
    tick_label_size = 20

    data = [results_hybrid, results_chemotaxis, results_random]

    # Find the minimum and maximum values across all lists
    min_value = min(min(l) for l in data)
    max_value = max(max(l) for l in data)

    # Normalize each value in each list
    normalized_data = []
    for data_list in data:
        normalized_list = [((x - min_value) / (max_value - min_value)) * 100 for x in data_list]
        normalized_data.append(normalized_list)

    plt.figure(figsize=(10, 7))
    violin_parts = plt.violinplot(normalized_data, showmeans=False, showmedians=False, showextrema=True, vert=False)

    # Color the violins
    for pc, color in zip(violin_parts['bodies'], colors):
        pc.set_facecolor(color)
        pc.set_edgecolor('black')
        pc.set_alpha(0.8)  # Adjust transparency

    plt.yticks([1, 2, 3], ['Hybrid', 'Chemotaxis', 'Random'])
    plt.xticks(np.arange(0, 101, 20))
    plt.title(title)
    plt.xlabel('Mean Effectiveness Percentage (%)', fontsize=tick_label_size)
    plt.ylabel('Navigation Mechanism', fontsize=tick_label_size)
    plt.tick_params(axis='both', which='major', labelsize=tick_label_size)
    plt.grid(axis='x')
    plt.tight_layout()
    plt.show()

def plot_tumor_evolution(total_tumors, time_steps, results_tumors, color, title):
    tumor_evolution_list = []
    for i in (range(int(time_steps / 10))):
        # print(i)
        tumor_sum = 0
        for tumor_list in results_tumors:
            if np.isnan(tumor_list[i]):
                tumor_list[i] = 0
            tumor_sum += tumor_list[i]
        tumor_evolution_list.append(tumor_sum / len(results_tumors))

    cured_voxels = [(total_tumors - x) / total_tumors for x in tumor_evolution_list]

    plt.figure(figsize=(10, 5))
    plt.plot(cured_voxels, linewidth=2, color=color)
    plt.title(title)
    plt.xlabel("Time Steps")
    plt.ylabel("Number of Cured Cancerous Voxels")
    plt.grid(True)
    plt.show()

def plot_compare_tumor_evolution(results_hybrid, results_chemotaxis, results_random, time_steps, tumors_sum, title, colors):
    tick_label_size = 24

    tumor_evolution_hybrid = []
    tumor_evolution_chemotaxis = []
    tumor_evolution_random = []
    for j in (range(int(time_steps))):
        i = j * 10

        tumor_sum = 0
        for tumor_list in results_hybrid:
            if np.isnan(tumor_list[i]):
                tumor_list[i] = 0
            tumor_sum += tumor_list[i]
        tumor_evolution_hybrid.append(tumor_sum / len(results_hybrid))

        tumor_sum = 0
        for tumor_list in results_chemotaxis:
            if np.isnan(tumor_list[i]):
                tumor_list[i] = 0
            tumor_sum += tumor_list[i]
        tumor_evolution_chemotaxis.append(tumor_sum / len(results_chemotaxis))

        tumor_sum = 0
        for tumor_list in results_random:
            if np.isnan(tumor_list[i]):
                tumor_list[i] = 0
            tumor_sum += tumor_list[i]
        tumor_evolution_random.append(tumor_sum / len(results_random))

    cured_voxels_hybrid = [(tumors_sum - x) / tumors_sum for x in tumor_evolution_hybrid]
    cured_voxels_chemotaxis = [(tumors_sum - x) / tumors_sum for x in tumor_evolution_chemotaxis]
    cured_voxels_random = [(tumors_sum - x) / tumors_sum for x in tumor_evolution_random]

    plt.figure(figsize=(10, 7))
    plt.plot(cured_voxels_hybrid, linewidth=5, color=colors[0], label="Hybrid")
    plt.plot(cured_voxels_chemotaxis, linewidth=5, color=colors[1], label="Chemotaxis")
    plt.plot(cured_voxels_random, linewidth=5, color=colors[2], label="Random")
    plt.title(title)
    plt.xlabel("Time (min)", fontsize=tick_label_size)
    plt.ylabel("Fraction of Cured Voxels", fontsize=tick_label_size)
    plt.xticks(fontsize=tick_label_size)
    plt.yticks(np.arange(0, 1.1, 0.2), fontsize=tick_label_size)
    plt.legend(fontsize=tick_label_size - 2)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_compare_tumor_evolution_limited(results_single, results_multiple, time_steps, tumors_sum_single, tumors_sum_multiple, title, colors):
    tumor_evolution_single = []
    tumor_evolution_multiple = []
    for j in (range(int(time_steps))):
        i = j * 10

        tumor_sum = 0
        for tumor_list in results_single:
            if np.isnan(tumor_list[i]):
                tumor_list[i] = 0
            tumor_sum += tumor_list[i]
        tumor_evolution_single.append(tumor_sum / len(results_single))
        # tumor_evolution_single.append(tumor_sum)


        tumor_sum = 0
        for tumor_list in results_multiple:
            if np.isnan(tumor_list[i]):
                tumor_list[i] = 0
            tumor_sum += tumor_list[i]
        tumor_evolution_multiple.append(tumor_sum / len(results_multiple))
        # tumor_evolution_multiple.append(tumor_sum)


    # cured_voxels_single = [(tumors_sum_single - x) / tumors_sum_single for x in tumor_evolution_single]
    # cured_voxels_multiple = [(tumors_sum_multiple - x) / tumors_sum_multiple for x in tumor_evolution_multiple]
    cured_voxels_single = [(tumors_sum_single - x) for x in tumor_evolution_single]
    cured_voxels_multiple = [(tumors_sum_multiple - x) for x in tumor_evolution_multiple]

    tick_label_size = 24
    plt.figure(figsize=(14, 8))
    plt.plot(cured_voxels_single, linewidth=5, color=colors[0], label='Single-tumor')
    plt.plot(cured_voxels_multiple, linewidth=5, linestyle='--', color=colors[1], label='Multiple-tumors')
    # plt.title(title)
    plt.xlabel("Time (min)", fontsize=tick_label_size)
    plt.ylabel("Number of Cured Voxels", fontsize=tick_label_size)
    plt.xticks(fontsize=tick_label_size)
    plt.yticks(np.arange(0, 201, 25), fontsize=tick_label_size)
    plt.grid(True)
    plt.legend(fontsize=tick_label_size)
    plt.tight_layout()
    plt.show()

def plot_hybrid_states(results_states, time_steps, ensemble_count):
    average_states_zero = []  # Chemotaxis = 0
    average_states_minus = [] # Chemotaxis = -1
    for i in range(int(time_steps // 10)):
        sum_states_zero = 0  # Chemotaxis = 0
        sum_states_minus = 0 # Chemotaxis = -1
        for state_list in results_states:
            if state_list[i] == 0:
                sum_states_minus += -1
            sum_states_zero += state_list[i]

        average_states_zero.append(sum_states_zero / len(results_states))
        average_states_minus.append(sum_states_minus / len(results_states))

    plt.figure(figsize=(8, 4))
    plt.plot(list(np.arange(1, len(average_states_zero) + 1)), average_states_zero, linewidth=2, color='blue', label='Chemotaxis = 0')
    plt.plot(list(np.arange(1, len(average_states_minus) + 1)), average_states_minus, linewidth=2, color='red', label='Chemotaxis = -1')
    plt.title("Mean Number of Entropy and Chemotaxis States over Time Steps")
    plt.xlabel("Time Steps")
    plt.ylabel("Entropy States")
    plt.grid(True)
    plt.legend()
    plt.show()

def box_plot(results_hybrid, results_chemotaxis, results_random, title, colors):
    # Ensure data are 1-D arrays
    groups = [
        np.asarray(results_hybrid).ravel(),
        np.asarray(results_chemotaxis).ravel(),
        np.asarray(results_random).ravel()
    ]

    # Global normalization (optional)
    gmin = min(g.min() for g in groups)
    gmax = max(g.max() for g in groups)
    rng = (gmax - gmin) if gmax != gmin else 1.0
    data = [((g - gmin) / rng) * 100.0 for g in groups]

    tick_label_size = 24
    plt.figure(figsize=(9, 6))
    box = plt.boxplot(
        data,
        vert=True,                 # now vertical
        patch_artist=True,
        boxprops=dict(linewidth=1.5),
        meanprops=dict(color='black', linewidth=2),
        whiskerprops=dict(linewidth=1.2),
        capprops=dict(linewidth=1.2),
        flierprops=dict(marker='o', markersize=3, alpha=0.5)
    )

    # Color each box
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_alpha(0.8)

    plt.xticks([1, 2, 3], ['Hybrid', 'Chemotaxis', 'Random'])
    plt.yticks(np.arange(0, 101, 20), fontsize=tick_label_size)
    plt.title(title, fontsize=tick_label_size)
    plt.ylabel('Mean Effectiveness (%)', fontsize=tick_label_size)
    plt.xlabel('Navigation Mechanism', fontsize=tick_label_size)
    plt.tick_params(axis='both', which='major', labelsize=tick_label_size)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.show()

def box_plot_compare(results_hybrid_single, results_hybrid_multiple,
                     results_chemotaxis_single, results_chemotaxis_multiple,
                     results_random_single, results_random_multiple, title, colors):
    # Ensure data are 1-D arrays
    groups_hybrid = [
        np.asarray(results_hybrid_single).ravel(),
        np.asarray(results_hybrid_multiple).ravel()
    ]

    groups_chemotaxis = [
        np.asarray(results_chemotaxis_single).ravel(),
        np.asarray(results_chemotaxis_multiple).ravel()
    ]

    groups_random = [
        np.asarray(results_random_single).ravel(),
        np.asarray(results_random_multiple).ravel()
    ]

    groups_single = [
        np.asarray(results_hybrid_single).ravel(),
        np.asarray(results_chemotaxis_single).ravel(),
        np.asarray(results_random_single).ravel()
    ]

    groups_multiple = [
        np.asarray(results_hybrid_multiple).ravel(),
        np.asarray(results_chemotaxis_multiple).ravel(),
        np.asarray(results_random_multiple).ravel()
    ]

    # Global normalization (optional)
    gmin_single = min(g.min() for g in groups_single)
    gmax_single = max(g.max() for g in groups_single)

    gmin_multiple = min(g.min() for g in groups_multiple)
    gmax_multiple = max(g.max() for g in groups_multiple)

    rng_single = (gmax_single - gmin_single) if gmax_single != gmin_single else 1.0
    rng_multiple = (gmax_multiple - gmin_multiple) if gmax_multiple != gmin_multiple else 1.0

    data_single = [((g - gmin_single) / rng_single) * 100.0 for g in groups_single]
    data_multiple = [((g - gmin_multiple) / rng_multiple) * 100.0 for g in groups_multiple]
    data = [
        data_single[2], data_multiple[2],
        data_single[1], data_multiple[1],
        data_single[0], data_multiple[0]
    ]

    tick_label_size = 24
    positions = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
    plt.figure(figsize=(14, 8))
    box = plt.boxplot(
        data,
        positions=positions,
        widths=0.3,
        vert=True,
        patch_artist=True,
        boxprops=dict(linewidth=2),
        showmeans=True,
        meanline=True,
        # meanprops=dict(color='black', linewidth=4),
        medianprops=dict(linewidth=0),
        meanprops=dict(linewidth=2, color='orange', linestyle='-'),
        # meanprops=dict(marker='-', markeredgecolor='orange', markerfacecolor='orange', markersize=10),
        whiskerprops=dict(linewidth=2),
        capprops=dict(linewidth=2),
        flierprops=dict(marker='o', markersize=5, alpha=0.5)
    )

    # Color each box
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_edgecolor('black')
        patch.set_alpha(0.8)

    handles = [
        box['boxes'][0],
        box['boxes'][2],
        box['boxes'][4],
    ]

    plt.xticks(positions, ['Single', 'Multiple', 'Single', 'Multiple', 'Single', 'Multiple'])
    plt.yticks(np.arange(0, 101, 20))
    plt.title(title, fontsize=tick_label_size, fontweight='bold')
    plt.ylabel('Mean Effectiveness (%)', fontsize=tick_label_size)
    plt.xlabel('Scenario', fontsize=tick_label_size, fontweight='bold', labelpad=5)
    plt.tick_params(axis='both', which='major', labelsize=tick_label_size)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    plt.tight_layout()
    plt.legend(handles, ['Random', 'Chemotaxis', 'Hybrid'], fontsize=tick_label_size-2, loc='lower right')
    plt.show()

    # # # # # # # # # # # # Function Calls # # # # # # # # # # # #


# Since you can generate many plots using the functions defined above, it is better to comment all the function calls,
# except the ones you are interested in. Additionally, there are two flags defined to enable or disable the corresponding
# function calls for sigle-tumor and multi-tumor scenarios.

load_single_data_flag   = True # The flag that utilized to enable the function calls for the function calls, which are relevant to the single-tumor scenario.
load_multiple_data_flag = True # The flag that utilized to enable the function calls for the function calls, which are relevant to the multi-tumor scenario.

# Defining colors for the points or curves, related to each evaluated navigation mechanism, in the plots
hybrid_color     = "purple"
chemotaxis_color = "green"
random_color     = "blue"
colors = [hybrid_color, chemotaxis_color, random_color]

dt         = 0.6        # Each time step in seconds
T_single   = 500 * 60   # Total simulation time for the single-tumor scenario
T_multiple = 2000 * 60  # Total simulation time for the multi-tumor scenario

time_steps_single   = T_single / dt    # Total time steps for the single-tumor scenario
time_steps_multiple = T_multiple / dt  # Total time steps for the multi-tumor scenario

E_single   = 1000 # Total number of ensemble simulations for the single-tumor scenario
E_multiple = 1000 # Total number of ensemble simulations for the multi-tumor scenario

# Single Tumor Plots
if load_single_data_flag:

    # # # # # # # # # # # # Load the Multi-tumor Scenario's Data # # # # # # # # # # # #

    single_tumor_mask = pd.read_csv(f"Data/tumor_mask_100.csv", header=None).to_numpy() # Load the initial single tumor mask
    total_single_tumors = np.sum(single_tumor_mask) # Total number of cancerous voxels in the single-tumor scenario

    results_hybrid_single_tumor = pd.read_csv(f"Results/Single Tumor/results_hybrid_single_tumor.csv", header=None) # Load effectiveness results obtained from the hybrid navigation mechanism run in the sigle-tumor scenario
    results_tumors_hybrid_single_tumor = pd.read_csv(f"Results/Single Tumor/results_tumors_hybrid_single_tumor.csv", header=None).to_numpy() # Load efficiency results obtained from the hybrid navigation mechanism run in the sigle-tumor scenario
    results_states_hybrid_single_tumor = pd.read_csv(f"Results/Multiple_Tumors/results_states_hybrid_multiple_tumors.csv", header=None).to_numpy() # Load state dynamics data of the hybrid navigation mechanism in the single-tumor scenario

    results_chemotaxis_single_tumor = pd.read_csv(f"Results/Single Tumor/results_chemotaxis_single_tumor.csv", header=None) # Load effectiveness results obtained from the chemotaxis navigation mechanism run in the sigle-tumor scenario
    results_tumors_chemotaxis_single_tumor = pd.read_csv(f"Results/Single Tumor/results_tumors_chemotaxis_single_tumor.csv", header=None).to_numpy() # Load efficiency results obtained from the chemotaxis navigation mechanism run in the sigle-tumor scenario

    results_random_single_tumor = pd.read_csv(f"Results/Single Tumor/results_random_single_tumor.csv", header=None) # Load effectiveness results obtained from the random navigation mechanism run in the sigle-tumor scenario
    results_tumors_random_single_tumor = pd.read_csv(f"Results/Single Tumor/results_tumors_random_single_tumor.csv", header=None).to_numpy() # Load efficiency results obtained from the random navigation mechanism run in the sigle-tumor scenario

    # # # # # # # # # # # # Mean Effectiveness Results # # # # # # # # # # # #

    # Report some information about the effectiveness results for each navigation mechanism run in the single-tumor scenario
    mean_effectiveness_hybrid_single = (np.mean(results_hybrid_single_tumor.values[0]) / total_single_tumors) * 100
    mean_effectiveness_chemotaxis_single = (np.mean(results_chemotaxis_single_tumor.values[0]) / total_single_tumors) * 100
    mean_effectiveness_random_single = (np.mean(results_random_single_tumor.values[0]) / total_single_tumors) * 100
    print("Total number of cancerous voxels = ", total_single_tumors)
    print_statistics(True, mean_effectiveness_hybrid_single, mean_effectiveness_chemotaxis_single, mean_effectiveness_random_single)

    # # Plot effectiveness results over the ensemble simulations for the hybrid navigation mechanism run in the single-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Hybrid Method over {E_single} Ensemble Simulations for Single Tumor Scenario"
    plot_effectiveness(list(results_hybrid_single_tumor.values[0]), total_single_tumors, title, hybrid_color)

    # # Plot effectiveness results over the ensemble simulations for the chemotaxis navigation mechanism run in the single-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Chemotaxis Method over {E_single} Ensemble Simulations for Single Tumor Scenario"
    plot_effectiveness(list(results_chemotaxis_single_tumor.values[0]), total_single_tumors, title, chemotaxis_color)

    # # Plot effectiveness results over the ensemble simulations for the random navigation mechanism run in the single-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Random Method over {E_single} Ensemble Simulations for Single Tumor Scenario"
    plot_effectiveness(list(results_random_single_tumor.values[0]), total_single_tumors, title, random_color)

    # # Comarison of the hybrid, chemotaxis, and random navigation mechanisms, run in the single-tumor scenario, using a scatter plot.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Navigation Mechanisms over {E_single} Ensemble Simulations for Single Tumor Scenario"
    plot_compare_effectiveness(results_hybrid_single_tumor.values[0], results_chemotaxis_single_tumor.values[0], results_random_single_tumor.values[0], total_single_tumors, title)

    # # Comarison of the hybrid, chemotaxis, and random navigation mechanisms, run in the single-tumor scenario, using a violin plot.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Navigation Mechanisms over {E_single} Ensemble Simulations for Single Tumor Scenario"
    plot_effectiveness_violin(results_hybrid_single_tumor.values[0], results_chemotaxis_single_tumor.values[0], results_random_single_tumor.values[0], title, colors)

    # # Comarison of the hybrid, chemotaxis, and random navigation mechanisms, run in the single-tumor scenario, using a box plot.
    # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Single Tumor Scenario"
    box_plot(results_hybrid_single_tumor, results_chemotaxis_single_tumor, results_random_single_tumor, title, colors)

    # # # # # # # # # # # # Tumor Treatment Evolution Results # # # # # # # # # # # #

    # # Mean fraction of cured voxels by the hybrid navigation mechanism over the ensemble simulation in the single-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Ensemble Simulation Results of Hybrid over {E_single} Ensemble Simulations, \n Mean Number of Cancerous Voxels Over {int((time_steps_single * 0.6) / 60)} Minutes for Each Run"
    plot_tumor_evolution(total_single_tumors, time_steps_single, results_tumors_hybrid_single_tumor, hybrid_color, title)

    # # Mean fraction of cured voxels by the chemotaxis navigation mechanism over the ensemble simulation in the single-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Ensemble Simulation Results of Chemotaxis over {E_single} Ensemble Simulations, \n Mean Number of Cancerous Voxels Over {int((time_steps_single * 0.6) / 60)} Minutes for Each Run"
    plot_tumor_evolution(total_single_tumors, time_steps_single, results_tumors_chemotaxis_single_tumor, chemotaxis_color, title)

    # # Mean fraction of cured voxels by the random navigation mechanism over the ensemble simulation in the single-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Ensemble Simulation Results of Random over {E_single} Ensemble Simulations, \n Mean Number of Cancerous Voxels Over {int((time_steps_single * 0.6) / 60)} Minutes for Each Run"
    plot_tumor_evolution(total_single_tumors, time_steps_single, results_tumors_random_single_tumor, random_color, title)

    # # Comarison of Mean number of cured voxels using the hybrid, chemotaxis, and random navigation mechanisms, run in the single-tumor scenario.
    # To see the results, uncomment the next two lines:
    title = f"Ensemble Simulation Results of Navigation Mechanism over {E_single} Ensemble Simulations, \n Mean Number of Cured Cancerous Voxels Over {int((time_steps_single * 0.6) / 60)} Minutes for Each Run"
    plot_compare_tumor_evolution(results_tumors_hybrid_single_tumor, results_tumors_chemotaxis_single_tumor, results_tumors_random_single_tumor, 500, total_single_tumors, title, colors)

    # # Mean number of occurrence of the chemotaxis-based and entropy-based states while the NMA is operating by the hybrid navigation mechanism in the single-tumor scenario over the ensemble simulations.
    # # To see the results, uncomment the next line:
    plot_hybrid_states(results_states_hybrid_single_tumor, time_steps_single, E_single)

# Multi-tumor Plots
if load_multiple_data_flag:

    # # # # # # # # # # # # Load the Multi-tumor Scenario's Data # # # # # # # # # # # #

    multiple_tumor_mask = pd.read_csv(f"Data/multiple_tumors_voxels.csv", header=None).to_numpy() # Load the initial multiple tumors mask
    total_multiple_tumors = np.sum(multiple_tumor_mask) # Total number of cancerous voxels in the multi-tumor scenario

    results_hybrid_multiple_tumors = pd.read_csv(f"Results/Multiple_Tumors/results_hybrid_multiple_tumors.csv", header=None) # Load effectiveness results obtained from the hybrid navigation mechanism run in the multi-tumor scenario
    results_tumors_hybrid_multiple_tumors = pd.read_csv(f"Results/Multiple_Tumors/results_tumors_hybrid_multiple_tumors.csv", header=None).to_numpy() # Load efficiency results obtained from the hybrid navigation mechanism run in the multi-tumor scenario
    results_states_hybrid_multiple_tumors = pd.read_csv(f"Results/Multiple_Tumors/results_states_hybrid_multiple_tumors.csv", header=None).to_numpy() # Load state dynamics data of the hybrid navigation mechanism in the multi-tumor scenario

    results_chemotaxis_multiple_tumors = pd.read_csv(f"Results/Multiple_Tumors/results_chemotaxis_multiple_tumors.csv", header=None) # Load effectiveness results obtained from the chemotaxis navigation mechanism run in the multi-tumor scenario
    results_tumors_chemotaxis_multiple_tumors = pd.read_csv(f"Results/Multiple_Tumors/results_tumors_chemotaxis_multiple_tumors.csv", header=None).to_numpy() # Load efficiency results obtained from the chemotaxis navigation mechanism run in the multi-tumor scenario

    results_random_multiple_tumors = pd.read_csv(f"Results/Multiple_Tumors/results_random_multiple_tumors.csv", header=None) # Load effectiveness results obtained from the random navigation mechanism run in the multi-tumor scenario
    results_tumors_random_multiple_tumors = pd.read_csv(f"Results/Multiple_Tumors/results_tumors_random_multiple_tumors.csv", header=None).to_numpy() # Load efficiency results obtained from the random navigation mechanism run in the multi-tumor scenario

    # # # # # # # # # # # # Mean Effectiveness Results # # # # # # # # # # # #

    # Report some information about the effectiveness results for each navigation mechanism run in the single-tumor scenario
    mean_effectiveness_hybrid_multiple = (np.mean(results_hybrid_multiple_tumors.values[0]) / total_multiple_tumors) * 100
    mean_effectiveness_chemotaxis_multiple = (np.mean(results_chemotaxis_multiple_tumors.values[0]) / total_multiple_tumors) * 100
    mean_effectiveness_random_multiple = (np.mean(results_random_multiple_tumors.values[0]) / total_multiple_tumors) * 100
    print_statistics(False, mean_effectiveness_hybrid_multiple, mean_effectiveness_chemotaxis_multiple, mean_effectiveness_random_multiple)

    # # Plot effectiveness results over the ensemble simulations for the hybrid navigation mechanism run in the multi-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Hybrid Method over {E_multiple} Ensemble Simulations for Multiple Tumors Scenario"
    plot_effectiveness(list(results_hybrid_multiple_tumors.values[0]), total_multiple_tumors, title, hybrid_color)

    # # Plot effectiveness results over the ensemble simulations for the chemotaxis navigation mechanism run in the multi-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Chemotaxis Method over {E_multiple} Ensemble Simulations for Multiple Tumors Scenario"
    plot_effectiveness(list(results_chemotaxis_multiple_tumors.values[0]), total_multiple_tumors, title, chemotaxis_color)

    # # Plot effectiveness results over the ensemble simulations for the random navigation mechanism run in the multi-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Random Method over {E_multiple} Ensemble Simulations for Multiple Tumors Scenario"
    plot_effectiveness(list(results_random_multiple_tumors.values[0]), total_multiple_tumors, title, random_color)

    # # Comarison of the hybrid, chemotaxis, and random navigation mechanisms, run in the multi-tumor scenario, using a scatter plot.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Navigation Mechanisms over {E_multiple} Ensemble Simulations for Multiple Tumors Scenario"
    plot_compare_effectiveness(results_hybrid_multiple_tumors.values[0], results_chemotaxis_multiple_tumors.values[0], results_random_multiple_tumors.values[0], total_multiple_tumors, title)

    # # Comarison of the hybrid, chemotaxis, and random navigation mechanisms, run in the multi-tumor scenario, using a violin plot.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Navigation Mechanisms over {E_multiple} Ensemble Simulations for Multiple Tumors Scenario"
    plot_effectiveness_violin(results_hybrid_multiple_tumors.values[0], results_chemotaxis_multiple_tumors.values[0], results_random_multiple_tumors.values[0], title, colors)

    # # Comarison of the hybrid, chemotaxis, and random navigation mechanisms, run in the multi-tumor scenario, using a box plot.
    # # To see the results, uncomment the next two lines:
    title = f"Mean Effectiveness of Multiple Tumors Scenario"
    box_plot(results_hybrid_multiple_tumors, results_chemotaxis_multiple_tumors, results_random_multiple_tumors, title, colors)

    # # # # # # # # # # # # Tumor Treatment Evolution Results # # # # # # # # # # # #

    # # Mean fraction of cured voxels by the hybrid navigation mechanism over the ensemble simulation in the multi-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Ensemble Simulation Results of Hybrid over {E_multiple} Ensemble Simulations, \n Mean Number of Cancerous Voxels Over {int((time_steps_multiple * 0.6) / 60)} Minutes for Each Run"
    plot_tumor_evolution(total_multiple_tumors, time_steps_multiple, results_tumors_hybrid_multiple_tumors, hybrid_color, title)

    # # Mean fraction of cured voxels by the chemotaxis navigation mechanism over the ensemble simulation in the multi-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Ensemble Simulation Results of Chemotaxis over {E_multiple} Ensemble Simulations, \n Mean Number of Cancerous Voxels Over {int((time_steps_multiple * 0.6) / 60)} Minutes for Each Run"
    plot_tumor_evolution(total_multiple_tumors, time_steps_multiple, results_tumors_chemotaxis_multiple_tumors, chemotaxis_color, title)

    # # Mean fraction of cured voxels by the random navigation mechanism over the ensemble simulation in the multi-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Ensemble Simulation Results of Random over {E_multiple} Ensemble Simulations, \n Mean Number of Cancerous Voxels Over {int((time_steps_multiple * 0.6) / 60)} Minutes for Each Run"
    plot_tumor_evolution(total_multiple_tumors, time_steps_multiple, results_tumors_random_multiple_tumors, random_color, title)

    # # Comarison of Mean number of cured voxels using the hybrid, chemotaxis, and random navigation mechanisms, run in the multi-tumor scenario.
    # # To see the results, uncomment the next two lines:
    title = f"Ensemble Simulation Results of Navigation Mechanism over {E_multiple} Ensemble Simulations, \n Mean Number of Cured Cancerous Voxels Over {int((time_steps_multiple * 0.6) / 60)} Minutes for Each Run"
    plot_compare_tumor_evolution(results_tumors_hybrid_multiple_tumors, results_tumors_chemotaxis_multiple_tumors, results_tumors_random_multiple_tumors, 2000, total_multiple_tumors, title, colors)

    # # Mean number of occurrence of the chemotaxis-based and entropy-based states while the NMA is operating by the hybrid navigation mechanism in the multi-tumor scenario over the ensemble simulations.
    # # To see the results, uncomment the next line:
    plot_hybrid_states(results_states_hybrid_multiple_tumors, time_steps_multiple, E_multiple)

# Comparison of effectiveness and state switching dynamic results for both single-tumor and multi-tumor scenarios
if load_single_data_flag and load_multiple_data_flag:
    # # Comarison of the hybrid, chemotaxis, and random navigation mechanisms, run in the single-tumor and multi-tumor scenarios over the ensemble simulations, using a box plot.
    # To see the results, uncomment the next 5 lines:
    title = f"Mean Effectiveness of Single & Multiple Tumors Scenario"
    color_list = [random_color, random_color, chemotaxis_color, chemotaxis_color, hybrid_color, hybrid_color]
    box_plot_compare(results_hybrid_single_tumor, results_hybrid_multiple_tumors,
                     results_chemotaxis_single_tumor, results_chemotaxis_multiple_tumors,
                     results_random_single_tumor, results_random_multiple_tumors, title, color_list)

    # # Comarison of Mean number of cured voxels using the hybrid, chemotaxis, and random navigation mechanisms over the ensemble simulations, run in the single-tumor and multi-tumor scenario.
    # To see the results, uncomment the next 4 lines:
    equalized_hybrid_multiple = results_tumors_hybrid_multiple_tumors[:, :int(time_steps_single/10)] # Selects the same number of the single-tumor results
    title = "Comparison of Hybrid over the same time duration"
    colors = [hybrid_color, hybrid_color]
    plot_compare_tumor_evolution_limited(results_tumors_hybrid_single_tumor, equalized_hybrid_multiple,
                                         int(T_single/60), total_single_tumors, total_multiple_tumors, title, colors)

print("All the requested plots generated successfully.")