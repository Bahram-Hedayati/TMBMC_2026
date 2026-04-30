################################################################################
# publicFns.py
# This file includes functions required by other files.
################################################################################
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.ticker import FormatStrFormatter

# This file includes

tick_axis = 20       # a global variable that specifies the tick step on each axis of the resulting plot.
tick_label_size = 38 # a global variable that specifies the font size of each tick step on axis of the resulting plot.

def plot_tumor_growth(info, grid_size, tumors, title, color_bar=True):
    """
    This function gets the following arguments as input to show the simulated tumor on the tissue:
    info is the resistance values,
    grid_size specifies the size of the 2D grid of voxels,
    tumors contains the mask specifying the position of each cancerous voxel on the grid,
    title is the title of figure,
    color_bar is a flag to decide whether the color bar should be shown next to the figure or not.
    """

    plt.figure(figsize=(8, 8))
    plt.imshow(info, cmap='Wistia')
    plt.title(title)

    if color_bar:
        cbar = plt.colorbar(shrink=0.8)                     # Make the colorbar 80% of its original height
        cbar.set_label('Resistance of Voxels', labelpad=10) # Increased labelpad for more space
        cbar.set_ticks(np.arange(0, 1.1, 0.1))              # Set colorbar ticks

    y_coords, x_coords = np.where(tumors == 1)                                               # Extract the positions of cancerous voxels from the tumor mask named as tumors.
    plt.scatter(x_coords, y_coords, color='black', label='Cancerous Voxel', marker='o', s=3) # Cancerous voxels will be illustrated on the plot as black filled circles with size s.
    plt.legend()
    plt.xticks(np.arange(0, grid_size+1, tick_axis)) # tick_axis is the global variable showing the steps of ticks
    plt.yticks(np.arange(0, grid_size+1, tick_axis)) # tick_axis is the global variable showing the steps of ticks
    plt.gca().invert_yaxis()                         # Enforce both axis to start from 0 at the bottom left hand side of the figure
    plt.show()

def plot_agent_trajectory(bgf, title, trajectory, grid_size):
    """
    This function gets the following arguments as input to show the voxels traversed by the NMA to detect the tumor:
    bgf is the biological gradient field values, like oxygen pressure levels, on the 2D grid,
    title is the title of figure,
    trajectory contains the position of voxels traveresed by the NMA on the grid,
    grid_size is the size of the squared grid.
    """
    plt.figure(figsize=(7, 7))
    plt.imshow(bgf, cmap="Blues", vmin=np.min(bgf), vmax=np.max(bgf))             # Shows the bgf values on the grid
    cbar = plt.colorbar(shrink=0.8)                                               # Make the colorbar 80% of its original height
    cbar.set_label('Oxygen Levels', fontsize=13, fontweight='bold', labelpad=10)  # Adjusted font size and weight
    plt.title(title)

    plt.scatter([trajectory[0][1]], [trajectory[0][0]], color='pink', label='Start') # Specify the voxel as the NMA's start posititon by pink filled circle
    plt.scatter([trajectory[-1][1]], [trajectory[-1][0]], color='red', label='End')  # Specify the voxel as the NMA's end posititon by red filled circle

    trajectory_x, trajectory_y = zip(*trajectory)                                               # Extract x and y values of the trajectory traversed by the NMA
    plt.plot(trajectory_y, trajectory_x, marker='.', color='purple', linewidth=1, markersize=1) # Shows each traversed voxel by a purple point on the 2D grid

    plt.xticks(np.arange(0, grid_size+1, tick_axis), fontsize=13) # tick_axis is the global variable showing the steps of ticks
    plt.yticks(np.arange(0, grid_size+1, tick_axis), fontsize=13) # tick_axis is the global variable showing the steps of ticks
    plt.gca().invert_yaxis()

def plot_multiple_trajectories(bgf, title, chemotaxis_path, entropy_path, tumor_sites, grid_size):
    """
    This function gets the following arguments as input to show the voxels traversed by the NMA at chemotaxis and entropy operating modes, separately:
    bgf is the biological gradient field values, like oxygen pressure levels, on the 2D grid,
    title is the title of figure,
    chemotaxis_path shows the trajectory traversed by the NMA while it had been in the chemotaxis mode,
    entropy_path shows the trajectory traversed by the NMA while it had been in the entropy mode,
    tumor_sites contains the position of cancerous voxels,
    grid_size is the size of the squared grid.
    """
    plt.figure(figsize=(8, 8))
    plt.imshow(bgf, cmap="Blues", vmin=np.min(bgf), vmax=np.max(bgf), zorder=0)               # Shows the bgf values on the grid
    cbar = plt.colorbar(shrink=0.99)                                                          # Make the colorbar 80% of its original height
    cbar.set_label('Oxygen Level', fontsize=tick_label_size, fontweight='bold', labelpad=10)  # Adjusted font size and weight
    cbar.ax.tick_params(labelsize=tick_label_size)                                            # tick_label_size is the global variable showing font size of the ticks on the color bar
    plt.title(title)

    cbar.formatter.set_useOffset(False)                           # Turn off offset display
    cbar.formatter.set_scientific(False)                          # Turn off scientific notation
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')) # Enforce the color bar ticks to be shown with only one floating point

    # Compute an appropriate tick_step value
    tick_min = np.min(bgf)
    tick_max = np.max(bgf)
    tick_step = (tick_max - tick_min) / 6.0
    cbar.set_ticks(np.arange(tick_min, tick_max, tick_step))

    # Extract x and y values of the trajectory traversed by the NMA while it had operated in chemotaxis mode
    if len(chemotaxis_path) > 0:
        path_x_chemotaxis, path_y_chemotaxis = zip(*chemotaxis_path)
    else:
        path_x_chemotaxis, path_y_chemotaxis = [], []

    # Extract x and y values of the trajectory traversed by the NMA while it had operated in entropy mode
    if len(entropy_path) > 0:
        path_x_entropy, path_y_entropy = zip(*entropy_path)
    else:
        path_x_entropy, path_y_entropy = [], []

    plt.plot(path_y_chemotaxis, path_x_chemotaxis, marker='.', color='yellow', linewidth=1, markersize=1, label='Chemotaxis', zorder=1) # Shows each traversed voxel by a yellow point on the 2D grid for the chemotaxis mode
    plt.plot(path_y_entropy, path_x_entropy, marker='.', color='orange', linewidth=1, markersize=1, label='Entropy', zorder=1)          # Shows each traversed voxel by an orange point on the 2D grid for the entropy mode

    plt.scatter([chemotaxis_path[0][1]], [chemotaxis_path[0][0]], color='pink', label='Start', zorder=2) # Specify the voxel as the NMA's start posititon by pink filled circle
    plt.scatter([chemotaxis_path[-1][1]], [chemotaxis_path[-1][0]], color='red', label='End', zorder=2)  # Specify the voxel as the NMA's end posititon by red filled circle

    y_coords, x_coords = np.where(tumor_sites == 1)                           # Extract the positions of cancerous voxels from the tumor mask named as tumors.
    plt.scatter(x_coords, y_coords, color='black', marker='o', s=1, zorder=3) # Cancerous voxels will be illustrated on the plot as black filled circles with size s.

    plt.xticks(np.arange(0, grid_size+1, tick_axis), fontsize=tick_label_size) # tick_label_size is the global variable showing font size of the ticks on the grid
    plt.yticks(np.arange(0, grid_size+1, tick_axis), fontsize=tick_label_size) # tick_label_size is the global variable showing font size of the ticks on the grid
    plt.gca().invert_yaxis()

def plot_oxygen_diffusion(info, grid_size, tumors, title, tumor_flag=True):
    """
    This function gets the following arguments as input to show how the biological gradient field (BGF) diffuses throughout the grid:
    info is the BGF values,
    grid_size specifies the size of the 2D grid of voxels,
    tumors contains the mask specifying the position of each cancerous voxel on the grid,
    title is the title of figure,
    tumor_flag is used to decide whether the cancerous voxels should be shown or not.
    All font sizes are set to the value of global variable tick_label_size.
    """
    plt.figure(figsize=(8, 8))                                # Create a new figure for each plot
    plt.imshow(info, cmap='Blues')                            # Shows the bgf values on the grid
    plt.title(title, fontsize=tick_label_size, fontweight='bold', pad=15)  # Shows the figure title

    cbar = plt.colorbar(shrink=0.99)                                                         # Make the colorbar 80% of its original height
    cbar.set_label('Oxygen Level', fontsize=tick_label_size, fontweight='bold', labelpad=10) # Increased labelpad for more space
    cbar.ax.tick_params(labelsize=tick_label_size)                                           # tick_label_size is the global variable showing font size of the ticks on the color bar

    # Turn off scientific notation
    cbar.formatter.set_useOffset(False)                           # Turn off offset display
    cbar.formatter.set_scientific(False)                          # Turn off scientific notation
    cbar.ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f')) # Enforce the color bar ticks to be shown with only one floating point

    # Compute an appropriate tick_step value
    tick_min = np.min(info)
    tick_max = np.max(info)
    tick_step = (tick_max - tick_min) / 5.0
    cbar.set_ticks(np.arange(tick_min, tick_max))
    cbar.set_ticks(np.arange(tick_min, tick_max, tick_step))

    # Font Adjustments
    # font = FontProperties()
    # font.set_size(tick_label_size = 14)
    # # font.set_weight('bold')
    # for label in cbar.ax.get_xticklabels():
    #     label.set_fontproperties(font)
    #     # label.set_fontsize(14) # Adjusted font size
    #
    # for label in cbar.ax.get_yticklabels():
    #     label.set_fontproperties(font)
    #     # label.set_fontsize(14) # Adjusted font size
    #
    #     # Set bold for tick values on the y-axis
    # for label in plt.gca().get_yticklabels():
    #     label.set_fontproperties(font)

    # Decide whether the cancerous voxels should be shown or not
    if tumor_flag:
        y_coords, x_coords = np.where(tumors == 1)
        plt.scatter(x_coords, y_coords, color='black', label="Cancerous Voxel", marker='o', s=3)
        # plt.legend()

        # Font Adjustments
        # font = FontProperties()
        # # font.set_weight('bold')
        # font.set_size(13)
        # # font.set_weight('bold')
        # for label in plt.gca().get_xticklabels():
        #     label.set_fontproperties(font)
        #
        # # Set bold for tick values on the y-axis
        # for label in plt.gca().get_yticklabels():
        #     label.set_fontproperties(font)

    plt.gca().invert_yaxis()
    plt.xticks(np.arange(0, grid_size+1, tick_axis), fontsize=tick_label_size)
    plt.yticks(np.arange(0, grid_size+1, tick_axis), fontsize=tick_label_size)
    plt.tick_params(axis='both', which='major', labelsize=tick_label_size) # You can adjust the labelsize value
    # plt.show(block=False) # If you want to see multiple figures consequently, you should comment this command

def show_result(navigation, results, tumor_distribution, ensemble_count):
    """
    This function gets the following arguments as input to show how the results numerically and visually:
    navigation is a text representing the NMA's navigation strategy,
    results contains the number of cured voxels at each round of ensemble simulations,
    tumor_distribution shows the tumor mask,
    ensemble_count determines the number of ensemble simulations.
    """
    print(f"Cured Voxels Results of {navigation}:")
    print(f"Mean: {np.mean(results)} out of {np.sum(tumor_distribution)}")                           # Shows the mean number of cured voxels out of the total number of cancerous voxels over the entire ensemble simulations
    print(f"Effectiveness: {round(float(np.mean(results) / np.sum(tumor_distribution)) * 100, 2)}%") # Reports the percentage of cancerous voxels detected by the NMA
    print(f"St. Deviation: {np.std(results)}")                                                       # Reports the standard deviation of cured voxels

    # Shows number of cured voxels at each run of ensemble simulations
    plt.figure(figsize=(12, 5))
    plt.scatter(list(range(len(results))), results, label=navigation, marker='o', s=3)
    plt.title(f"Ensemble Simulations, E = {ensemble_count} of TDD using {navigation} Navigation")
    plt.grid(axis='y')
    plt.show()