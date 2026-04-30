import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from publicFns import *

# This file loads the single-tumor and multi-tumor tumor masks and visualize them

def plot_tumor_development(info, grid_size, tumors, title, color_bar=True):
    tick_label_size = 28
    plt.figure(figsize=(8, 8))  # Create a new figure for each plot
    plt.imshow(info, cmap='Wistia')
    if color_bar:
        cbar = plt.colorbar(shrink=0.99)  # Make the colorbar 80% of its original height
        cbar.set_label('Resistance of Voxels', labelpad=10, fontsize=tick_label_size) # Increased labelpad for more space
        # Set colorbar ticks
        cbar.set_ticks(np.arange(0, 1.1, 0.2))
        cbar.ax.tick_params(labelsize=tick_label_size)  # Increase tick font size

    y_coords, x_coords = np.where(tumors == 1)
    plt.scatter(x_coords, y_coords, color='black', label='Cancerous Voxel', marker='o', s=3)
    # plt.legend(fontsize=tick_label_size-8)
    plt.xticks(np.arange(0, grid_size+1, tick_axis), fontsize=tick_label_size)
    plt.yticks(np.arange(0, grid_size+1, tick_axis), fontsize=tick_label_size)
    plt.gca().invert_yaxis()
    plt.show()

single_tumor_mask = pd.read_csv(f"Data/tumor_mask_100.csv", header=None).to_numpy()            # Load the developed single tumor from the corresponding csv file
multiple_tumor_mask = pd.read_csv(f"Data/multiple_tumors_voxels.csv", header=None).to_numpy()  # Load the developed multiple tumors from the corresponding csv file

grid_size = 100
resistance = np.random.rand(grid_size, grid_size)  # random resistance field for visual background

# Plot the single-tumor and multi-tumor TMEs
plot_tumor_development(resistance, grid_size, single_tumor_mask, title="Resistance Values of Voxels at initial time", color_bar=True)
plot_tumor_development(resistance, grid_size, multiple_tumor_mask, title="Resistance Values of Voxels at initial time", color_bar=True)