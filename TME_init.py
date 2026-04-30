################################################################################
# TME_init.py (restored from notebook)
################################################################################
import numpy as np
import pandas as pd

# This file includes two functions that initialize the lattice for each single-tumor and multi-tumor environments.

def tme_multi_tumor_single_agent_init():
    # This function loads and returns grid size, cancerous voxel (tumor_mask), and the biological gradient field, like Oxygen, for multiple tumors scenario.
    gridSize = 100
    tumor_mask = pd.read_csv(f"Data/multiple_tumors_voxels.csv", header=None)
    substrate = pd.read_csv(f"Data/oxygen_multi_tumor.csv", header=None)
    tumor_mask = tumor_mask.to_numpy()
    substrate = substrate.to_numpy()
    return gridSize, tumor_mask, substrate

def tme_single_tumor_single_agent_init():
    # This function loads and returns grid size, cancerous voxel (tumor_mask), and the biological gradient field, like Oxygen, for the single tumor scenario.
    gridSize = 100
    tumor_mask = pd.read_csv(f"Data/tumor_mask_100.csv", header=None)
    substrate = pd.read_csv(f"Data/oxygen_single_tumor.csv", header=None)
    tumor_mask = tumor_mask.to_numpy()
    substrate = substrate.to_numpy()
    return gridSize, tumor_mask, substrate