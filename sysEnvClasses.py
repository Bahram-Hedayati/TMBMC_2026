################################################################################
# sysEnvClasses.py
################################################################################
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import heapq
from matplotlib.font_manager import FontProperties
from scipy.ndimage import gaussian_filter
from scipy.ndimage import laplace
import csv
from matplotlib.ticker import FormatStrFormatter

# The System_Class class dedicated to NMA with considerations of scalability
class System_Class:
    def __init__(self, x, y, receptor_list):
        self.x = x                          # The NMA's position on the x-axis of the grid
        self.y = y                          # The NMA's position on the y-axis of the grid
        self.bgf = -1                       # The biological gradient field (BGF)
        self.direction = -1                 # The NMA's direction
        self.last_bgf = -1                  # The last BGF value sensed by the NMA
        self.last_direction = -1            # The last direction that the NMA has traveresed
        self.find_tumor = False             # A flag to determine whether the NMA has found a cancerous voxel or not in each time step
        self.receptors = receptor_list      # List of NMA's receptors
        self.released_drug = 0              # The number of times that the NMA detects a cancerous voxel and consequently releases its therapeutic payload to cure that voxel
        self.first_cure = 0                 # The time step at which the NMA finds the first cancerous voxel
        self.last_cure = 0                  # The time step at which the NMA finds the last cancerous voxel
        self.speed = 1                      # The number of voxels that the NMA can traverse at each time step
        self.state = "Chemotaxis"           # Determines the NMA's state. Possible values are Chemotaxis and Entropy
        self.time_in_state = 0.0            # Represents the time that the NMA remains at each state over a simulation

    def get_directions_probabilities(self, bgf_list):
        """
        This function gets the oxygen (bgf) values of the neighboring voxels and returns a probability list containing
        the probability of choosing each neighboring voxel by the NMA at the next step according to the highest hypoxia:
        The maximum oxygen value assumed for a voxel is 38.0
        """
        directions_probabilities = []
        total = np.sum(bgf_list)
        for bgf in bgf_list:
            # prob = (1 - bgf) / total
            prob = (38.0 - bgf) / total
            directions_probabilities.append(prob)
        return directions_probabilities

    def get_directions_entropy(self, directions_probabilities):
        """
        This function gets a probability list containing the probability of choosing each neighboring voxel
        by the NMA at the next step according to the highest hypoxia, which was computed by
        the get_directions_probabilities function. Then, it computes the entropy for each single value.
        """
        directions_entropy = []
        for p in directions_probabilities:
            p = max(p, 0.00001)
            directions_entropy.append(-p * np.log2(p))
        return directions_entropy

    def get_directions_entropy_2(self, directions_probabilities):
        """
        This function gets a probability list containing the probability of choosing each neighboring voxel
        by the NMA at the next step according to the highest hypoxia, which was computed by
        the get_directions_probabilities function. Then, it computes the Shannon entropy from the given
        probability distribution.
        """
        entropy_value = 0
        for p in directions_probabilities:
            p = max(p, 0.00001)
            entropy_value += -p * np.log2(p)
        return entropy_value

    def choose_direction_tumble(self, bgf_list):
        """
        This function gets the substrate values inside the 8-connected voxels as the input and chooses the next direction,
        forms a probability distribution based on those values, and chooses the next direction of the NMA randomly but based on
        the computed probability distribution to simply simulate the tumble mode of e. coli bacteria.
        """
        total = np.sum(bgf_list)
        directions_probabilities = []

        for bgf in bgf_list:
            prob = (38.0 - bgf) / total
            directions_probabilities.append(prob)

        # Normalize probabilities
        probabilities_sum = np.sum(directions_probabilities)
        if probabilities_sum > 0:
            directions_probabilities = [p / probabilities_sum for p in directions_probabilities]
        else:
            # Handle the case where all probabilities are zero (e.g., uniform probability)
            directions_probabilities = [1.0 / len(bgf_list)] * len(bgf_list)

        # Select an index randomly based on the probabilities
        d = np.random.choice(len(bgf_list), p=directions_probabilities)
        if d == 0:
            direction = (-1,-1)
        elif d == 1:
            direction = (-1, 0)
        elif d == 2:
            direction = (-1, 1)
        elif d == 3:
            direction = (0,-1)
        elif d == 4:
            direction = (0, 1)
        elif d == 5:
            direction = (1,-1)
        elif d == 6:
            direction = (1, 0)
        elif d == 7:
            direction = (1, 1)

        return direction

    def choose_direction(self, bgf_list):
        """
        This function gets a list as input that contains oxygen values sensed by the NMA from all connected neighboring
        voxels. After computing the probability of NMA's movement to each voxel according to the hypoxia values,
        the entropy for each voxel is computed and the direction of a voxel with the highest entropy will return
        as the next NMA's direction.
        """
        directions_probabilities = self.get_directions_probabilities(bgf_list)           # Compute the probability of NMA's movement to the neighboring voxel at its next step, based on the BGF's gradient at each neighboring voxel
        directions_entropy       = self.get_directions_entropy(directions_probabilities) # Compute the entropy of each neighboring voxel and selects a direction with the highest value of entropy

        d = np.argmax(directions_entropy) # Extract the next direction
        if d == 0:
            direction = (-1,-1) # Left and Bottom
        elif d == 1:
            direction = (-1, 0) # Left
        elif d == 2:
            direction = (-1, 1) # Left and Above
        elif d == 3:
            direction = (0,-1) # Bottom
        elif d == 4:
            direction = (0, 1) # Above
        elif d == 5:
            direction = (1,-1) # Right and Bottom
        elif d == 6:
            direction = (1, 0) # Right
        elif d == 7:
            direction = (1, 1) # Right and Above

        return direction

    def choose_direction_entropy(self, bgf_list_chunks):
        """
        This function gets a list as input that contains oxygen values sensed by the NMA from all connected neighboring
        voxels. After computing the probability of NMA's movement to each voxel according to the hypoxia values,
        the entropy for each voxel is computed and the direction of a voxel with the highest entropy will return
        as the next NMA's direction.
        """
        entropy_values = []
        for chunk in bgf_list_chunks:
            directions_probabilities = self.get_directions_probabilities(chunk)
            entropy = self.get_directions_entropy_2(directions_probabilities)
            entropy_values.append(entropy)

        d = np.argmax(entropy_values) # Extract the next direction

        if d == 0:
            direction = (-1,-1) # Left and Bottom
        elif d == 1:
            direction = (-1, 0) # Left
        elif d == 2:
            direction = (-1, 1) # Left and Above
        elif d == 3:
            direction = (0,-1) # Bottom
        elif d == 4:
            direction = (0, 1) # Above
        elif d == 5:
            direction = (1,-1) # Right and Bottom
        elif d == 6:
            direction = (1, 0) # Right
        elif d == 7:
            direction = (1, 1) # Right and Above

        return direction

class Environment_Class:
    """
    The Environment_Class includes the following properties which are assigned by their values in the constructor of the class:
    width determines the width of the 2D grid,
    height specifies the height of the 2D grid,
    bgf_values contains the substrate diffused throughout the grid,
    tumors argument represents the position of cancerous voxels,
    cancer_voxels is a list containing the number of cancerous voxel at each iteration.
    """

    def __init__(self, width, height, bgf_values, tumors):
        self.width  = width       # Width of the lattice
        self.height = height      # Height of the lattice
        self.tumors = tumors      # Tumor mask
        self.bgf    = bgf_values  # Oxygen pressure values
        self.cancer_voxels = []

    def check_reflection(self, dx, dy, x, y, N):
        """
        This function gets dx and dy as the length of the NMA's next step in horizontal and vertical directions, respectively,
        and avoid the NMA located in x and y position to hit the edge of the grid, which has the height and width equal to N."
        """
        # if the NMA hits the edge of the grid in its next horizontal step, it will be reflected back to the same direction.
        if (x + dx) <= 0 or (x + dx) >= N-3:
            dx = -dx

        # if the NMA hits the edge of the grid in its next vertical step, it will be reflected back to the same direction.
        if (y + dy) <= 0 or (y + dy) >= N-3:
            dy = -dy

        return (dx, dy)

    def get_voxel_neighbors(self, x, y):
        """
        This function gets x and y coordinates of a voxel as the inputs to return the 8-connected neighbor voxels.
        It is called when the NMA is operating in the chemotaxis-based state of the Hybrid navigation mechanism, or while it uses the Chemotaxis naviagtion mechanism.
        """
        bgf_values = []
        bgf_values.append(self.bgf[x-1][y-1])
        bgf_values.append(self.bgf[x-1][y])
        bgf_values.append(self.bgf[x-1][y+1])
        bgf_values.append(self.bgf[x][y-1])
        bgf_values.append(self.bgf[x][y+1])
        bgf_values.append(self.bgf[x+1][y-1])
        bgf_values.append(self.bgf[x+1][y])
        bgf_values.append(self.bgf[x+1][y+1])

        return bgf_values

    def get_voxel_neighbors_entropy(self, x, y):
        """
        This function gets x and y coordinates of a reference voxel as the inputs to return the neighbor voxels.
        It is called, when the NMA is operating in the entropy-based state of the Hybrid navigation mechanism.
        For each 8 possible directions, the NMA senses a 9 * 9 patch of neighboring voxels.
        """
        bgf_list_chunks = [] # All the BGF values of the neighboring voxels, within a radius of 2 voxels from the reference voxel, are stored in this list.

        bgf_values = [] # Down Left
        bgf_values.append(self.bgf[x-2][y-2])
        bgf_values.append(self.bgf[x-2][y-1])
        bgf_values.append(self.bgf[x-2][y])
        bgf_values.append(self.bgf[x-1][y-2])
        bgf_values.append(self.bgf[x-1][y-1])
        bgf_values.append(self.bgf[x-1][y])
        bgf_values.append(self.bgf[x][y-2])
        bgf_values.append(self.bgf[x][y-1])
        bgf_values.append(self.bgf[x][y])
        bgf_list_chunks.append(bgf_values)

        bgf_values = [] # Left
        bgf_values.append(self.bgf[x-2][y-1])
        bgf_values.append(self.bgf[x-2][y])
        bgf_values.append(self.bgf[x-2][y+1])
        bgf_values.append(self.bgf[x-1][y-1])
        bgf_values.append(self.bgf[x-1][y])
        bgf_values.append(self.bgf[x-1][y+1])
        bgf_values.append(self.bgf[x][y-1])
        bgf_values.append(self.bgf[x][y])
        bgf_values.append(self.bgf[x][y+1])
        bgf_list_chunks.append(bgf_values)

        bgf_values = [] # Up Left
        bgf_values.append(self.bgf[x-2][y+2])
        bgf_values.append(self.bgf[x-2][y+1])
        bgf_values.append(self.bgf[x-2][y])
        bgf_values.append(self.bgf[x-1][y+2])
        bgf_values.append(self.bgf[x-1][y+1])
        bgf_values.append(self.bgf[x-1][y])
        bgf_values.append(self.bgf[x][y+2])
        bgf_values.append(self.bgf[x][y+1])
        bgf_values.append(self.bgf[x][y])
        bgf_list_chunks.append(bgf_values)

        bgf_values = [] # Down
        bgf_values.append(self.bgf[x-1][y-2])
        bgf_values.append(self.bgf[x-1][y-1])
        bgf_values.append(self.bgf[x-1][y])
        bgf_values.append(self.bgf[x][y-2])
        bgf_values.append(self.bgf[x][y-1])
        bgf_values.append(self.bgf[x][y])
        bgf_values.append(self.bgf[x+1][y-2])
        bgf_values.append(self.bgf[x+1][y-1])
        bgf_values.append(self.bgf[x+1][y])
        bgf_list_chunks.append(bgf_values)

        bgf_values = [] # Up
        bgf_values.append(self.bgf[x-1][y+2])
        bgf_values.append(self.bgf[x-1][y+1])
        bgf_values.append(self.bgf[x-1][y])
        bgf_values.append(self.bgf[x][y+2])
        bgf_values.append(self.bgf[x][y+1])
        bgf_values.append(self.bgf[x][y])
        bgf_values.append(self.bgf[x+1][y+2])
        bgf_values.append(self.bgf[x+1][y+1])
        bgf_values.append(self.bgf[x+1][y])
        bgf_list_chunks.append(bgf_values)

        bgf_values = [] # Down Right
        bgf_values.append(self.bgf[x+2][y-2])
        bgf_values.append(self.bgf[x+2][y-1])
        bgf_values.append(self.bgf[x+2][y])
        bgf_values.append(self.bgf[x+1][y-2])
        bgf_values.append(self.bgf[x+1][y-1])
        bgf_values.append(self.bgf[x+1][y])
        bgf_values.append(self.bgf[x][y-2])
        bgf_values.append(self.bgf[x][y-1])
        bgf_values.append(self.bgf[x][y])
        bgf_list_chunks.append(bgf_values)

        bgf_values = [] # Right
        bgf_values.append(self.bgf[x+2][y-1])
        bgf_values.append(self.bgf[x+2][y])
        bgf_values.append(self.bgf[x+2][y+1])
        bgf_values.append(self.bgf[x+1][y-1])
        bgf_values.append(self.bgf[x+1][y])
        bgf_values.append(self.bgf[x+1][y+1])
        bgf_values.append(self.bgf[x][y-1])
        bgf_values.append(self.bgf[x][y])
        bgf_values.append(self.bgf[x][y+1])
        bgf_list_chunks.append(bgf_values)

        bgf_values = [] # Up Right
        bgf_values.append(self.bgf[x+2][y+2])
        bgf_values.append(self.bgf[x+2][y+1])
        bgf_values.append(self.bgf[x+2][y])
        bgf_values.append(self.bgf[x+1][y+2])
        bgf_values.append(self.bgf[x+1][y+1])
        bgf_values.append(self.bgf[x+1][y])
        bgf_values.append(self.bgf[x][y+2])
        bgf_values.append(self.bgf[x][y+1])
        bgf_values.append(self.bgf[x][y])
        bgf_list_chunks.append(bgf_values)

        return bgf_list_chunks

