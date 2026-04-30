try:
    # Local imports for shared modules
    from sysEnvClasses import *
    from publicFns import *
    from TME_init import *

except ImportError:
    # Fallback to relative imports if used as a package
    from .publicFns import *
    from .sysEnvClasses import *
    from TME_init import *

# Read the size of lattice (grid_size), as well as tumor mask values and substrate values, which are distributed throughout the lattice.
grid_size, tumor_mask, substrate = tme_single_tumor_single_agent_init()

# Show the substrate distribution throughout the lattice at the steady-state.
plot_oxygen_diffusion(substrate, grid_size, tumor_mask, f"Tumor Microenvironment at Steady State", tumor_flag=True)

# Delta pairs, including delta_x and delta_y, defined for the NMA to move based on them.
directions = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

# List of the NMA's receptors. The NMA is able to sense the substrate values of its 8-connected neighboring voxels,
# as well as the host voxel. To define some levels of noise, we can distort this list by reducing the number of receptors.
receptors = [1, 2, 3, 4, 5, 6, 7, 8, 9]

# Simulation Parameters
dx = 100.0                 # Length of each voxel in μm
dt = 0.6                   # Duration of each time step in second
T = 500 * 60               # Total simulation time per run in seconds
D = 3333.3333              # Oxygen diffusion coefficient in μm²/s
lambda_decay = 0.0008333   # Oxygen uptake rate in sec⁻¹
time_steps = int(T / dt)   # Total number of time steps
E = 1000                   # The number of ensemble simulations

results_chemotaxis = []        # The list containing the effectiveness of chemotaxis navigation mechanism measured at the end of each run
results_tumors_chemotaxis = [] # The list containing the number of cancerous voxels remained at the end of each run

for e in np.arange(1, E + 1): # This loop handles ensemble simulations
    # For each ensemble simulation, the initial location of the NMA is choosen randomly w.r.t a margin of 3 voxels
    initial_x = random.randint(3, 96)
    initial_y = random.randint(3, 96)

    # Instantiate an NMA object and a TME object from the System and Environment classes
    sc = System_Class(initial_x, initial_y, receptors)                                          # Define the system (NMA) object
    env_chemotaxis = Environment_Class(grid_size, grid_size, substrate.copy(), tumor_mask.copy()) # Define the environment (TME) object

    # Report the mean effectiveness every 10 ensmble simutations
    if e % 10 == 0:
        mean_effectiveness = (np.mean(results_chemotaxis) / np.sum(tumor_mask)) * 100
        print(f"e = {e}, Mean Effectiveness = {mean_effectiveness}")

    for t in np.arange(1, time_steps + 1): # This loop handles oxygen diffusion and the NMA's movement
        # Environmental Dynamics, including oxygen diffusion according to the presence of cancerous voxels
        laplacian = laplace(env_chemotaxis.bgf, mode='constant', cval=38.0)      # Laplacian operator represents the diffusion process
        diffusion_term = D * laplacian / dx ** 2                                 # Diffusion term w.r.t each voxel
        decay_term = -lambda_decay * env_chemotaxis.bgf * env_chemotaxis.tumors  # Decay oxygen pressure only in tumor voxels
        env_chemotaxis.bgf += dt * (diffusion_term + decay_term)                 # Explicit update
        env_chemotaxis.bgf = np.clip(env_chemotaxis.bgf, 0, 38.0)                # Clamp to physical range of oxygen pressure

        # Applying the substrate diffusion time scale. The NMA's movement time scale is 10 times slower than that of the substrate diffusion
        if t % 10 > 0:
            continue

        bgf_list = env_chemotaxis.get_voxel_neighbors(sc.x, sc.y) # Gets the substrate values inside the 8-connected neighboring voxels
        sc.direction = sc.choose_direction_tumble(bgf_list)       # Chooses the next direction of the NMA based on the tumble mode of the chemotaxis-based state

        # Movement
        # Apply boundary conditions for the NMA: it reflects if it reaches the boundary of the lattice
        sc.direction = env_chemotaxis.check_reflection(sc.direction[0], sc.direction[1], sc.x, sc.y, grid_size)

        sc.x = sc.x + sc.direction[0] # Update the position of the NMA along the x-axis
        sc.y = sc.y + sc.direction[1] # Update the position of the NMA along the y-axis

        # Check the voxel
        x, y = sc.x, sc.y                                # Uses temporary variables x and y to avoid repetition of "sc."
        sc.bgf = env_chemotaxis.bgf[x][y]                # Senses the substrate value inside the current voxel
        sc.find_tumor = env_chemotaxis.tumors[x][y] == 1 # Checks whether the current voxel is cancerous or not

        if sc.find_tumor: # If the NMA recognizes a cancerous voxel
            sc.released_drug += 1           # Releases the drug (counts the number of treated cancerous voxel)
            env_chemotaxis.tumors[x][y] = 0 # Removes the diagnosed cancerous voxel from the tumor mask
            env_chemotaxis.bgf[x][y] = 38.0 # Rising the oxygen value inside the voxel to the maximum, due to releasing drug

        # Stores the number of remained cancerous voxels at the end of the current iteration
        env_chemotaxis.cancer_voxels.append(np.sum(env_chemotaxis.tumors))

    results_chemotaxis.append(sc.released_drug)                     # Stores the number of treated voxels over each ensemble simulation
    results_tumors_chemotaxis.append(env_chemotaxis.cancer_voxels)  # Stores the number of remained cancerous voxels over each ensemble simulation

# Save results into csv files
import csv

# Save results_chemotaxis to a CSV file
with open('results_chemotaxis_single_tumor.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(results_chemotaxis)

# Save results_tumors_chemotaxis to a CSV file
with open('results_tumors_chemotaxis_single_tumor.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(results_tumors_chemotaxis)

print("Results saved to the corresponding csv files.")