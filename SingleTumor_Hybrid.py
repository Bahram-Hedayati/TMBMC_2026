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

# The number of NMA's receptors. The NMA is able to sense the substrate values of its 8-connected neighboring voxels,
# as well as the host voxel. To define some levels of noise, we can distort this list by reducing the number of receptors.
receptors = [1, 2, 3, 4, 5, 6, 7, 8, 9]

# Simulation Parameters
dx = 100.0                 # Length of each voxel in µm
dt = 0.6                   # Duration of each time step in second
T = 500 * 60               # Total simulation time per run (in seconds)
D = 3333.3333              # Oxygen diffusion coefficient in μm²/s
lambda_decay = 0.0008333   # Oxygen uptake rate in sec⁻¹
time_steps = int(T / dt)   # Total number of time steps
E = 10                   # The number of ensemble simulations

results_hybrid = []        # Stores the effectiveness value of the NMA for each ensemble simulation
results_tumors_hybrid = [] # Stores the number of cancerous voxels remained at the end of each ensemble simulation
results_hybrid_states = [] # For each ensemble simulation, it stores the list of the NMA's state over each ensemble simulation

M = 10 # Controls the switching dynamics between the chemotaxis-based and entropy-based states

for e in np.arange(1, E + 1): # This loop handles ensemble simulations
    # For each ensemble simulation, the initial location of the NMA is choosen randomly w.r.t a margin of 3 voxels
    initial_x = random.randint(3, 96)
    initial_y = random.randint(3, 96)

    # Instantiate an NMA object and a TME object from the System and Environment classes
    sc = System_Class(initial_x, initial_y, receptors)                                        # Define the system (NMA) object
    env_hybrid = Environment_Class(grid_size, grid_size, substrate.copy(), tumor_mask.copy()) # Define the environment (TME) object

    results_states = [] # clear the list of states at the beginning of each ensemble simulation

    # Report the mean effectiveness every 10 ensmble simutations
    if e % 10 == 0:
        mean_effectiveness = (np.mean(results_hybrid) / np.sum(tumor_mask)) * 100
        print(f"e = {e}, Mean Effectiveness = {mean_effectiveness}")

    for t in np.arange(1, time_steps + 1): # This loop handles oxygen diffusion and the NMA's movement
        # Environmental Dynamics, including oxygen diffusion according to the presence of cancerous voxels
        laplacian = laplace(env_hybrid.bgf, mode='constant', cval=38.0)  # Laplacian operator represents the diffusion process
        diffusion_term = D * laplacian / dx ** 2                         # Diffusion term for further usage
        decay_term = -lambda_decay * env_hybrid.bgf * env_hybrid.tumors  # Decay oxygen pressure only in tumor voxels
        env_hybrid.bgf += dt * (diffusion_term + decay_term)             # Explicit update
        env_hybrid.bgf = np.clip(env_hybrid.bgf, 0, 38.0)                # Clamp to physical range of oxygen pressure

        # Applying the substrate diffusion time scale. The NMA's movement time scale is 10 times slower than that of the substrate diffusion
        if t % 10 > 0:
            continue

        # Depending on the current state of the NMA (either chemotaxis-based or entropy-based),
        # the NMA senses the related neighboring voxels and chooses its direction accordingly.
        if sc.state == "Chemotaxis": # The chemotaxis-based state
            bgf_list = env_hybrid.get_voxel_neighbors(sc.x, sc.y)  # Gets the substrate values inside the 8-connected neighboring voxels
            sc.direction = sc.choose_direction_tumble(bgf_list)    # Chooses the next direction of the NMA based on the tumble mode of the chemotaxis-based state

        else: # The entropy-based state
            # A sensing set of 5 * 5 neighboring voxel w.r.t the reference voxel, including the current voxel is selected.
            # These values are stored according to each 8 possible direction in eight patches of voxels separately.
            bgf_list_chunks = env_hybrid.get_voxel_neighbors_entropy(sc.x, sc.y)
            sc.direction = sc.choose_direction_entropy(bgf_list_chunks)

        # Movement
        # Applying boundary conditions for the NMA: it reflects if it reaches the boundary of the lattice
        sc.direction = env_hybrid.check_reflection(sc.direction[0], sc.direction[1], sc.x, sc.y, grid_size)

        sc.x = sc.x + sc.direction[0] # Update the position of the NMA along the x-axis
        sc.y = sc.y + sc.direction[1] # Update the position of the NMA along the y-axis

        # Check the voxel
        x, y = sc.x, sc.y                             # Uses temporary variables x and y to avoid repetition of "sc."
        sc.bgf = env_hybrid.bgf[x][y]                 # Senses the substrate value inside the current voxel
        sc.find_tumor = env_hybrid.tumors[x][y] == 1  # Checks whether the current voxel is cancerous or not

        if sc.find_tumor: # If the NMA recognizes a cancerous voxel
            sc.state = "Entropy"         # The NMA switches to the entropy-state
            m = 0                        # Resets the counter that controls the switching dynamics between chemotaxis-based and entropy-based states
            sc.released_drug += 1        # Releases the drug (counts the number of treated cancerous voxel)
            env_hybrid.tumors[x][y] = 0  # Removes the diagnosed cancerous voxel from the tumor mask
            env_hybrid.bgf[x][y] = 38.0  # Rising the oxygen value inside the voxel to the maximum, due to releasing drug

        else: # If the NMA recognizes the current voxel as normal
            if sc.state == "Entropy":
                m += 1                      # counts the number of times that the NMA in the entropy-state and has not find a cancerous voxel (stored in m)
                if m >= M:                  # if the m value reaches the maximum number of time defined as M for the NMA, it switches to the chemotaxis-based state
                    sc.state = "Chemotaxis"

        # Stores the state of the NMA (1 and 0 denotes the entropy-based and chemotaxis-based states, respectively) in the current iteration
        if sc.state == "Entropy":
            results_states.append(1)
        else:
            results_states.append(0)

        # Stores the number of remained cancerous voxels at the end of the current iteration
        env_hybrid.cancer_voxels.append(np.sum(env_hybrid.tumors))

    results_hybrid.append(sc.released_drug)                 # Stores the number of treated voxels over each ensemble simulation
    results_tumors_hybrid.append(env_hybrid.cancer_voxels)  # Stores the number of remained cancerous voxels over each ensemble simulation
    results_hybrid_states.append(results_states)            # Stores the switching dynamics between the chemotaxis-based and the entropy-based states

# # Save results into csv files
# import csv
#
# # Save results_hybrid to a CSV file
# with open('results_hybrid_single_tumor_Entropy.csv', 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerow(results_hybrid)
#
# # Save results_tumors_hybrid to a CSV file
# with open('results_tumors_hybrid_single_tumor_Entropy.csv', 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerows(results_tumors_hybrid)
#
# # Save results_hybrid_states to a CSV file
# with open('results_states_hybrid_single_tumor_Entropy.csv', 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile)
#     writer.writerows(results_hybrid_states)

print("Results saved to the corresponding csv files.")