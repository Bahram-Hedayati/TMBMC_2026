# Bio-Inspired Navigation for Intelligent Nanoscale Drug Delivery Systems based on Chemotaxis and Entropy

This repository contains the implementation of a simulation framework for evaluating navigation strategies of a **Nanoscale Medical Agent (NMA)** in a **Tumor Microenvironment (TME)**.

The project is developed based on our research on hybrid navigation mechanisms combining **chemotaxis** and **information-theoretic (entropy-based)** decision-making.

---

## 📌 Overview

The NMA operates as an autonomous agent within a discretized TME and aims to locate and treat tumor regions by sensing oxygen concentration (hypoxia).

The framework evaluates navigation performance under different environmental conditions and strategies.

---

## 🧪 Scenarios

The project includes two types of tumor environments:

1. **Single-Tumor TME**
2. **Multi-Tumor TME**

---

## 🚀 Navigation Mechanisms

Each scenario is evaluated using three navigation strategies:

1. **Hybrid (Chemotaxis + Entropy)**
2. **Chemotaxis-based**
3. **Random Walk**

---

## 🧠 Key Concept

- **Chemotaxis**: Moves toward lower oxygen concentration (hypoxic regions)
- **Entropy-based navigation**: Explores regions with higher uncertainty
- **Hybrid mechanism**: Dynamically switches between exploration and exploitation

---

## 📂 Project Structure

### Core Modules

- **`sysEnvClasses.py`**  
  Defines the main system entities (NMA and TME) using object-oriented design.

- **`publicFns.py`**  
  Utility functions for visualization and reporting.

---

### Environment & Initialization

- **`Tumor_Development.py`**  
  Loads and visualizes tumor masks (single and multi-tumor).

- **`TME_init.py`**  
  Initializes the lattice-based TME environments.

---

### Single-Tumor Simulations

- **`SingleTumor_Hybrid.py`**  
- **`SingleTumor_Chemotaxis.py`**  
- **`SingleTumor_Random.py`**  

Each script runs ensemble simulations for the respective navigation strategy.

---

### Multi-Tumor Simulations

- **`MultipleTumor_Hybrid.py`**  
- **`MultipleTumor_Chemotaxis.py`**  
- **`MultipleTumor_Random.py`**

---

### Visualization

- **`Visual_Ensemble.py`**  
  Generates plots and evaluation metrics across simulation runs.

---

## 📊 Data

The `Data/` folder contains:

- Tumor masks (single & multi-tumor)
- Steady-state oxygen distributions

All data are stored in CSV format and are used to initialize the simulation environments.

---

## ⚙️ Simulation Details

The simulation framework is designed as a lattice-based agent-based model (ABM) to capture the interaction between a nanoscale medical agent (NMA) and the tumor microenvironment (TME).

### Environment Representation
- The TME is discretized into a **2D lattice of voxels**, where each voxel represents a biologically meaningful spatial unit.
- Each voxel stores:
  - Oxygen concentration (hypoxia level)
  - Tumor status (healthy/cancerous)

### Tumor Modeling
- Tumor growth is generated using an **Invasion Percolation (IP)** model.
- Two configurations are considered:
  - Single-tumor (centralized growth)
  - Multi-tumor (spatially distributed tumors with comparable total size)

### Oxygen Diffusion
- Oxygen distribution is obtained by solving a **reaction–diffusion model**:
  - Diffusion across the lattice
  - Uptake by cancerous voxels
- The system is evolved until **steady state**, which is verified via RMSE convergence analysis.
- The resulting steady-state oxygen field is used as input for navigation.

### Time Scales
- Two separate time scales are used:
  - **Diffusion time step**: fine-grained (for numerical stability)
  - **Movement time step**: coarser (for agent navigation)
- This separation improves computational efficiency while preserving physical realism.

### Agent Model (NMA)
- The NMA is modeled as a **single autonomous agent** operating on the lattice.
- At each time step, the NMA:
  1. Senses local oxygen values
  2. Forms a probability distribution over candidate directions
  3. Selects the next movement direction based on the navigation policy

### Navigation Strategies
- **Chemotaxis**: biased movement toward hypoxic regions
- **Entropy-based**: exploration based on local uncertainty
- **Hybrid**: adaptive switching between the two strategies

### Entropy Computation
- For each direction, a **local 3×3 voxel set** is evaluated
- Oxygen values are transformed into weights and normalized to form a probability distribution
- **Shannon entropy** is computed to quantify uncertainty

### Evaluation Protocol
- Each experiment is repeated over **1000 ensemble runs**
- Metrics include:
  - Tumor detection effectiveness
  - Number of treated voxels
  - Navigation dynamics and trajectories
  - State-switching behavior (hybrid method)

## 👥 Contributors
- Bahram Hedayati — Implementation and research
- Dr. Karthik Reddy Gorla — Scientific guidance and supervision  

## 🎓 Academic Supervision
- Prof. Maurizio Magarini - Supervisor
  