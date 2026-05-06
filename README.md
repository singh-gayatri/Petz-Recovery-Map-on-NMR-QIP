# Realizing the Petz Recovery Map on an NMR Quantum Processor
This repository contains the experimental tomography data for and Quantum state tomography MATLAB code for the implementation of Petz Recovery Maps on a three-qubit Nuclear Magnetic Resonance (NMR) quantum information processor. The full theoretical framework and experimental results are available in "[https://journals.aps.org/pra/accepted/10.1103/xd6k-swv7](https://journals.aps.org/pra/abstract/10.1103/xd6k-swv7)"


## Research Overview: 
We demonstrated experimental realization of damping channel and associated Petz recovery maps using duality quantum computing (DQC) algorithm. 
#### Noise Models & Reference States
We investigate two paradigmatic single-qubit noise models with parametrized reference state $\sigma$ for Petz map: 
1. **Amplitude Damping (AD) channel**: $\sigma=(1-\epsilon)\ket{0}\bra{0}+\epsilon\ket{1}\bra{1}$
2. **Phase Damping (PD) channel**: $\sigma=(1-\epsilon)\ket{+}\bra{+}+\epsilon\ket{-}\bra{-}$
The Petz recovery operation was implemented for reference states parameter $\epsilon \in \\{0.2,0.5,0.8\\}$.
#### Experimental Input States
 1. **For AD channel** : $\ket{0},\ket{1},\ket{+}=\frac{\ket{0}+\ket{1}}{\sqrt{2}}$ and $0.9268 \ket{0}+0.3754\imath \ket{1}$.
 2. **For PD channel**: $\ket{\pm}=\frac{\ket{0}\pm\ket{1}}{\sqrt{2}},\ket{0}$ and $0.9268 \ket{0} + 0.3754\imath \ket{1}$. 
## Data Organization
The experimental data is organized into two primary folders: **AD** and **PD**. Each directory contains subfolders named after the specific input state used, such as **zero, one, plus, minus**, or **extra**.

Each data file contains data points representing the full experimental sequence: 
1. Index 0: Pseudo-Pure State as a reference.
2. Indices 1–22: Data for the Damping Channel only (before recovery).
3. Indices 23–44: Petz recovery map with $\epsilon = 0.2$.
4. Indices 45–66: Petz recovery map with $\epsilon = 0.5$.
5. Indices 67–88: Petz recovery map with $\epsilon = 0.8$.

The **Tomography_code** folder contains the main MATLAB script and necessary helper functions for performing Quantum State Tomography (QST) and calculating fidelity with respect to the theoretical state:

## How to Use
To reconstruct a density matrix and calculate the QST fidelity:
1. Navigate to the Tomography_code folder.
2. Within the script, define the basePath (currently set to "Petz recovery map\final data\") where your primary data folders are stored.
3. Execute the main MATLAB script.
4. The script is set to process a default number of files (nFiles=4, can be modified). On execution, it create an output text file **results.txt** containing the mean and standard deviation of the fidelity in the following format:

```
Channel: [Name] | State: [Name]
--------------------------------------------------------------------------------
damping strength |                Petz recovery map with ϵ              
      p          |        0.2               0.5               0.8        
--------------------------------------------------------------------------------
    [data]              [data]            [data]            [data]
```    
The reconstruction follows the convex optimization method described by Gaikwad et al., QIP 20, 19 (2021).

