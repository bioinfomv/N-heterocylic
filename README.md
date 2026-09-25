# N-heterocylic
Study of N-heterocylic
# MD Simulation Analysis

A collection of scripts for analyzing molecular dynamics (MD) simulation trajectories and characterizing protein–ligand interactions, conformational states, and structural dynamics.

### Analyses Included

* **PCA:** Principal component analysis of MD trajectories to characterize major conformational motions and compare conformational distributions.
* **Clustering:** Clustering of MD trajectory conformations to identify representative structural states.
* **Aromatic Interaction Analysis:** Calculation of aromatic ring centroid distances and interaction occupancy between protein and ligand aromatic rings.
* **Hydrophobic Contact Analysis:** Calculation of protein hydrophobic residue–ligand contact occupancy during the MD trajectory.
* **Scaffold PCA:** PCA based on selected scaffold/contact-associated protein Cα residues to compare structural conformations across systems.

### Requirements

```text
Python
MDTraj
MDAnalysis
NumPy
Pandas
scikit-learn
Matplotlib
PyEMMA
RDKit
memory_profiler
VMD
```

### Output

The scripts generate trajectory analysis files containing:

* Contact occupancies
* Aromatic interaction distances and occupancies
* Cluster assignments/representative structures
* PCA projections (PC1/PC2)
* PCA plots for comparison of different systems

The scripts are designed to be adapted by changing the input files, selections, and analysis parameters for different MD systems.
