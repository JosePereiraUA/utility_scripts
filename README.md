# Utility Scripts

This repository serves as a backup for utility scripts created by José Pereira, on the scope of Protein Conformation Prediction.
Most scripts are deployed using a TKinter GUI. It was last updated on 26 April 2018. The following sequence of scripts facilitates the process of:
 1. Downloading PT trajectory data from Mongo;
 2. Run clusteriazation;
 3. Transform the selected coarse-grained structure to 'all'-atoms model;
 4. Run complete molecular dynamics simulations;
 5. Analyze the generated data.
***


## Atomizer 1.1

This script automates the process of downloading data from MONGO database as XTC trajectory files.
The XTC trajectory frames are then clustered based on RMSD, and, if requested, the plot depicting Energy values vs RMSD is displayed to the user for Cluster extraction selection. The selected clusters' representative structure is then retrieved and transformed into an all-atom model for MD analysis.

### Prerequisites and dependencies

The script requires the pre-installment of the following dependencies:
 * Tkinter (including tkMessageBox)
 * ttkThemes
 * MatPlotLib
 * MDTraj
 * PyMongo
 * ProtoSyn (Not publicly avaliable)
 * Numpy
 * Gui (Avaliable from this repository)
 
 For the full usage of this script, the following files are required:
 * PDB structure file of the original protein for RMSD comparison
 * A Parallel Tempering Database at a MONGO server
 
 
### Tool output
The script outputs the following files to the current working directory ('traj' is the default name):
 * __traj_all_atoms.pdb__ > Representative structure of the selected cluster in all_atoms form
 * __traj.log__           > Log file depicting all the clusters formed
 * __clusters.png__       > Figure of the RMSD vs Energy graph
 * __traj_clusters.pdb__  > Trajectory where the PDB frames are the representative structures of each cluster (without 'X's)
 * __traj_repr_str.pdb__  > Representative stucture of the selected cluster in caterpillar form (with 'X')
 * __others__             > Mostly leftover files, but can be usefull in certain cases: traj.xtc, traj.xvg, rmsd-dist.xvg, traj.xpm

![alt text](https://raw.githubusercontent.com/JosePereiraUA/utility_scripts/master/atomizer_py_schematics.png)
***

## MD Maker 1.0

This script automates the process of running and MD simulation.
for each simulation step, the user can select the number of steps (or simulation time by multiplying by dt = 0.002), aswell as the frequency of data storage.
The MD process is comprised of 7 automated steps:
 1. Preparation (Adds bounding box and creates topology files);
 2. Minimization 1 (Without solvent);
 3. Solvation (Adds spc216 water molecules);
 4. Minimization 2 (With solvent);
 5. NVT equilibration;
 6. NPT equilibration;
 7. Data collection.

### Prerequisites and dependencies

The script requires the pre-installment of the following dependencies:
 * Tkinter (including tkMessageBox)
 * ttkThemes
 * Gui (Avaliable from this repository)
 
 For the full usage of this script, the following files are required:
 * An .mdp storage folder (Avaliable from this repository)
 * A starting .pdb structure
***

## Analyzer 1.0

This script intends to automize the process of analysing GROMACS MD results. Based on user requests, it can perform up 4 tasks:
* Plot RMSD vs Frame 0 of the MD trajectory
* Plot RMSD vs Original protein
* Plot Energy potential as a function of simulation time
* Plot Gyration degrees (total and in discrete axis) as a function of simulation time

### Prerequisites and dependencies

The script requires the pre-installment of the following dependencies:
 * Tkinter (including tkMessageBox)
 * ttkThemes
 * MatPlotLib
 * Gui (Avaliable from this repository)
 
For the full usage of the script, the following GROMACS files are required:
 * TPR topology
 * XTC trajectory
 * EDR energy file
 * PDB structure file of the original protein for RMSD comparison
