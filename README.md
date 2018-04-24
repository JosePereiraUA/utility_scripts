# Utility Scripts

This repository serves as a backup for utility scripts created by José Pereira, on the scope of Protein Conformation Prediction.
Most scripts are deployed using a TKinter GUI. It was last updated on 20 April 2018.

## Atomizer 1.0

This script automates the process of downloading data from MONGO database as XTC trajectory files.
The XTC trajectory frames are then clustered based on RMSD, and, if requested, the plot depicting Energy values vs RMSD is displayed to the user for Cluster extraction selection. The selected clusters' representative structure is then retrieved and transformed into an all-atom model for MD analysis.

### Prerequisites and dependencies

The script requires the pre-installment of the following dependencies:
 * Tkinter (including tkMessageBox)
 * ttkThemes
 * MatPlotLib
 * MDTraj
 * PyMongo
 * ProtoSyn (Not pubicly avaliable)
 * Numpy
 * Gui (Avaliable from this repository)
 
 For the full usage of this script, the following files are required:
 * PDB structure file of the original protein for RMSD comparison
 * A Parallel Tempering Database at a MONGO server


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
