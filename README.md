# Utility Scripts

This repository serves as a backup for utility scripts created by José Pereira.
Most scripts are deployed using a TKinter GUI.

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
 
For the full usage of the script, the following GROMACS files are required:
 * TPR topology
 * XTC trajectory
 * EDR energy file
 * PDB structure file of the original protein for RMSD comparison

