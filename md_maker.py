#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Tkinter import *
from ttk import *
from ttkthemes import ThemedStyle
from ttkthemes import themed_tk as tk
import tkMessageBox as alert
from gui import *
import os

class Default(Step):
    """
    Holds default values for various variables.
    """

    #Application window settings
    TITLE          = "MD Maker"
    FOOTER         = "José Pereira @ 2018"
    VERSION        = 1.1
    RESIZABLE      = False

    #Input defaults
    INPUT_DIR      = os.getcwd()
    MDP_DIR        = "/home/jpereira/Desktop/utility_scripts/mdps"
    INPUT_FILE     = "input.pdb"
    ALL_H_FILE     = INPUT_FILE[:-4] + '_all_h.pdb'
    BOXED_FILE     = INPUT_FILE[:-4] + '_boxed.pdb'
    SOLVATED_FILE  = INPUT_FILE[:-4] + '_solvt.pdb'
    IONIZED_TPR    = INPUT_FILE[:-4] + '_ion.tpr'
    IONIZED_FILE   = INPUT_FILE[:-4] + '_ion.pdb'
    BOX_SHAPE      = "cubic"
    BOX_SIZE       = 4
    ION_COUNT      = 2
    DEBUG_MODE     = True
    STEP_TYPES     = ['minimization', 'equilibration', 'collection']
    N_STEPS        = {'minimization': 50000, 'equilibration': 500000, 'collection': 250000}
    PRINT_EVERY    = {'minimization': 1000, 'equilibration': 2000, 'collection': 20000}
    WATER_TYPE     = "spc216"
    DT             = 0.002
    STEPS          = [
                        Step('minimization_1'),
                        Step('solvate'),
                        Step('minimization_2'),
                        Step('equilibration_nvt'),
                        Step('equilibration_npt'),
                        Step('collection')
                    ]

    #Application theme settings
    PAD            = 10
    LABEL_WIDTH    = 170
    THEME          = "arc"
    BACKGROUND     = "white"
    HOVER_COLOR    = "steel blue"
    FONT_SIZE      = "10"
    FONT           = "System"


class Utils(BasicGUIFunctions):

    def __init__(self, master):
        """
        Creates a new application window with corrects settings and GUI layout.
        """
        
        #Set application's window settings and title
        self.config_master(master, Default.TITLE, Default.VERSION, Default.RESIZABLE, Default.BACKGROUND)  
        self.n_steps     = {'minimization': None, 'equilibration': None, 'collection': None}
        self.print_every = {'minimization': None, 'equilibration': None, 'collection': None}

        self.input     = SelectFile(master, "Input PDB:", default_value = Default.INPUT_DIR + '/' + Default.INPUT_FILE,
            row = 0, label_width = Default.LABEL_WIDTH)
        self.mdps      = SelectFile(master, "MDPs directory:", default_value = Default.MDP_DIR,
            row = 1, label_width = Default.LABEL_WIDTH, search_directory_only = True)
        
        master_frame_1 = Frame(master, background = Default.BACKGROUND)
        master_frame_1.grid(row = 2, sticky = 'w', pady = Default.PAD)
        self.box_size  = Input(master_frame_1, "Box size (nm):", Default.BOX_SIZE, label_width = Default.LABEL_WIDTH, width = 5)
        self.ion_count = Input(master_frame_1, "Number of ions:", Default.ION_COUNT, label_width = Default.LABEL_WIDTH, width = 5, column = 1)

        master_frame_2 = Frame(master, background = Default.BACKGROUND)
        master_frame_2.grid(row = 3)
        for column, step in enumerate(Default.STEP_TYPES):
            frame = Container(master_frame_2, default_text = " {step} settings".format(step = step.title()), column = column)
            self.n_steps[step]     = Input(frame(), "N steps:", Default.N_STEPS[step], label_width = 120, width = 7)
            self.print_every[step] = Input(frame(), "Print every:", Default.PRINT_EVERY[step], row = 1, label_width = 120, width = 7)
        
        self.debug_mode = CheckBox(master, "Debug Mode", default_status = Default.DEBUG_MODE, row = 4)
        
        self.add_basic_utilities(master, 5, pad = Default.PAD, background = Default.BACKGROUND,
            hover_color = Default.HOVER_COLOR, footer = Default.FOOTER)


    def process(self):
        """
        Process control. Generate data.
        """

        #1. Set progress bar and debug settings
        self.progress_bar['maximum'] = self.to_do()
        self.progress_bar['value'] = 0
        self.debug = '' if self.debug_mode() else '2>&1'

        #2. Verify input
        self.verify_input()
        
        self.create_topology_files(self.input())

        self.add_bounding_box(self.box_size())

        self.apply_step_settings()
        previous_step_file = Default.BOXED_FILE
        for step in Default.STEPS:
            if step.title == 'solvate':
                previous_step_file = self.solvate(previous_step_file)
            else:
                previous_step_file = self.md_run(step, previous_step_file)

        self.clean()

    
    def solvate(self, input_file, water_type = Default.WATER_TYPE, solvated_file = Default.SOLVATED_FILE, 
        tpr_file = Default.IONIZED_TPR, output_file = Default.IONIZED_FILE):
        """
        Adds solvent molecules to the simulation box,
        correctly labeling the added molecules in the topology file.
        """

        self.update_terminal("Adding {water_type} water molecules to simulation ..."\
            .format(water_type = water_type))
        self.run("gmx solvate -cp {input_file} -cs {water_type}.gro -o {output_file} -p topol.top"\
                .format(input_file = input_file, water_type = water_type, output_file = solvated_file))
        
        self.update_terminal("Adding ions to system ...")
        self.run("gmx grompp -f {mdp} -c {input_file} -p topol.top -o {output_file}"\
            .format(mdp = self.mdps() + '/ions.mdp', input_file = solvated_file, output_file = tpr_file),
            update_progress_bar = False)

        self.run("echo '13' | gmx genion -s {input_file} -p topol.top -o {output_file} -pname NA -np {ion_count}"\
            .format(input_file = tpr_file, output_file = output_file, ion_count = self.ion_count()))

        return output_file


    def md_run(self, step, input_file):
        """
        Runs an MD simulation and stores the produced data.
        """
                
        #1.
        self.update_terminal("Configuring {mdp} and pre-processing topology ...".format(mdp = step.mdp))
        step.config_mdp(self.mdps())
        self.run("gmx grompp -f {mdp} -c {input_file} -p topol.top -o {output_file} -maxwarn 1"\
            .format(mdp = step.mdp, input_file = input_file, output_file = step.tpr))

        #2.
        self.update_terminal("Running {title} ...".format(title = step.title))
        self.run("gmx mdrun -deffnm {title} -nt 1 -s {tpr} -v"\
            .format(title = step.title, tpr = step.tpr))
        
        #3.
        self.update_terminal("Extracting trajectory's last frame to {title}.pdb ...".format(title = step.title))
        last_frame = self.verify_output(step)
        self.run("echo '0 \n 0' | gmx trjconv -f {title}.trr -o {title}.pdb -s {title}.tpr -pbc mol -conect -center -b {frame}"\
            .format(title = step.title, frame = last_frame))

        #4.
        self.update_terminal("Extracting trajectory to {title}_traj.pdb ...".format(title = step.title))
        self.run("echo '0 \n 0' | gmx trjconv -f {title}.trr -o {title}_traj.pdb -s {title}.tpr -pbc mol -conect -center"\
            .format(title = step.title))

        #5.
        self.update_terminal("Storing data in {folder} ...".format(folder = '/' + step.title))
        self.store_data(step.title)

        return step.title + '/' + step.title + '.pdb'


    def apply_step_settings(self, steps = Default.STEPS, step_types = Default.STEP_TYPES):
        """
        Applies input settings regarding number of steps and frequncy of data storage to each step.
        """

        for step_type in step_types:
            for step in steps:
                if step.title.startswith(step_type):
                    step.n_steps     = self.n_steps[step_type]()
                    print step.n_steps
                    print type(step.n_steps)
                    step.print_every = self.print_every[step_type]()


    def store_data(self, title):
        """
        Creates a new folder and moves all necessary files for correct data keeping.
        """

        self.run("mkdir {title}".format(title = title), update_progress_bar = False)
        self.run("mv {title}* mdout.mdp {title}".format(title = title))


    def create_topology_files(self, input_file, output_file = Default.ALL_H_FILE):
        """
        1) Creates a topology file using GMX PDB2GMX;
        2) Edits the produced topology file to include solvent and ions forcefield.
        """

        self.update_terminal("Creating topology files ...")
        self.run("echo '6 \n 7' | gmx pdb2gmx -f {input_file} -o {output_file}"\
            .format(input_file = input_file, output_file = output_file))
        
        self.update_terminal("Adding solvent data to topology file ...")
        lines = []
        with open('topol.top', 'r') as file_in:
            for line in file_in:
                elem = line.split()
                if len(elem) > 3 and elem[3] == 'restraint':
                    lines.append('; Include water topology\n#include "amber99sb-ildn.ff/tip3p.itp"\n\n; Include topology for ions\n#include "amber99sb-ildn.ff/ions.itp"\n\n')
                    lines.append(line)
                else:
                    lines.append(line)
        with open('topol.top', 'w') as file_out:
            for line in lines:
                file_out.write(line)
        self.update_progress_bar()

        #self.update_terminal("Cleaning unnecessary files ...")
        #self.run("rm -rf posre.itp")


    def add_bounding_box(self, size, input_file = Default.ALL_H_FILE,
        output_file = Default.BOXED_FILE, shape = Default.BOX_SHAPE):
        """
        Adds a bounding box to the input structure, of the requested size, using GMX EDITCONF
        """
        
        self.update_terminal("Boxing protein with {shape} box of size {size} ..."\
            .format(shape = shape, size = size))
        self.run("gmx editconf -f {input_file} -o {output_file} -c -box {size} -bt {shape}"\
            .format(input_file = input_file, output_file = output_file, size = size, shape = shape))


    def verify_output(self, step):
        """
        1) Verifies if the minimization steps converged in the input number of steps;
        2) Returns the frame number/time of the last strucuture in the trajectory.
        """

        with open(step.title + '.log', 'r') as file_in:
            for line in file_in:
                if step.title.startswith('minimization'):
                    if line.startswith("Steepest Descents did not converge"):
                        self.error("Convergence Error", line, auto_exit = False)
                        self.store_data(title)
                        self.clean()
                        exit(1)
                    elif line.startswith("Steepest Descents converged"):
                        elems = line.split()
                        return elems[7]
                else:
                    return float(step.n_steps) * Default.DT


    def verify_input(self):
        """
        1) Verifies the existance of the necessary .mdp files;
        2) Verifies the existance of previous MD data in the current work directory.
        """

        for step in Default.STEPS:
            if step.title == 'solvate':
                continue
            file_name = self.mdps() + '/' + step.title + '.mdp'
            if not os.path.isfile(file_name):
                self.error("Missing File", "{file} not found".format(file = file_name))

        for step in Default.STEPS:
            if os.path.isdir(os.getcwd() + '/' + step.title):
                result = alert.askquestion("Please confirm directory deletion",
                    "One or more directories with previous MD run data have been found. Delete and replace them?",
                    icon='warning')
                if result == "yes":
                    self.clean_directories()
                    break
                else:
                    exit(1)


    def clean_directories(self):
        """
        Delete pre-existent folders of previous MD data.
        """

        for step in Default.STEPS:
            self.run("rm -rf {title}".format(title = step.title), update_progress_bar = False)


    def to_do(self):
        """
        Count the total number of tasks requested to proper set the progress bar.
        """

        to_do  = 1                            #clean
        to_do += 2                            #create_topology_files
        to_do += 1                            #add_bounding_box
        to_do += 5 * (len(Default.STEPS) - 1) #MD Steps
        to_do += 2                            #Solvate

        return to_do
    

#*****************************************
#*                 MAIN                  *
#*****************************************

#1. Instantiate window and apply themed style
root = tk.ThemedTk()
root.set_theme(Default.THEME)

#2. Instantiate GUI and run main loop
app = Utils(root)
root.mainloop()