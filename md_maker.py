#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tkinter import *
# from ttk import *
from ttkthemes import ThemedStyle
from ttkthemes import themed_tk as tk
from tkinter import messagebox as alert
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
    INPUT_FILE     = "md.pdb"
    INPUT_ITP      = "ligand.itp"
    ALL_H_FILE     = INPUT_FILE[:-4] + '_all_h.pdb'
    BOXED_FILE     = INPUT_FILE[:-4] + '_boxed.pdb'
    SOLVATED_FILE  = INPUT_FILE[:-4] + '_solvt.pdb'
    IONIZED_TPR    = INPUT_FILE[:-4] + '_ion.tpr'
    IONIZED_FILE   = INPUT_FILE[:-4] + '_ion.pdb'
    BOX_SHAPE      = "cubic"
    BOX_SIZE       = 7
    ION_COUNT      = 0
    DEBUG_MODE     = True
    LIGAND_ONLY    = False
    STEP_TYPES     = ['minimization', 'equilibration', 'collection']
    N_STEPS        = {'minimization': 50000, 'equilibration': 500000, 'collection': 10000000}
    PRINT_EVERY    = {'minimization': 1000, 'equilibration': 2000, 'collection': 20000}
    TIME_STEP      = {'minimization': 0.0, 'equilibration': 0.002, 'collection': 0.002}
    WATER_TYPE     = "spc216"
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
        self.n_steps         = {'minimization': None, 'equilibration': None, 'collection': None}
        self.print_every     = {'minimization': None, 'equilibration': None, 'collection': None}
        self.time_step       = {'minimization': None, 'equilibration': None, 'collection': None}
        self.simulation_time = {'minimization': None, 'equilibration': None, 'collection': None}

        self.input     = SelectFile(master, "Input PDB:", default_value = Default.INPUT_DIR + '/' + Default.INPUT_FILE,
            row = 0, label_width = Default.LABEL_WIDTH)
        self.itp     = SelectFile(master, "Add Ligand ITP:", default_value = Default.INPUT_DIR + '/' + Default.INPUT_ITP,
            row = 1, label_width = Default.LABEL_WIDTH)
        self.mdps      = SelectFile(master, "MDPs directory:", default_value = Default.MDP_DIR,
            row = 2, label_width = Default.LABEL_WIDTH, search_directory_only = True)
        
        master_frame_1 = Frame(master, background = Default.BACKGROUND)
        master_frame_1.grid(row = 3, sticky = 'w', pady = Default.PAD)
        self.box_size  = Input(master_frame_1, "Box size (nm):", Default.BOX_SIZE, label_width = Default.LABEL_WIDTH, width = 5)
        self.ion_count = Input(master_frame_1, "Number of ions:", Default.ION_COUNT, label_width = Default.LABEL_WIDTH, width = 5, column = 1)

        master_frame_2 = Frame(master, background = Default.BACKGROUND)
        master_frame_2.grid(row = 4)
        for column, step in enumerate(Default.STEP_TYPES):
            frame = Container(master_frame_2, default_text = " {step} settings".format(step = step.title()), column = column)
            self.n_steps[step]       = Input(frame(), "N steps:", Default.N_STEPS[step], label_width = 120, width = 7)
            self.print_every[step]   = Input(frame(), "Print every:", Default.PRINT_EVERY[step], row = 1, label_width = 120, width = 7)
            if Default.TIME_STEP[step] != 0.0:
                self.time_step[step] = Input(frame(), "Time Step (ps):", Default.TIME_STEP[step], row = 2, label_width = 120, width = 7)
                footer = Label(frame(), text = "Simulation time: " + str((float(self.n_steps[step]()) * float(self.time_step[step]())) * 0.001) + " ns",
                    foreground = Default.HOVER_COLOR, background = Default.BACKGROUND)
                footer.grid(row = 3)
        
        self.debug_mode = CheckBox(master, "Debug Mode", default_status = Default.DEBUG_MODE, row = 5)
        self.ligand_only = CheckBox(master, "Ligand Only", default_status = Default.LIGAND_ONLY, row = 6)
        
        self.add_basic_utilities(master, 7, default_text = "Start_MD", pad = Default.PAD,
            background = Default.BACKGROUND, hover_color = Default.HOVER_COLOR, footer = Default.FOOTER)


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

        if self.ligand_only():
            self.run("echo '4' | gmx genion -s {input_file} -p topol.top -o {output_file} -pname NA -np {ion_count}"\
                .format(input_file = tpr_file, output_file = output_file, ion_count = self.ion_count()))
        else:
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
                    step.n_steps       = self.n_steps[step_type]()
                    step.print_every   = self.print_every[step_type]()
                    if self.time_step[step_type] != None:
                        step.time_step = self.time_step[step_type]()


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

        def extract_topology_section(file_handler):
            lines = []
            title = ""
            while True:
                line = file_handler.readline()
                elem = line.split()
                if len(elem) == 3 and elem[0] == "[" and elem[2] == "]":
                    title = elem[1]
                lines.append(line)
                if len(elem) == 0:
                    break
            return title, lines

        self.update_terminal("Creating topology files ...")
        # 1) Read the input ITP (if it exists) and gather information about the
        # ligand
        to_extract   = False
        residue_name = None
        itp_exists   = False
        itp_atoms    = []

        # The given ITP file should be the output of the tleap + acpype method.
        # Therefore, two parts are taken from this file, the ligand itp and the
        # required atomtypes. These will be stored in the correct files and 
        # later "included" into the final topology.
        new_itp, new_atomtypes = [], []
        s = ["moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals"]
        if os.path.exists(self.itp()) and not self.ligand_only():
            itp_exists = True
            itp = open(self.itp(), "r")
            while True:
                next_title, next_section = extract_topology_section(itp)
                if next_title == "atomtypes":
                    new_atomtypes += next_section
                if next_title in s:
                    new_itp += next_section

                # Extract atom records from the atom section
                if next_title == "atoms":
                    for line in next_section:
                        if not line.startswith("[") or not line.startswith(";"):
                            elem = line.split()
                            if len(elem) < 4 or elem[0] == ";":
                                continue
                            itp_atoms.append(elem[4])
                            if residue_name == None:
                                residue_name = elem[3]
                            if elem[3] != residue_name:
                                print(elem[3], "!=", residue_name)
                                print("More than 1 ligand molecule found in itp!")
                                exit(0)
                if next_title == "" and len(next_section) == 1:
                    break

            # Write atomtypes file for later to be "included" into the final
            # topology
            with open("atomtypes.itp", "w") as atomtypes:
                for line in new_atomtypes:
                    atomtypes.write(line)

            itp.close()
        else:
            print(" Defined itp not found: continuing ...\n")
        
        if itp_exists:

            # 2) Read PDB and try to identify any ligand that matches the given
            # itp, by residue name
            pdb_atoms, protein_lines, ligand_lines = [], [], []
            with open(input_file) as pdb:
                for line in pdb:
                    elem = line.split()
                    if len(elem) < 4:
                        continue
                    if elem[3] == residue_name:
                        pdb_atoms.append(elem[2])
                        ligand_lines.append(line)
                    elif elem[0] == "ATOM":
                        protein_lines.append(line)

            # Note: itp_atoms and pdb_atoms both refer to the ligand atoms.
            # itp_atoms are records extracted from the provided itp while
            # pdb_atoms are records from the input PDB. These records are only
            # the atom names. Therefore, in order for a ligand in the PDB to be
            # recognized as the ligand in the itp file, the files must show:
            #  . Same number of ligand atoms
            #  . Same order of ligand atoms
            #  . Same atom names
            # On the other hand, protein_lines and ligand_lines contain
            # information from the PDB complete line: only of all other system
            # atoms (normally a protein) that are not part of the ligand, in the
            # case of protein_lines, and all ligand atoms in the case of
            # ligand_lines
            if itp_atoms != pdb_atoms:
                print("Ligand itp atoms do not match the protein ligand atoms!")
                exit(0)

            # 3) Create correct ITP file that will later be included in the
            # final topology
            with open(self.itp()[:-4] + "_auto.itp", "w") as new_itp_fh:
                for line in new_itp:
                    new_itp_fh.write(line)

            # 4) Create temporary PDB of just protein for pdb2gmx topology
            input_name = "%s_auto_no_ligand.pdb" % (input_file[:-4])
            with open(input_name, "w") as tmp:
                for line in protein_lines:
                    tmp.write(line)

            # 5) pdb2gmx on protein only
            output_name = "%s_auto_no_ligand_corrected.pdb" % (input_file[:-4])
            self.run("echo '6 \n 7 \n 0' | gmx pdb2gmx -f {input_file} -o {output_file} -his"\
                .format(input_file = input_name, output_file = output_name))
            # Note: In the above command, histidines are automatically set to be
            # of type HID. These are the possible histidine types:
            # 0) H on HD1  -> HID
            # 1) H on NE2  -> HIE
            # 2) H on both -> HIP
            # 3) No H      -> HIS (?)
            # You should manually check this and change the above line!
            # Todo: Automatically check the histodine protonation state
            
            # 6) Add initial ligand information to the corrected_protein_only
            # PDB file, with correct numeration
            joined_pdb = []
            last_line  = None
            with open(output_name, "r") as output_name_in:
                for line in output_name_in:
                    if line == "TER\n":
                        for ligand_line in ligand_lines:
                            joined_pdb.append(ligand_line)
                    joined_pdb.append(line)
                    last_line = line

            with open(output_file, "w") as joined_out:
                for line in joined_pdb:
                    joined_out.write(line)                    

        elif not self.ligand_only():
            self.run("echo '6 \n 7 \n 0' | gmx pdb2gmx -f {input_file} -o {output_file} -his"\
                .format(input_file = input_file, output_file = output_file))

        if self.ligand_only():
            self.run("cp %s topol.top" % (self.itp()))
            self.run("cp %s %s" % (input_file, output_file))

        self.update_terminal("Adding solvent and ligand data to topology file ...")
        lines = []
        with open('topol.top', 'rb') as file_in:
            for line in file_in:
                try:
                    elem = [e.decode("utf-8") for e in line.split()]
                except:
                    title = os.path.basename(input_file[:-4])
                    lines.append("%s\n".encode() % (title).capitalize().encode())
                    continue
                if len(elem) == 3 and elem[1] == "atomtypes":
                    if self.ligand_only():
                        lines.append('; Include forcefield parameters\n#include "amber99sb-ildn.ff/forcefield.itp"\n\n'.encode())
                    lines.append(line)
                elif len(elem) == 2 and elem[1][-15:-1] == "forcefield.itp":
                    lines.append(line)
                    if itp_exists:
                        lines.append('#include "atomtypes.itp"\n'.encode())
                elif len(elem) == 3 and elem[1] == 'system':
                    lines.append('; Include water topology\n#include "amber99sb-ildn.ff/tip3p.itp"\n\n; Include topology for ions\n#include "amber99sb-ildn.ff/ions.itp"'.encode())
                    if itp_exists:
                        lines.append('\n\n; Include ligand topology\n#include "%s"'.encode() % (self.itp()[:-4] + "_auto.itp").encode())
                    lines.append("\n\n".encode())
                    lines.append(line)
                elif len(elem) == 2 and elem[1] == '1':
                    lines.append(line)
                    if itp_exists:
                        lines.append("%-19s 1\n".encode() % (residue_name).encode())
                else:
                    lines.append(line)

        with open('topol.top', 'wb') as file_out:
            for line in lines:
                file_out.write(line)
        self.update_progress_bar()

        if itp_exists:
            self.run("cp %s %s" % (input_file, output_file))
            self.run("rm -rf %s %s" % (input_name, output_name))

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
                        if elems[4] == "machine":
                            return elems[7]
                        else:
                            return elems[8]
                else:
                    return float(step.n_steps) * float(step.time_step)


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
