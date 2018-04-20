#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Tkinter import *
from ttk import *
from ttkthemes import ThemedStyle
from ttkthemes import themed_tk as tk
from gui import *
import os
from os import popen as instantiate
import tkMessageBox as alert
from matplotlib import pyplot as plt

class Default:
    """
    Holds default values for various variables.
    """

    #Application window settings
    TITLE         = "MD Analyzer"
    THEME         = "arc"
    BACKGROUND    = "white"
    VERSION       = 1.0
    RESIZABLE     = False
    
    #Application runtime settings
    PAD           = 10
    LABEL_SIZE    = 12
    PRINT_EVERY   = 100
    DEBUG_MODE    = False
    ORIGINAL_DIR  = os.getcwd()
    ORIGINAL_FILE = '1ctf.pdb'
    HOVER_COLOR   = "steel blue"
    PLOT_COLOR    = ["#4682b4", "#D74B4B", "#F39D41", "#77a047"]
    FOOTER        = "José Pereira @ 2018"


class Utils:

    def __init__(self, master):
        """
        Creates a new application window with corrects settings and GUI layout.
        Directory and file selection defaults to current directory.
        """

        #Set application's window settings and title
        master.configure(background = Default.BACKGROUND)
        master.title("{title} v{version}".format(title = Default.TITLE, version = Default.VERSION))
        if not Default.RESIZABLE: master.resizable(0, 0)
        self.master = master

        #Set GUI elements
        self.dir           = SelectFile(master, "Select a working directory:", search_directory_only = True, default_value = os.getcwd(), row = 0)
        self.file_name     = Input(master, "Select the file prefix:", "md", row = 1)
        self.original      = SelectFile(master, "Original protein for RMSD comparison:", default_value = Default.ORIGINAL_DIR + '/' + Default.ORIGINAL_FILE, row = 2)
        self.skip          = Input(master, "Print every:", Default.PRINT_EVERY, row = 3, width = 5)
        Title(master, "Choose the analysis to perform:", row = 4, sticky = 'w')
        frame = Frame(master, background = Default.BACKGROUND)
        frame.grid(row = 5, sticky = 'w')
        self.rmsd_frame_0  = CheckBox(frame, "RMSD vs Frame 0",  column = 0)
        self.rmsd_original = CheckBox(frame, "RMSD vs Original", column = 1)
        self.potential_ene = CheckBox(frame, "Potential energy", column = 2)
        self.gyration_degr = CheckBox(frame, "Gyration Degrees", column = 3)
        self.debug_mode    = CheckBox(frame, "Debug Mode", row = 6, column = 0, default_status = Default.DEBUG_MODE)
        Button(master, text="Analyze (Necessary files: .trp, .xtc, .edr)", width = 100, command= self.analyze,
            activebackground = Default.HOVER_COLOR, relief = 'flat', background = Default.BACKGROUND)\
            .grid(row = 7, columnspan = 5, padx = Default.PAD)
        self.progress_bar  = Progressbar(master, orient = 'horizontal', mode = 'determinate', length = 700)
        self.progress_bar.grid(columnspan = 8, row = 8, pady = Default.PAD)
        self.terminal      = Terminal(master, row = 9, columnspan = 5)
        self.footer        = Footer(master, Default.FOOTER, row = 10, columnspan = 5)


    def analyze(self):
        """
        Generate data based on query requests. Analyze generated files and plot data.
        """
        
        #1. Set progress bar with the maximum amount of query requests
        self.progress_bar['maximum'] = self.to_do() + 1
        self.progress_bar['value'] = 0
        
        #2. If query requests exist, generate and plot data
        if self.progress_bar['maximum'] > 1:
            #2.1 Generate files using GROMACS
            self.generate_data()
            #2.2 Read generated files and plot data using MATPLOTLIB
            self.analyse_data()
            #2.3 Delete unnecessary files
            self.update_terminal("Cleaning unnecessary files ...")
            self.clean()
        #3. Finish process
        self.update_progress_bar()
        self.update_terminal("All tasks successefully performed.")


    def clean(self):
        """
        Deletes unnecessary files.
        """

        self.update_progress_bar()
        self.run("rm -rf *#")


    def check_file(self, file_name):
        """
        Checks the existance of a file in the requested directory.
        Raises a error and exits if it does not exist.
        """

        self.update_terminal("Checking for the existance of file {file} ..."\
            .format(file = file_name))
        if not os.path.isfile(file_name):
            alert.showerror("Missing File", "{file} not found".format(file = file_name))
            exit(1)
        self.update_progress_bar()


    def update_progress_bar(self):
        """
        Increases the completness of the progress bar by one step.
        """

        self.progress_bar['value'] = self.progress_bar['value'] + 1
        root.update_idletasks()


    def update_terminal(self, message):
        """
        Shows a message in the GUI terminal.
        """

        self.terminal(message)
        root.update_idletasks()


    def to_do(self):
        """
        Counts the total amount of tasks to be completed, given the current query request.
        A step usually comprises 3 tasks:
        1) Check file existance;
        2) Create necessary files using GROMACS;
        3) Read data and plot it using MATPLOTLIB;
        """

        to_do = 0
        if self.rmsd_frame_0(): to_do += 6
        if self.rmsd_original():
            if not self.rmsd_frame_0():
                #Note: If the trajectory pdb wasn't written, this task is passed onto this step.
                to_do += 3
            to_do += 3
        if self.potential_ene(): to_do += 3
        if self.gyration_degr(): to_do += 3
        return to_do


    def run(self, command):
        """
        Instantiate a new terminal command (GROMACS)
        """

        instantiate(command)
        self.update_progress_bar()


    def generate_data(self):
        """
        Verify the existance of necessary files.
        Run GROMACS instances to convert and produce requested files.
        """
        
        #1. Define necessary variables
        skip = int(self.skip())
        file_name = self.dir() + '/' + self.file_name()
        debug = '' if self.debug_mode() else '2>&1'

        #2. RMSD vs frame 0 (ECHO 0 = all atoms)
        if self.rmsd_frame_0():

            #2.1 Check the existance of necessary files
            self.check_file("{file}.xtc".format(file = file_name))
            self.check_file("{file}.tpr".format(file = file_name))

            #2.2 Extract Frame 0 from .xtc trajectory
            self.update_terminal("Extracting {file}_frame_0.pdb from {file}.xtc ..."\
                .format(file = self.file_name()))
            self.run("echo '1' | gmx trjconv -e 1 -f {file}.xtc -s {file}.tpr -o {file}_frame_0.pdb {debug}"\
                .format(file = self.file_name(), debug = debug))
            
            #2.3 Extract trajectory as .pdb from .xtc, every {skip} frames, and applying PBC conditions
            self.update_terminal("Extracting {file}_frames.pdb from {file}.xtc (every {skip} frames) ..."\
                .format(file = self.file_name(), skip = skip))
            self.run("echo '1' | gmx trjconv -f {file}.xtc -s {file}.tpr -o {file}_frames.pdb -pbc mol -skip {skip} {debug}"\
                .format(file = file_name, skip = skip, debug = debug))
            
            #2.4 Generate RMSD data using 'gmx rms'
            self.update_terminal("Generating RMSD data in {file}_vs_frame_0.xvg ..."\
                .format(file = self.file_name()))
            self.run("echo '0 \n 0' | gmx rms -s {file}_frame_0.pdb -f {file}_frames.pdb -o {file}_vs_frame_0.xvg {debug}"\
                .format(file = file_name, debug = debug))
        
        #3 RMSD vs orgiginal protein (ECHO 0 = all atoms)
        if self.rmsd_original():

            #3.1 If the .pdb trajectory was not written, do so in this step
            if not self.rmsd_frame_0():

                #3.1.1 Check the existance of necessary files
                self.check_file("{file}.xtc".format(file = file_name))
                
                #3.1.2 Extract trajectory as .pdb from .xtc, every {skip} frames, and applying PBC conditions
                self.update_terminal("Extracting {file}_frames.pdb from {file}.xtc (every {skip} frames) ..."\
                    .format(file = self.file_name(), skip = skip))
                self.run("echo '0' | gmx trjconv -f {file}.xtc -s {file}.tpr -o {file}_frames.pdb -pbc mol -skip {skip} {debug}"\
                    .format(file = file_name, skip = skip, debug = debug))
            
            #3.2 Check the existance of necessary files
            self.check_file(self.original())
            
            #3.3 Generate RMSD data using 'gmx rms'
            self.update_terminal("Generating RMSD data in {file}_vs_original.xvg ..."\
                .format(file = self.file_name()))
            self.run("echo '0 \n 0' | gmx rms -s {original} -f {file}_frames.pdb -o {file}_vs_original.xvg -skip {skip} {debug}"\
                .format(original = self.original(), file = file_name, skip = skip, debug = debug))
        
        #4. Potential energy
        if self.potential_ene():

            #4.1 Check the existance of necessary files
            self.check_file("{file}.edr".format(file = file_name))

            #4.2 Generate potential energy data using 'gmx energy'
            self.update_terminal("Generating ENERGY data in {file}_potential.xvg ..."\
                .format(file = self.file_name()))
            self.run("echo '11 \n \n' | gmx energy -f {file}.edr -o {file}_potential.xvg -skip {skip} {debug}"\
                .format(file = file_name, skip = skip, debug = debug))
        
        #5. Gyration degrees
        if self.gyration_degr():

            #5.1 Check the existance of necessary files
            self.check_file("{file}.tpr".format(file = file_name))

            #5.2 Generate gyration data using 'gmx gyrate'
            self.update_terminal("Generating GYRATE data in {file}_gyrate.xvg ..."\
                .format(file = self.file_name()))
            self.run("echo 1 | gmx gyrate -f {file}.xtc -s {file}.tpr -o {file}_gyrate.xvg -fitfn exp {debug}"\
                .format(file = file_name, skip = skip, debug = debug))


    def analyse_data(self):
        """
        Read produced files and plot recovered data.
        """

        #1. Define necessary variables
        file_name = self.dir() + '/' + self.file_name()
        f = plt.figure()
        
        #2. Define requested graphs position in the final layout
        graphs = []
        index_start = 101
        if self.rmsd_frame_0():
            graphs.append('_vs_frame_0')
            index_start += 10
        if self.rmsd_original():
            graphs.append('_vs_original')
            index_start += 10
        if self.potential_ene():
            graphs.append('_potential')
            index_start += 10
        if self.gyration_degr():
            graphs.append('_gyrate')
            index_start += 10
        if len(graphs) > 2:
            index_start = 221

        #3. Plot data
        plots = {}
        for i in xrange(0, len(graphs)):
            self.update_terminal("Plotting data from {file}{graph}.xvg ..."\
                .format(file = self.file_name(), graph = graphs[i]))
            
            #3.1 Read file data
            data = self.parse_file(file_name + graphs[i] + '.xvg')
            
            #3.2 Instantiate plot
            index = index_start + i
            plots[graphs[i]] = f.add_subplot(index)
            plot = plots[graphs[i]]
            
            #3.3 Plot recovered data
            plot.plot(data['x'], data['y'],
                label = "{key}".format(key = graphs[i][1:]),
                color = Default.PLOT_COLOR[0])
            
            #3.4 If axis information regarding gyration degrees exist, plot the data
            for i, key in enumerate(['X', 'Y', 'Z']):
                if len(data[key]) > 0:
                    plot.plot(data['x'], data[key],
                    label = "{key} axis".format(key = key),
                    color = Default.PLOT_COLOR[i + 1])
            
            #3.5 Define plot title, legend and axis labels
            plot.legend(frameon = False)
            if 'title' in data.keys():
                plot.set_title(data['title'], fontsize = Default.LABEL_SIZE)
            if 'xlabel' in data.keys():
                plot.set_xlabel(data['xlabel'],  fontsize = Default.LABEL_SIZE)
            if 'ylabel' in data.keys():
                plot.set_ylabel(data['ylabel'],  fontsize = Default.LABEL_SIZE)
            self.update_progress_bar()
        
        #4. Show ploted data
        f.tight_layout()
        plt.show()

        
    def parse_file(self, file_name):
        """
        Read files and extract requested data.
        """

        #1. Define output dictionary
        output = {'x': [], 'y': [], 'X': [], 'Y': [], 'Z': []}
        
        #2. Parse file
        with open(file_name, 'r') as data:
            for line in data:
                elem = line.split()
                if len(elem) <= 1: continue
                
                #2.1 Recover title and axis labels
                if elem[1] == 'xaxis': output['xlabel'] = ' '.join(elem[3:])[1:-1]
                if elem[1] == 'yaxis': output['ylabel'] = ' '.join(elem[3:])[1:-1]
                if elem[1] == 'subtitle' or elem[1] == 'title': output['title'] = ' '.join(elem[2:])[1:-1]
                
                #2.2 Recover X and Y data
                if not line.startswith('@') and not line.startswith('#'):
                    output['x'].append(float(elem[0]))
                    output['y'].append(float(elem[1]))
                    
                    #2.3 If existent, recover gyrate axis data
                    if file_name.endswith("_gyrate.xvg"):
                        output['X'].append(float(elem[2]))
                        output['Y'].append(float(elem[3]))
                        output['Z'].append(float(elem[4]))
        return output


#*****************************************
#*                 MAIN                  *
#*****************************************

#1. Instantiate window and apply themed style
root = tk.ThemedTk()
root.set_theme(Default.THEME)

#2. Instantiate GUI and run main loop
app = Utils(root)
root.mainloop()