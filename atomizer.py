#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Tkinter import *
from ttk import *
from ttkthemes import ThemedStyle
from ttkthemes import themed_tk as tk
from gui import *
import os
import tkMessageBox as alert
from tkSimpleDialog import askinteger as prompt
from mdtraj.formats import XTCTrajectoryFile
import pymongo
import numpy as np
from matplotlib import pyplot as plt
from protosyn.builder import grow
import protosyn.mcmc.caterpillar.caterpillar as cat

class Default:
    """
    Holds default values for various variables.
    """

    #Application window settings
    TITLE          = "Atomizer"
    FOOTER         = "José Pereira @ 2018"
    VERSION        = 1.0
    RESIZABLE      = False
    INPUT_ERROR    = "Input error"
    
    #Input defaults
    DATABASE_NAME  = '1ctf_19_4_2018'
    DATABASE_HOST  = '192.168.179.114'
    DATABASE_PORT  = 27017
    DATABASE_T_ID  = 6
    FROM           = 0
    TO             = 0
    TRAJECTORY     = "traj.xtc"
    OUTPUT         = "traj"
    METHOD         = "linkage"
    STANDARD_DIR   = os.getcwd()
    STANDARD       = "model.pdb"
    PLOT_NAME      = "clusters.png"
    CUT_OFF        = 0.5
    PRINT_EVERY    = 10
    FRAME          = 1
    DEBUG_MODE     = False
    ONLY_PRINT_XTC = False
    PLOT_ENERGIES  = True
    ANNOTATE       = True

    #Application theme settings
    PAD            = 10
    FONT_SIZE      = 8
    LABEL_WIDTH    = 170
    BLOB_SCALE     = 10
    THEME          = "arc"
    BACKGROUND     = "white"
    HOVER_COLOR    = "steel blue"
    PLOT_COLOR     = "#4682b4"


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

        self.db_name    = Input(master, "Database name:", Default.DATABASE_NAME,  row = 0, label_width = Default.LABEL_WIDTH)
        self.db_host    = Input(master, "Database host:", Default.DATABASE_HOST,  row = 1, label_width = Default.LABEL_WIDTH)
        self.db_port    = Input(master, "Database port:", Default.DATABASE_PORT,  row = 2, label_width = Default.LABEL_WIDTH)
        self.db_t_id    = Input(master, "Temperature ID:",Default.DATABASE_T_ID,  row = 3, label_width = Default.LABEL_WIDTH)
        self.crop       = Slice(master, "From:", Default.FROM, "To:", Default.TO, row = 4, label_width = Default.LABEL_WIDTH)
        self.every      = Input(master, "Print every:",   Default.PRINT_EVERY,    row = 5, label_width = Default.LABEL_WIDTH)
        Title(master, "Selecting 'To:' as 0 outputs the trajectory until the end.", row = 6, font_size = 5)
        self.cut_off    = Input(master, "Clustering cut-off:",       Default.CUT_OFF,        row = 7,
            label_width = Default.LABEL_WIDTH)
        self.standard   = SelectFile(master, "Original protein:", default_value = Default.STANDARD_DIR + '/' + Default.STANDARD,
            row = 8, label_width = Default.LABEL_WIDTH)
        frame = Frame(master, background = Default.BACKGROUND)
        frame.grid(row = 9, sticky = 'w')
        self.debug_mode = CheckBox(frame, "Debug Mode", default_status = Default.DEBUG_MODE)
        self.print_xtc  = CheckBox(frame, "Only print XTC", column = 1, default_status = Default.ONLY_PRINT_XTC)
        self.plot_energ = CheckBox(frame, "Plot energies vs RMSD", column = 2, default_status = Default.PLOT_ENERGIES)
        self.annotate   = CheckBox(frame, "Anotate energy plots", column = 3, default_status = Default.ANNOTATE)
        Button(master, text="Atomize to {output}_all_atoms.pdb".format(output = Default.OUTPUT),
            width = 50, command= self.process, activebackground = Default.HOVER_COLOR,
            relief = 'flat', background = Default.BACKGROUND)\
            .grid(row = 10, padx = Default.PAD, columnspan = 2)
        self.progress_bar = Progressbar(master, orient = 'horizontal', mode = 'determinate', length = 400)
        self.progress_bar.grid(columnspan = 2, row = 11, pady = Default.PAD)
        self.terminal     = Terminal(master, row = 12, columnspan = 2)
        self.footer       = Footer(master, Default.FOOTER, row = 13, columnspan = 2)


    def process(self):
        """
        Process control. Generate and analyze data.
        """

        #1. Set progress bar settings
        self.progress_bar['maximum'] = self.to_do()
        self.progress_bar['value'] = 0

        #2. Verify input
        self.verify_input()

        #3. Retrieve data from database as XTC file.
        self.dump_xtc()

        #4. If "Only print xtc" option is selected, exit runtime
        if self.print_xtc():
            self.clean()
            return

        #5. Cluster using GMX CLUSTER
        self.cluster(standard = self.standard())

        #6. If "Plot energies" option is selected, select a single cluster from plot examination
        if self.plot_energ():
            f = self.max_cluster + 1
            while f > self.max_cluster:
                f = prompt("Pick a cluster", "Select the cluster number to extract structure (Default: 1)")
                if f > self.max_cluster:
                    alert.showerror(Default.INPUT_ERROR, "{f} is above maximum number of clusters {max}"\
                        .format(f = f, max = self.max_cluster))
        if f == None: f = Default.FRAME

        #7. Extract selected frame from trajectory
        frame_xyz = self.extract_frame(frame = f)

        #8. Convert selected frame to all_atoms
        self.all_atoms(frame_xyz)

        #9. Delete unnecessary files and exit runtime
        self.clean()


    def clean(self):
        """
        Delete unnecessary files.
        """

        self.update_terminal("Cleaning unnecessary files ...")
        self.run("rm -rf *#")
        self.update_terminal("All tasks successefully performed.")


    def to_do(self):
        """
        Count the total number of tasks requested to proper set the progress bar.
        """
        to_do = 1        #clean
        to_do += 3       #verify input
        to_do += 2       #dump_xtc
        if self.print_xtc():
            return to_do
        to_do += 1       #cluster
        to_do += 1       #extract_frame
        to_do += 4       #all_atoms
        if self.plot_energ():
            to_do += 2   #plot_energies
        return to_do


    def all_atoms(self, frame_xyz, output_name = Default.OUTPUT + '_all_atoms.pdb'):
        """
        Grows a Catterpillar model from the original protein sequence.
        Applies the desired frame XYZ position to the grown model.
        Applies the desired frame dihedral angles to the grown model.
        Export data as an all atom PDB file.
        """
        
        #1. Load orginal protein molecule from PDB file
        self.update_terminal("Growing model from original protein sequence ...")
        ref = cat.Molecule.LoadFromFile(self.standard())
        #2. Extract sequence from original protein molecule and apply it to a new caterpillar model
        mol = grow('all_atoms', ref.get_sequence())
        #3. Create new model from grown molecule
        seq, xyz = cat.get_model_from_molecule(mol)
        model = cat.get_molecule_from_model(seq, xyz)
        self.update_progress_bar()
        
        #4. Apply XYZ positions
        self.update_terminal("Applying the selected frame atom positions ...")
        frame_xyz = np.asarray(frame_xyz)
        model.set_coordinates(frame_xyz)
        self.update_progress_bar()

        #5. Apply Dihedral angles
        self.update_terminal("Applying the selected frame dihedral angles ...")
        for rtarget, rmodel in zip(mol.residues, model.residues):
            for dihd in rmodel.dihedrals:
                rtarget.dihedrals[dihd].degrees = rmodel.dihedrals[dihd].degrees
        self.update_progress_bar()

        #6. Export molecule as PDB
        self.update_terminal("Writting all_atoms model to {output} ...".format(output = output_name))
        with open(output_name, 'w') as file_out:
            file_out.write(mol.as_pdb(include_bonds=True))
        self.update_progress_bar()
    
    
    def extract_frame(self, cluster = Default.OUTPUT + '_clusters.pdb', frame = 1):
        """
        Read clusters.pdb file and extract selected frame
        """
        
        #1. Define necessary variables
        extract = False
        xyz = []
        if self.debug_mode():
            print "FRAME", frame

        #2. Read PDB and extract the desired frame
        self.update_terminal("Extracting the selected model ({model}) from clusters ..."\
            .format(model = frame))
        with open(cluster, 'r') as file_in:
            for line in file_in:
                elem = line.split()
                if len(elem) <= 0: continue
                if elem[0] == 'MODEL' and elem[1] == str(frame):
                    extract = True
                    continue
                if elem[0] == 'TER' and extract:
                    break
                if extract:
                    xyz.append([float(elem[6]), float(elem[7]), float(elem[8])])
        self.update_progress_bar()
        return xyz


    def cluster(self, trajectory = Default.TRAJECTORY, output = Default.OUTPUT,
        standard = Default.STANDARD_DIR + '/' + Default.STANDARD, method = Default.METHOD):
        """
        Use GMX CLUSTER to create clusters and necessary files.
        """

        #1. Define necessary variables
        debug = '' if self.debug_mode() else '2>&1'
        cut_off = float(self.cut_off())

        #2. Clusterization (ECHO 0 = all atoms)
        self.update_terminal("Clustering trajectory frames to {output_name}_clusters.pdb with cut off {cut_off} ..."\
                .format(output_name = output, cut_off = cut_off))
        self.run("echo '0 \n 0' | gmx cluster -f {trajectory} -s {standard} -o {output}.xpm -g {output}.log -cl {output}_clusters.pdb -cutoff {cut_off} -method {method} -dist {debug}"\
            .format(trajectory = trajectory, standard = standard, output = output,
                cut_off = cut_off, method = method, debug = debug))
        self.update_progress_bar()

        #3. If energy plot is requested, run plot_energies()
        if self.plot_energ():
            self.plot_energies("{output}.log".format(output = output), standard = self.standard())


    def plot_energies(self, log_file, trajectory = Default.OUTPUT, output = Default.OUTPUT + '.xvg',
        standard = Default.STANDARD_DIR + '/' + Default.STANDARD,):
        """
        Run GMX RMS to calculate correct RMSD values for the reference structure of each cluster
        vs the original standard protein.
        Read LOG file and extrat clusters (ABOVE DEFINED THRESHOLD) size, energies and labels.
        Plot RMSD vs ENERGY as a scatter plot in MATPLOTLIB. Annotate plot if requested.
        """

        #1. Define necessary variables
        debug = '' if self.debug_mode() else '2>&1'
        rmsd, total_rmsd   = [], []
        energy = []
        if self.annotate():
            labels = []

        #2. Run GMX RMS to calculate correct RMSD vs original protein (ECHO 0 = all atoms)
        self.run("echo '0 \n 0' | gmx rms -f {trajectory}_clusters.pdb -s {standard} -o {output} {debug}"\
            .format(trajectory = trajectory, standard = standard, output = output, debug = debug))
        #3. Read produced .xvg file and extract RMSD values
        with open(output, 'r') as file_in:
            for line in file_in:
                if not line.startswith('@') and not line.startswith('#'):
                    elem = line.split()
                    total_rmsd.append(float(elem[1])*10)
                
        #4. Read log file and extract clusters' above defined threshold energies and cluster sizes
        self.update_terminal("Reading clusters from {log_file} ...".format(log_file = log_file))
        with open(log_file, 'r') as file_in:
            for line in file_in:
                elem = line.split()
                if len(elem) <= 0: continue
                if len(elem) > 9 and not elem[0] == '|' and not elem[0] == 'cl.':
                    if self.debug_mode():
                        print line
                    rmsd.append(total_rmsd[int(line[0:4]) - 1])
                    energy.append(self.energy[int(line[19:25])])
                    if self.annotate():
                        labels.append((elem[0], elem[2], int(line[19:25])))
        self.max_cluster = len(rmsd) + 1
        self.update_progress_bar()

        #5. Plot RMSD vs energy
        self.update_terminal("Plotting RMSD vs energy ...")
        plt.scatter(rmsd, energy, s = [int(l[1]) * Default.BLOB_SCALE for l in labels], color = Default.PLOT_COLOR)
        plt.xlabel("RMSD (Angstrom)"), plt.ylabel("Energy (KJ/mol)")
        sizes = [int(label[1]) for label in labels]
        avg_size, std_size = sum(sizes) / len(sizes), np.std(np.asarray(sizes))
        plt.title("Average cluster size: {avg_size} +/- {std_size}"\
            .format(avg_size = avg_size, std_size = std_size))
        #6. If requested, add labels to data points
        if self.annotate():
            for i, label in enumerate(labels):
                a = plt.annotate("Cluster {i}\nSize {j}\nStructure {n:>3d}"\
                    .format(i = label[0], j = label[1], n = int(label[2]) + int(self.crop.start())),
                    (rmsd[i],energy[i]), size = Default.FONT_SIZE, ha = "left")
        #7. Save figure
        plt.savefig(Default.PLOT_NAME)
        plt.show()
        self.update_progress_bar()


    def dump_xtc(self, output_name = Default.TRAJECTORY):
        """
        Extract XYZ data from MONGO database and save as XTC trajectory file.
        """
        #1. Connect to MONGO @ host:port
        with pymongo.MongoClient(self.db_host(), int(self.db_port())) as server,\
            XTCTrajectoryFile(output_name, 'w') as xtc:

            #2. Define necessary variables
            temperature = 'structures{t_id}'.format(t_id = self.db_t_id())
            data = server[self.db_name()][temperature]
            
            #3. Recover necessary data from MONGO
            self.update_terminal("Recovering data from database '{name}' @ {host}:{port} for temperature ID {t_id} ..."\
                .format(name = self.db_name(), host = self.db_host(), port = self.db_port(), t_id = self.db_t_id()))
            query = {'xyz': 1, '_id': 0}
            if self.plot_energ():
                query['energy'] = 1
                query['step'] = 1
                self.energy = {}
            frames = data.find({}, query)

            #4. Apply FROM and TO (crop data)
            if int(self.crop.start()) > 0:
                if int(self.crop.start()) > frames.count():
                    self.error(Default.INPUT_ERROR, "Minimum crop value is {max}".format(max = frames.count()))
                frames.skip(int(self.crop.start()))
            if int(self.crop.end()) > 0:
                if int(self.crop.end()) > frames.count():
                    self.error(Default.INPUT_ERROR, "Maximum crop value is {max}".format(max = frames.count()))
                frames.limit(int(self.crop.end()))
            self.update_progress_bar()

            #5. Write XTC file, applying PRINT EVERY settings
            self.update_terminal("Writting XTC trajectory to file {output_name} ..."\
                .format(output_name = output_name))
            size = frames.count(with_limit_and_skip = True)
            for n, frame in enumerate(frames):
                if self.debug_mode():
                    print "{:>5d}/{:>5d}".format(n, size)
                if n % int(self.every()) == 0:
                    xtc.write(0.1*np.array(frame['xyz']), step=n, time=n)
                    #5.1 Save energy values if energy plots is requested
                    if self.plot_energ():
                        self.energy[int(frame['step']) - int(self.crop.start())] = float(frame['energy'])
            self.update_progress_bar()


    def verify_input(self):
        """
        Verifies input:
        1) Checks if dabate exists in the requested server
        2) Checks if the temperature ID exists in the selected database
        3) Checks if cut-off is possible
        4) Checks if STANDARD protein pdb file for RMSD comparison exists
        """

        with pymongo.MongoClient(self.db_host(), int(self.db_port())) as server:
            
            self.update_terminal("Connecting to database {name} @ {host}:{port} ..."\
                .format(name = self.db_name(), host = self.db_host(), port = self.db_port()))
            database_list = [str(db) for db in server.database_names()]
            if not self.db_name() in database_list:
                self.error("Missing database", "{name} not found @ {host}:{port}"\
                    .format(name = self.db_name(), host = self.db_host(), port = self.db_port()))
            self.update_progress_bar()

            self.update_terminal("Fetching temperature ID {t_id} from database {name} @ {host}:{port} ..."\
                .format(t_id = self.db_t_id(), name = self.db_name(), host = self.db_host(), port = self.db_port()))
            temperature = 'structures{t_id}'.format(t_id = self.db_t_id())
            temperature_list = [str(temp) for temp in server[self.db_name()].collection_names()]
            if not temperature in temperature_list:
                self.error("Missing temperature ID", "Temperature ID {t_id} not found {name} database"\
                    .format(t_id = self.db_t_id(), name = self.db_name()))
            self.update_progress_bar()

        if self.cut_off() <= 0:
            self.error(Default.INPUT_ERROR, "Cut off value must be greater than 0")
        self.check_file(self.standard())


    def check_file(self, file_name):
        """
        Checks the existance of a file in the requested directory.
        Raises a error and exits if it does not exist.
        """

        self.update_terminal("Checking for the existance of file {file} ..."\
            .format(file = file_name))
        if not os.path.isfile(file_name):
            self.error("Missing File", "{file} not found".format(file = file_name))
        self.update_progress_bar()


    def error(self, error_title, message):
        """
        Displays a pop-up alert to user and exits program.
        """

        alert.showerror(error_title, message)
        exit(1)


    def run(self, command):
        """
        Instantiate a new terminal command (GROMACS)
        """

        os.popen(command)
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

#*****************************************
#*                 MAIN                  *
#*****************************************

#1. Instantiate window and apply themed style
root = tk.ThemedTk()
root.set_theme("arc")

#2. Instantiate GUI and run main loop
app = Utils(root)
root.mainloop()
