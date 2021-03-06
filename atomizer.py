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
from itertools import cycle

class Default:
    """
    Holds default values for various variables.
    """

    #Application window settings
    TITLE          = "Atomizer"
    FOOTER         = "José Pereira @ 2018"
    VERSION        = 1.1
    RESIZABLE      = False
    INPUT_ERROR    = "Input error"
    
    #Input defaults
    DATABASE_NAME  = '1ctf_7_5_2018'
    DATABASE_HOST  = '192.168.179.114'
    DATABASE_PORT  = 27017
    DATABASE_T_ID  = 0
    FROM           = 0
    TO             = 0
    TRAJECTORY     = "traj.xtc"
    OUTPUT         = "traj"
    METHOD         = "gromos"
    STANDARD_DIR   = os.getcwd()
    STANDARD       = "model.pdb"
    PLOT_NAME      = "clusters.png"
    CUT_OFF        = 0.5
    PRINT_EVERY    = 100
    FRAME          = 1
    DEBUG_MODE     = True
    ONLY_PRINT_XTC = False
    PLOT_ENERGIES  = True
    ANNOTATE       = True

    #Model specific defaults
    PDB_FORMAT     = "ATOM {n:>6d}  {symbol:<3s} {residue:3s} {chain:<2s} {residue_n:>2d} {x:>11.3f} {y:>7.3f} {z:>7.3f}  1.00  0.00\n"
    ATOMS_PER_RESI = 6
    ATOM_NAMES     = cycle(['H', 'N', 'CA', 'C', 'O', 'X'])

    #Application theme settings
    PAD            = 10
    FONT_SIZE      = 8
    LABEL_WIDTH    = 170
    BLOB_SCALE     = 10
    THEME          = "arc"
    BACKGROUND     = "white"
    HOVER_COLOR    = "steel blue"
    PLOT_COLOR     = "#4682b4"


class Utils(BasicGUIFunctions):

    def __init__(self, master):
        """
        Creates a new application window with corrects settings and GUI layout.
        Directory and file selection defaults to current directory.
        """
        
        #Set application's window settings and title
        self.config_master(master, Default.TITLE, Default.VERSION, Default.RESIZABLE, Default.BACKGROUND)

        #Set GUI elements
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
        self.add_basic_utilities(master, 10, default_text = "Atomize to {output}_all_atoms.pdb"\
            .format(output = Default.OUTPUT), pad = Default.PAD, background = Default.BACKGROUND,
            hover_color = Default.HOVER_COLOR, footer = Default.FOOTER)


    def process(self):
        """
        Process control. Generate and analyze data.
        """

        #1. Set progress bar settings
        self.progress_bar['maximum'] = self.to_do()
        self.progress_bar['value'] = 0
        self.debug = '' if self.debug_mode() else '2>&1'

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
        f = None
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
        to_do += 4       #extract_frame
        to_do += 4       #all_atoms
        if self.plot_energ():
            to_do += 3   #plot_energies
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
    
    
    def extract_frame(self, log = Default.OUTPUT + '.log', frame = 1,
        output_name = Default.OUTPUT + '_repr_str.pdb'):
        """
        Read clusters.log file and extract the frame number of the representative
        structure of the selected cluster;
        Extract the frame XYZ from MONGO;
        Recover aminoacid sequence from standard pdb structure;
        Use the extracted XYZ and aminoacid sequence to write the PDB of the representative structure
        of the selected cluster;

        Returns the representative structure XYZ.
        """
        
        #1. Open cluster.log and extract the actual structure number for the representative frame
        # of the selected cluster
        self.update_terminal("Finding the selected representative structure from cluster {cluster} in {log} ..."\
            .format(cluster = frame, log = log))
        with open(log, 'r') as file_in:
            for line in file_in:
                elem = line.split()
                if len(elem) <= 0:
                    continue
                if not elem[0] == '|' and elem[1] == '|' and not line.startswith('cl.') and int(elem[0]) == frame:
                    structure_number = int(self.crop.start()) + int(float(elem[5]))
                    break
        self.update_progress_bar()

        #2. Query MONGO database for the found representative structure number and recover XYZ
        self.update_terminal("Extracting the selected structure ({model}) from MONGO database ..."\
            .format(model = structure_number))
        with pymongo.MongoClient(self.db_host(), int(self.db_port())) as server:
            temperature = 'structures{t_id}'.format(t_id = self.db_t_id())
            data = server[self.db_name()][temperature]
            query = {'xyz': 1, '_id': 0}
            xyz = data.find({'step': structure_number}, query)
        self.update_progress_bar()

        #3. Recover aminoacid sequence from the standard PDB structure.
        self.update_terminal("Recovering aminoacid sequence from {standard} ..."\
            .format(standard = self.standard()))
        aa_list = []
        with open(self.standard(), 'r') as file_in:
            for line in file_in:
                elem = line.split()
                if elem[0] == 'ATOM' and elem[2] == 'H':
                    aa_list.append(elem[3])
        self.update_progress_bar()

        #4. Use the recovered XYZ and aminoacid sequence to write the representative structure PDB
        self.update_terminal("Writting {output_file} with the selected structure ..."\
            .format(output_file = output_name))
        naa = 0
        with open(output_name, 'w') as file_out:
            atoms = xyz[0]['xyz']
            for i, atom in enumerate(atoms, start = 1):
                symbol = next(Default.ATOM_NAMES)
                file_out.write(
                    Default.PDB_FORMAT.format(
                        n = i, symbol = symbol, residue = aa_list[naa], chain = 'A',
                        residue_n = naa + 1, x = atom[0], y = atom[1], z = atom[2]
                    )
                )
                if symbol == 'X':
                    naa += 1
        self.update_progress_bar()
        
        return xyz[0]['xyz']


    def cluster(self, trajectory = Default.TRAJECTORY, output = Default.OUTPUT,
        standard = Default.STANDARD_DIR + '/' + Default.STANDARD, method = Default.METHOD):
        """
        Use GMX CLUSTER to create clusters and necessary files.
        """

        #1. Define necessary variables
        cut_off = float(self.cut_off())

        #2. Clusterization (ECHO 0 = all atoms)
        self.update_terminal("Clustering trajectory frames to {output_name}_clusters.pdb with cut off {cut_off} ..."\
                .format(output_name = output, cut_off = cut_off))
        self.run("echo '0 \n 0' | gmx cluster -f {trajectory} -s {standard} -o {output}.xpm -g {output}.log -cl {output}_clusters.pdb -cutoff {cut_off} -method {method} -dist {debug}"\
            .format(trajectory = trajectory, standard = standard, output = output,
                cut_off = cut_off, method = method, debug = self.debug))
        #self.update_progress_bar()

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
        rmsd, total_rmsd   = [], []
        energy = []
        labels = []

        #2. Run GMX RMS to calculate correct RMSD vs original protein (ECHO 0 = all atoms)
        self.run("echo '0 \n 0' | gmx rms -f {trajectory}_clusters.pdb -s {standard} -o {output} {debug}"\
            .format(trajectory = trajectory, standard = standard, output = output, debug = self.debug))
        #3. Read produced .xvg file and extract RMSD values
        with open(output, 'r') as file_in:
            for line in file_in:
                if not line.startswith('@') and not line.startswith('#'):
                    elem = line.split()
                    total_rmsd.append(float(elem[1])*10)
                
        #4. Read log file and extract clusters' above defined threshold energies and cluster sizes
        self.update_terminal("Reading clusters from {log_file} ...".format(log_file = log_file))
        max_value, min_value = None, None
        with open(log_file, 'r') as file_in:
            for line in file_in:
                elem = line.split()
                if len(elem) <= 0: continue
                if len(elem) > 9 and not elem[0] == '|' and not elem[0] == 'cl.':
                    if max_value == None:
                        min_value = line.find(elem[5])
                        max_value = line.find(elem[5]) + len(elem[5])
                    rmsd.append(total_rmsd[int(line[0:4]) - 1])
                    try:
                        energy.append(self.energy[int(line[min_value:max_value])])
                    except ValueError:
                        self.error("Input Error", 'Cut-Off is too big, the resulting RMSD values between structures of the same cluster exceeds 10 Angstrom.')
                    labels.append((elem[0], elem[2], int(line[min_value:max_value])))
        self.max_cluster = len(rmsd) + 1
        self.update_progress_bar()

        #5. Plot RMSD vs energy
        self.update_terminal("Plotting RMSD vs energy ...")
        plt.scatter(rmsd, energy, s = [int(l[1]) * Default.BLOB_SCALE for l in labels], color = Default.PLOT_COLOR)
        size_l = [int(l[1]) * Default.BLOB_SCALE for l in labels]
        with open("cluster_info.dat", "w") as file_out:
            for i in xrange(0, len(rmsd)):
                #print type(size_l[i]), type(energy[i]), type(rmsd[i])
                file_out.write("{:>7.3f} {:>7.3f} {:>7d}\n".format(rmsd[i], energy[i], size_l[i]))

        plt.xlabel("RMSD (Angstrom)"), plt.ylabel("Energy (KJ/mol)")
        sizes = [int(label[1]) for label in labels]
        if len(sizes) == 0:
            self.error('Input Error', 'Cut-Off is too small or not enough structures: no clusters were built.')
        avg_size, std_size = sum(sizes) / len(sizes), np.std(np.asarray(sizes))
        plt.title("Average cluster size: {avg_size} +/- {std_size:<.3f}"\
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
                frames.limit(int(self.crop.end()) - int(self.crop.start()))
            self.update_progress_bar()

            #5. Write XTC file, applying PRINT EVERY settings, and WITHOUT X's
            self.update_terminal("Writting XTC trajectory to file {output_name} ..."\
                .format(output_name = output_name))
            size = frames.count(with_limit_and_skip = True)
            for n, frame in enumerate(frames):
                if self.debug_mode():
                    print "{:>5d}/{:>5d}".format(n, size)
                if n % int(self.every()) == 0:
                    #5.1 Extract X's
                    frame_without_Xs = []
                    for index, atom in enumerate(frame['xyz'], start = 1):
                        if index % Default.ATOMS_PER_RESI == 0 and index > 1:
                            continue
                        frame_without_Xs.append(atom)
                    #5.2 Write XTC file
                    xtc.write(0.1*np.array(frame_without_Xs), step=n, time=n)
                    #5.3 Save energy values if energy plots is requested
                    if self.plot_energ():
                        self.energy[int(frame['step']) - int(self.crop.start())] = float(frame['energy'])
            
            self.update_progress_bar()


    def verify_input(self):
        """
        Verifies input:
        1) Checks if dabate exists in the requested server
        2) Checks if the temperature ID exists in the selected database
        3) Checks if XYZ data exists for the selected temperature ID
        4) Checks if cut-off is possible
        5) Checks if STANDARD protein pdb file for RMSD comparison exists
        """

        with pymongo.MongoClient(self.db_host(), int(self.db_port())) as server:
            
            self.update_terminal("Connecting to database {name} @ {host}:{port} ..."\
                .format(name = self.db_name(), host = self.db_host(), port = self.db_port()))
            
            #1. Check if database exists
            database_list = [str(db) for db in server.database_names()]
            if not self.db_name() in database_list:
                self.error("Missing database", "{name} not found @ {host}:{port}"\
                    .format(name = self.db_name(), host = self.db_host(), port = self.db_port()))
            self.update_progress_bar()

            self.update_terminal("Fetching temperature ID {t_id} from database {name} @ {host}:{port} ..."\
                .format(t_id = self.db_t_id(), name = self.db_name(), host = self.db_host(), port = self.db_port()))
            
            #2. Check if temperature ID exists
            temperature = 'structures{t_id}'.format(t_id = self.db_t_id())
            temperature_list = [str(temp) for temp in server[self.db_name()].collection_names()]
            if not temperature in temperature_list:
                self.error("Missing temperature ID", "Temperature ID {t_id} not found {name} database"\
                    .format(t_id = self.db_t_id(), name = self.db_name()))
            
            #3. Check if existent temperature has xyz info
            try:
                checker = server[self.db_name()][temperature].find_one({}, {'xyz': 1, '_id': 0})['xyz']
            except KeyError:
                self.error("No XYZ data", "Temperature ID {t_id} does not contain any XYZ data"\
                    .format(t_id = self.db_t_id()))
            self.update_progress_bar()

        #4. Check if cut-off is a possible value
        if self.cut_off() <= 0:
            self.error(Default.INPUT_ERROR, "Cut off value must be greater than 0")
        
        #5. Check if standard file for RMSD comparison exists
        self.check_file(self.standard())


#*****************************************
#*                 MAIN                  *
#*****************************************

#1. Instantiate window and apply themed style
root = tk.ThemedTk()
root.set_theme("arc")

#2. Instantiate GUI and run main loop
app = Utils(root)
root.mainloop()
