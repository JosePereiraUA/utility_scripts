#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Tkinter import *
from ttk import *
from ttkthemes import ThemedStyle
from ttkthemes import themed_tk as tk
from gui import *
import pymongo
from itertools import cycle
import json

class Default(Global):
    """
    Holds default values for various variables.
    """

    #Application window settings
    TITLE          = "Extractor"
    FOOTER         = "José Pereira @ 2018"
    VERSION        = "_DEV"
    RESIZABLE      = False

    #Input defaults
    PROTEIN_SEQ    = "EFDVILKAAGANKVAVIKAVRGATGLGLKEAKDLVESAPAALKEGVSKDDAEALKKALEEAGAEVEVK"
    DATABASE_NAME  = '1ctf_7_5_2018'
    DATABASE_HOST  = '192.168.179.114'
    DATABASE_PORT  = 27017
    DATABASE_T_ID  = 0
    FROM           = 30000
    TO             = 30001
    EVERY          = 1
    DEBUG_MODE     = True
    REMOVE_Xs      = True
    OUTPUT_TYPE    = 0    #0 - As PDB; 1 - As XYZ;

    #Runtime defaults
    OUTPUT         = "teste"
    ATOMS_PER_RESI = 6
    ATOM_NAMES     = ['H', 'N', 'CA', 'C', 'O', 'X']
    PDB_FORMAT     = "ATOM {n:>6d}  {symbol:<3s} {residue:3s} {chain:<2s} {residue_n:>2d} {x:>11.3f} {y:>7.3f} {z:>7.3f}  1.00  0.00\n"


    #Application theme settings
    THEME          = "arc"
    LABEL_WIDTH    = 170


class Utils(BasicGUIFunctions):

    def __init__(self, master):
        """
        Creates a new application window with corrects settings and GUI layout.
        Directory and file selection defaults to current directory.
        """

        #Set application's window settings and title
        self.config_master(master, Default.TITLE, Default.VERSION, Default.RESIZABLE, Default.BACKGROUND)

        #Set GUI elements
        self.add_database_picker(master, row = 0, default_database_name = Default.DATABASE_NAME,
            default_database_host = Default.DATABASE_HOST, default_database_port = Default.DATABASE_PORT,
            default_database_t_id = Default.DATABASE_T_ID, default_from = Default.FROM,
            default_to = Default.TO, default_every = Default.EVERY,
            label_width = Default.LABEL_WIDTH, background = Default.BACKGROUND)
        
        self.sequence    = Input(master, "Protein sequence:", row = 1, default_value = Default.PROTEIN_SEQ,
            width = 40, label_width = Default.LABEL_WIDTH, columnspan = 2)
        self.output_type = RadioButtonGroup(master, row = 2, button_texts = ["as PDB", "as XYZ"],
            default_state = Default.OUTPUT_TYPE)

        frame = Frame(master, background = Default.BACKGROUND)
        frame.grid(row = 3, sticky = 'w')
        self.debug_mode = CheckBox(frame, "Debug Mode", default_status = Default.DEBUG_MODE)
        self.remove_Xs  = CheckBox(frame, "Remove X's", default_status = Default.REMOVE_Xs, column = 1)

        self.add_basic_utilities(master, row = 4, default_text = "Extract to {output}"\
            .format(output = Default.OUTPUT), pad = Default.PAD, background = Default.BACKGROUND,
            hover_color = Default.HOVER_COLOR, footer = Default.FOOTER)

    def process(self):
        
        self.progress_bar['maximum'] = self.to_do()
        self.progress_bar['value'] = 0
        self.debug = '' if self.debug_mode() else '2>&1'

        #1.
        self.verify_input(Default.OUTPUT)
        
        #2.
        self.extract_output_file()

        #3.
        self.clean()
    

    def extract_output_file(self, output_name = Default.OUTPUT):
        
        with pymongo.MongoClient(self.db_host(), int(self.db_port())) as server:
            temperature = 'structures{t_id}'.format(t_id = self.db_t_id())
            data = server[self.db_name()][temperature]

            self.update_terminal("Recovering data from database '{name}' @ {host}:{port} for temperature ID {t_id} ..."\
                .format(name = self.db_name(), host = self.db_host(), port = self.db_port(), t_id = self.db_t_id()))
            query = {'xyz': 1, '_id': 0}
            frames = data.find({}, query)
            if int(self.crop.start()) > 0:
                if int(self.crop.start()) > frames.count():
                    self.error(Default.INPUT_ERROR, "Minimum crop value is {max}".format(max = frames.count()))
                frames.skip(int(self.crop.start()))
            if int(self.crop.end()) > 0:
                if int(self.crop.end()) > frames.count():
                    self.error(Default.INPUT_ERROR, "Maximum crop value is {max}".format(max = frames.count()))
                frames.limit(int(self.crop.end()) - int(self.crop.start()))
            self.update_progress_bar()

            self.update_terminal("Writting trajectory to file {output_name} ..."\
                .format(output_name = output_name))
            size = frames.count(with_limit_and_skip = True)
            
            if self.remove_Xs:
                atom_names = cycle(Default.ATOM_NAMES[:-1])
                atoms_per_resi = Default.ATOMS_PER_RESI - 1
            else:
                cycle(Default.ATOM_NAMES)
                atoms_per_resi = Default.ATOMS_PER_RESI

            if self.output_type() == 0:
                sequence = self.load_sequence()
            
            i = 0
            for n, frame in enumerate(frames):
                if self.debug_mode():
                    print "{:>5d}/{:>5d}".format(n, size)
                if n % int(self.every()) == 0:
                    i += 1
                    #1 Extract X's
                    if self.remove_Xs():
                        frame_without_Xs = []
                        for index, atom in enumerate(frame['xyz'], start = 1):
                            if index % atoms_per_resi == 0 and index > 1:
                                continue
                            frame_without_Xs.append(atom)
                    else:
                        frame_without_Xs = frame

                    if self.output_type() == 0:
                        #1 Write PDB file
                        with open(output_name + '.pdb', 'a') as file_out:
                            file_out.write("TITLE {n}\n".format(n = n + int(self.crop.start())))
                            file_out.write("MODEL {i}\n".format(i = i))
                            residue_n = 0
                            residue = next(sequence)
                            for index, atom in enumerate(frame_without_Xs):
                                if index % atoms_per_resi == 0:
                                    residue = next(sequence)
                                    residue_n += 1
                                file_out.write(Default.PDB_FORMAT.format(
                                    n = index + 1, symbol = next(atom_names), residue = residue,
                                    chain = 'A', residue_n = residue_n,
                                    x = atom[0], y = atom[1], z = atom[2]
                                ))

                    elif self.output_type() == 1:
                        #2 Write XYZ file
                        with open(output_name + '.xyz', 'a') as file_out:
                            file_out.write("{size}\n".format(size = len(frame_without_Xs)))
                            file_out.write("Frame {n}\n".format(n = n + int(self.crop.start())))
                            for atom in frame_without_Xs:
                                file_out.write("{symbol:>2s} {x:>8.3f} {y:>8.3f} {z:>8.3f}\n"\
                                    .format(symbol = next(atom_names), x = atom[0], y = atom[1], z = atom[2]))
            self.update_progress_bar()

    def load_sequence(self): #to_do += 0
        loc = os.path.dirname(os.path.realpath(__file__))
        conversion_table_location = loc + "/static/aa_conversion.json"
        self.check_file(conversion_table_location, update_progress_bar = False)
        with open(conversion_table_location) as json_file:
            conversion_table = json.load(json_file)
        aa_3 = []
        for aa_1 in self.sequence():
            aa_3.append(conversion_table[aa_1])
        return cycle(aa_3)

    
    def to_do(self):
        """
        Count the total number of tasks requested to proper set the progress bar.
        """
        to_do = 1        #clean
        to_do += 3       #verify input
        to_do += 2       #extract output

        return to_do


    def verify_input(self, output_name):
        """
        Verifies input:
        1) Checks if dabate exists in the requested server
        2) Checks if the temperature ID exists in the selected database
        3) Checks if XYZ data exists for the selected temperature ID
        4) Checks if OUTPUT file already exists
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
        
        #4. Check if output file already exists
        if self.output_type() == 0:
            self.check_file_delete(output_name + '.pdb')
        elif self.output_type() == 1:
            self.check_file_delete(output_name + '.xyz')





#*****************************************
#*                 MAIN                  *
#*****************************************

#1. Instantiate window and apply themed style
root = tk.ThemedTk()
root.set_theme(Default.THEME)

#2. Instantiate GUI and run main loop
app = Utils(root)
root.mainloop()