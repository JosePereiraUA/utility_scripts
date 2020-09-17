#!/usr/bin/env python
from Tkinter import Tk, Label, Button, Entry
from mdtraj.formats import XTCTrajectoryFile
from scipy import signal
import matplotlib.pyplot as plt
import Tkinter as tk
import tkMessageBox
import numpy as np
import pymongo
import ttk
import os

"""________________________________________________________________________________
> Default values <
"""
class Default:
    VERSION        = '1.1.0'
    NAME           = 'gaussian'
    HOST           = 'localhost'
    PORT           = '27017'
    INPUT_ERROR    = "Input Error"
    OUTPUT_ERROR   = "Output Error"
    RAW            = "raw_cdf.dat"
    PDB            = "1ctf.pdb"
    RMSD_PREFIX    = "t0"
    CDF_H          = "cdf.h"
    TEMPERATURE_ID = '0'
    ENERGY_TYPES   = ['eBond', 'eH', 'eSol', 'eIJ', 'energy']
    BATCH_STEP     = 15000
    PRINT_EVERY    = 1
    BINS           = 180
    MIN_KERNEL     = 1.0
    MAX_KERNEL     = 20.0
    MIN_TEMP       = 0.2
    MAX_TEMP       = 0.6
    EXPONENTIAL    = 0.15
    COUNT          = 12
    THRESHOLD      = 0
    ALPHA          = 0.01
    BASELINE       = 0.2
    N_POINTS       = 360
    SIZE           = 20
    FROM           = 0
    TO             = 0
    CENTERS        = [(-105, 130), (45, 25), (-85, -30)]
    RADIUS         = 5

"""________________________________________________________________________________
> TOOL FUNCTIONS <
"""
def warning(error_type = 'Error', error_body = 'Something went wrong!'):
    print "[{:>20s} ] : {}".format(error_type, error_body)
    tkMessageBox.showinfo(error_type, error_body)

def find_file_length(file_name):
    count = 0
    with open(file_name, 'r') as f:
        for line in f: count += 1
    return count

def find_database_length(name, host, port, temperature_id):
    structures = 'structures' + str(temperature_id)
    with pymongo.MongoClient(host, port) as mongo:
        return mongo[name][structures].find().count()

def calculate_dihedral(x1, x2, x3, x4):
    v12 = np.asarray(x2) - np.asarray(x1)
    v23 = np.asarray(x3) - np.asarray(x2)
    v34 = np.asarray(x4) - np.asarray(x3)
    v123 = np.cross(v12,v23)
    v234 = np.cross(v23,v34)

    angle = np.arctan2(np.dot(np.cross(v123,v234),v23)/np.sqrt(np.dot(v23,v23)),
                       np.dot(v123,v234))
    return angle

def validate_input_data(name, host, port, temperature):
        with pymongo.MongoClient(host, port) as mongo:
            #Check if database exists
            if not name in mongo.database_names():
                warning(Default.INPUT_ERROR, "Database name '{} not found in Mongo @ {}:{}".format(name, host, port))
                return False
            #Check if temperature is available
            if temperature < 0:
                warning(Default.INPUT_ERROR, "Temperature {} should be an integer equal or greater than 0".format(temperature))
            available_temperatures = mongo[name]['iterations'].distinct("temperature")
            if temperature >= len(available_temperatures):
                warning(Default.INPUT_ERROR, "Temperature {} not found: Maximum temperature ID is {}".format(temperature, len(available_temperatures)-1))
                return False
        return True

def load_raw(fname):
    data = {}
    with open(fname) as fin:
        for line in fin:
            psi, phi, value = map(float, line.split())
            if not data.has_key(phi):
                data[phi] = {}
            data[phi][psi] = value
    for k, v in data.iteritems():
        data[k] = [v[q] for q in sorted(v.keys())]
    return np.array([data[k] for k in sorted(data.keys())])

def query(every = 1, add = {}):
    if not every == 1: a = {'step': {'$mod': [every, 0]}}
    else: a = {}
    if not add == {}: a.update(add)
    return a

def validate_file(output_file):
    if output_file.endswith('xtc'): output_file = output_file.split('.')[0] + '_1.xtc'
    if os.path.isfile(output_file):
        result = tkMessageBox.askquestion(Default.OUTPUT_ERROR, "{} file found. Delete and continue?".format(output_file), icon='warning')
        if result == 'no': return False
    with open(output_file, 'w') as filout: filout.write('')
    return True

def read_rmsd_file(input_file):
    output = []
    with open(input_file, 'r') as f:
        for line in f:
            elem = line.split()
            if len(elem) == 2 and not elem[0].startswith('@') and not elem[0].startswith('#'):
                output.append(float(elem[1]))
    return output

def read_original_pdb(input_file):
    output = []
    with open(input_file, 'r') as original_pdb:
        for i, line in enumerate(original_pdb):
            elem = line.split()
            if i % 5 == 0:
                output.append(elem[3])
    return output

"""________________________________________________________________________________
> GUI CLASSES <
"""
class Input:
    def __init__(self, master, label, default, row_n, column_n = 0, pad_y = 5, width_n = 15, secondary_label = None):
        Label(master, text = label).grid(row = row_n, column = column_n, padx = 10, pady = pad_y, sticky = 'e')
        f = tk.Frame(master)
        f.grid(row = row_n, column = column_n + 1, sticky = 'w')
        self.input = Entry(f, width = width_n)
        self.input.grid(row = row_n, column = 0, pady = pad_y, sticky = 'w')
        self.input.insert(tk.END, default)
        if not secondary_label == None:
            Label(f, text = secondary_label).grid(row = row_n, column = 1, pady = pad_y, sticky = 'w')

    def get(self):
        return self.input.get()

class Range:
    def __init__(self, master, label, start_label, start_default, end_label, end_default, row_n, pad_y = 5, secondary_label = None, column_n = 0, **kwargs):
        if not label == "":
            self.status = tk.BooleanVar(value = False)
            tk.Radiobutton(master, text=label, variable=self.status, value=True).grid(row = row_n, column = 0, padx = 10, pady = pad_y, sticky = 'ne')
            column_n += 1
        f = tk.Frame(master)
        f.grid(row = row_n, column = column_n, columnspan = 4 + column_n, sticky = 'w')
        self.start = Input(f, start_label, start_default, 0, **kwargs)
        self.end   = Input(f, end_label,   end_default,   0, column_n = 2, **kwargs)
        if not secondary_label == None:
            Label(f, text = secondary_label, font=("Helvetica", 8)).grid(row = row_n + 1, column = 0, columnspan = 4 + column_n, sticky = 'w')

    def get_from(self, as_float = False):
        if as_float: return float(self.start.get())
        return int(self.start.get())

    def get_to(self, as_float = False):
        if as_float: return float(self.end.get())
        return int(self.end.get())

    def get_status(self):
        return self.status.get()

class CheckBox:
    def __init__(self, master, label, row_n, default_status = 0, col_span = 2):
        self.status = tk.IntVar(value = default_status)
        tk.Checkbutton(master, text=label, variable=self.status).grid(row = row_n, columnspan = col_span, padx = 15, pady = 5, sticky = 'w')

    def get_status(self):
        return self.status.get()

class Utils:
    """____________________________________________________________________________
    > GUI DEFINITION <
    """
    def __init__(self, master):
        #Set windows size
        self.master = master
        master.title("Utils {}".format(Default.VERSION))
        master.geometry("600x500")
        master.resizable(0, 0)

        #set tabs
        n            = ttk.Notebook(master)
        download_xyz = ttk.Frame(n)
        ramachandran = ttk.Frame(n)
        cdf          = ttk.Frame(n)
        temp         = ttk.Frame(n)
        graph        = ttk.Frame(n)
        n.add(download_xyz, text = 'Trajectory')
        n.add(ramachandran, text = 'Ramachandran')
        n.add(cdf,          text = 'cdf.h')
        n.add(temp,         text = 'Temperatures')
        n.add(graph,        text = 'Energy graphs')
        n.pack(expand = True, fill = tk.BOTH)

        #Tab1 - Download_xyz
        self.type = tk.IntVar(value = 0)
        self.name_x        = Input(download_xyz, "Database name:" , Default.NAME           , 0)
        self.host_x        = Input(download_xyz, "Database host:" , Default.HOST           , 1)
        self.port_x        = Input(download_xyz, "Database port:" , Default.PORT           , 2)
        self.temperature_x = Input(download_xyz, "Temperature ID:", Default.TEMPERATURE_ID , 3)
        self.output_file   = Input(download_xyz, "Output file:"   , Default.NAME           , 4, secondary_label = '.xyz/xtc/pdb')
        as_xyz = tk.Radiobutton(download_xyz, text="As xyz", variable=self.type, value=1)
        as_xyz.grid(row = 5, pady = (10, 0), sticky = 'w')
        as_xyz.select()
        tk.Radiobutton(download_xyz, text="As xtc", variable=self.type, value=2).grid(row = 6, column = 0, pady = (10, 0), sticky = 'w')
        tk.Radiobutton(download_xyz, text="As pdb", variable=self.type, value=3).grid(row = 7, column = 0, pady = (10, 0), sticky = 'w')
        self.pdb_x         = Input(download_xyz, "Requires original pdb file:", Default.PDB , 7, column_n = 1)
        self.slice = Range(download_xyz, "Truncate file?", "From:", Default.FROM, "To:", Default.TO, 8, secondary_label = "If 'To:' is set to 0, the trajectory is read until the end.")
        self.every_x       = Input(download_xyz, "Print every:"   , Default.PRINT_EVERY    , 9)
        b = Button(download_xyz, text = "Write trajectory", command = self.produce_xyz)
        b.grid(columnspan = 3, pady = 5, sticky = 'ew',                                row = 10)
        self.pb_x = ttk.Progressbar(download_xyz, orient = 'horizontal', mode = 'determinate')
        self.pb_x.grid(columnspan = 3, sticky = 'ew',                                  row = 11)

        #Tab2 - Produce Ramachandran
        self.v = tk.IntVar(value = 0)
        from_xyz = tk.Radiobutton(ramachandran, text="From XYZ file:", variable=self.v, value=1)
        from_xyz.grid(pady = (10, 0), sticky = 'w')
        self.input_file = Entry(ramachandran)
        self.input_file.grid(row = 0, column = 1, pady = (10, 0))
        self.input_file.insert(tk.END, Default.NAME)
        Label(ramachandran, text = ".xyz").grid(row = 0, column = 2, pady = (10, 0), sticky = 'w')
        tk.Radiobutton(ramachandran, text="From Mongo Database:", variable=self.v, value=2).grid(row = 2, sticky = 'w')
        self.name_r        = Input(ramachandran, "Database name:",  Default.NAME           , 3)
        self.host_r        = Input(ramachandran, "Database host:",  Default.HOST           , 4)
        self.port_r        = Input(ramachandran, "Database port:",  Default.PORT           , 5)
        self.temperature_r = Input(ramachandran, "Temperature ID:", Default.TEMPERATURE_ID , 6)
        b = Button(ramachandran, text="Plot Ramachandran Map", command=self.produce_ramachandran)
        b.grid(columnspan = 3, pady = 5, sticky = 'ew',                                row = 7)
        self.pb_r = ttk.Progressbar(ramachandran, orient = 'horizontal', mode = 'determinate')
        self.pb_r.grid(columnspan = 3, sticky = 'ew',                                  row = 8)

        #Tab3 - Make cdf.h
        Label(cdf, text = "Kernel parameters:").grid(pady = 5, sticky = 'e')
        self.min_cdf    = Input(cdf, "Minimum probability:", Default.MIN_KERNEL  , 3)
        self.max_cdf    = Input(cdf, "Maximum probability:", Default.MAX_KERNEL  , 4)
        self.alpha_cdf  = Input(cdf, "Alpha:"              , Default.ALPHA       , 5)
        self.size_cdf   = Input(cdf, "Size:"               , Default.SIZE        , 6)
        self.type_c = tk.IntVar(value = 0)
        from_raw = tk.Radiobutton(cdf, text="From raw data:", variable=self.type_c, value=1)
        from_raw.grid(row = 7, column = 0, pady = (0, 0), sticky = 'w')
        from_raw.select()
        self.raw_cdf    = Input(cdf, 'File:', Default.RAW                        , 7 , column_n = 1)
        tk.Radiobutton(cdf, text="Simplified (gaussian):", variable=self.type_c, value=2).grid(row = 8, column = 0, pady = (0, 0), sticky = 'w')
        self.radius     = Input(cdf, 'Radius:'             , Default.RADIUS      , 8 , column_n = 1)
        self.n_points   = Input(cdf, 'Number of points:'   , Default.N_POINTS    , 9 , column_n = 1)
        self.baseline   = Input(cdf, 'Baseline:'           , Default.BASELINE    , 10, column_n = 1)
        self.visualize  = CheckBox(cdf, "Visualize CDF.H"                        , 11)
        b = Button(cdf, text="Create cdf.h", command=self.produce_cdf)
        b.grid(columnspan = 3, pady = 5, sticky = 'ew',                      row = 12)

        #Tab4 - Create temperatures
        self.exp_t       = Input(temp, "Exponential:"           , Default.EXPONENTIAL, 1, width_n = 10)
        self.count_t     = Input(temp, "Number of temperatures:", Default.COUNT      , 2, width_n = 10)
        self.temps       = Range(temp, "", "MIN temperature:"   , Default.MIN_TEMP, "MAX temperature:", Default.MAX_TEMP, 3, width_n = 10)
        self.threshold_t = Input(temp, "Threshold:"             , Default.THRESHOLD  , 4, width_n = 10)
        self.add_temp_t  = Input(temp, "Add temperatures:"      , ""                 , 5, width_n = 30)
        Label(temp, text = "Add multiple temperatures separated by ','", font=("Helvetica", 8)).grid(row = 6, column = 0, columnspan = 2, sticky = 'e')
        b = Button(temp, text="Get temperature range", command=self.produce_temperatures)
        b.grid(columnspan = 2, pady = 5, sticky = 'ew',                          row = 7)
        self.text_box_t = tk.Text(temp, width = 1, height = 5, wrap = 'word', exportselection=1)
        self.text_box_t.grid(columnspan = 2, pady = 5, sticky = 'ew',            row = 8)

        #Tab5 - Make energy graphs
        self.name_g        = Input(graph, "Database name:" , Default.NAME           , 0)
        self.host_g        = Input(graph, "Database host:" , Default.HOST           , 1)
        self.port_g        = Input(graph, "Database port:" , Default.PORT           , 2)
        self.temperature_g = Input(graph, "Temperature ID:", 'all'                  , 3, secondary_label = " Use 'all' to plot all temperatures")
        self.every_g       = Input(graph, "Print every:"   , Default.PRINT_EVERY    , 4)
        self.rmsd          = CheckBox(graph, "Plot RMSD:"                           , 5)
        self.rmsd_prefix   = Input(graph, ""               , Default.RMSD_PREFIX    , 5)
        self.contributions = CheckBox(graph, "Plot energy contributions"            , 6)
        b = Button(graph, text="Plot graphs", command=self.produce_graphs)
        b.grid(columnspan = 2, pady = 5, sticky = 'ew',                         row = 7)
        self.text_box_g    = tk.Text(graph, width = 1, height = 5, wrap = 'word', exportselection=1)
        self.text_box_g.grid(columnspan = 2, pady = 5, sticky = 'ew',           row = 8)

    """____________________________________________________________________________
    > PRIMARY FUNCTIONS <
    """
    def produce_graphs(self):
        if self.contributions.get_status() and self.temperature_g.get() == 'all':
            warning(Default.INPUT_ERROR, "Please define a single temperature index to analize energy contributions")
            return
        if self.rmsd.get_status(): f1,(ax1,ax2,ax3, ax4) = plt.subplots(4)
        else: f1,(ax1) = plt.subplots(1)
        name, host, port= self.name_g.get(), self.host_g.get(), int(self.port_g.get())
        every = int(self.every_g.get())
        with pymongo.MongoClient(host, port) as mongo:
            self.temperature_count = len(mongo[name]['iterations'].distinct("temperature"))
            station_list = mongo[name]['iterations'].distinct("station")
            temperature_indexes = np.arange(0, self.temperature_count)
            if self.temperature_g.get() == 'all':
                f1.canvas.set_window_title('System status for all temperatures')
                for temp in temperature_indexes:
                    structures = 'structures'+str(temp)
                    energies = [float(item['energy']) for item in mongo[name][structures].find(query(every), {'energy': 1, '_id': 0})]
                    # ax3.plot(energies, label = str(temp))
                    h, edges = np.histogram(energies, bins = Default.BINS)
                    bin_centers = 0.5*(edges[:-1] + edges[1:])
                    # ax2.plot(bin_centers, h, label = str(temp))
                    # ax2.legend()
            else:
                f1.canvas.set_window_title('System status for temperature {}'.format(self.temperature_g.get()))
                structures = 'structures' + self.temperature_g.get()
                energies = [float(item['energy']) for item in mongo[name][structures].find(query(every), {'energy': 1, '_id': 0})]
                # ax3.plot(energies, label = self.temperature_g.get())
                h, edges = np.histogram(energies, bins = Default.BINS)
                bin_centers = 0.5*(edges[:-1] + edges[1:])
                # ax2.plot(bin_centers, h, label = self.temperature_g.get())
                # ax2.legend()
            for station in station_list:
                ax1.plot([float(item['temperature']) for item in mongo[name]['iterations'].find({'station': station}, {'temperature': 1, '_id': 0})])
                with open("%s.dat" % (station), "w") as file_out:
                    for item in mongo[name]['iterations'].find({'station': station}, {'temperature': 1, '_id': 0}):
                        for value in [str(item['temperature'])]:
                            file_out.write("\n%s" % (value))
            if self.rmsd.get_status():
                input_prefix = self.rmsd_prefix.get()
                i, current_rmsd = 1, 0
                input_file = input_prefix + '_' + str(i) + '_rmsd.xvg'
                while os.path.isfile(input_file): 
                    rmsd = read_rmsd_file(input_file)
                    self.text_box_g.insert("end", "Structure {structure}: minimum = {min_rmsd:.3f}nm | average = {avg_rmsd:.3f}\n".format(structure = rmsd.index(min(rmsd)) + current_rmsd, min_rmsd = min(rmsd), avg_rmsd = sum(rmsd)/float(len(rmsd))))
                    self.text_box_g.see("end")
                    root.update_idletasks()
                    rmsd_x = np.arange(current_rmsd, current_rmsd + len(rmsd))
                    ax4.plot(rmsd_x, rmsd, color = '#1f77b4')
                    i += 1
                    current_rmsd += len(rmsd)
                    input_file = input_prefix + '_' + str(i) + '.xvg'
                ax4.ylabel = '(nm)'
            f1.tight_layout()
            if self.contributions.get_status():
                f2 = plt.figure(2)
                f2.canvas.set_window_title('Energy contributions for temperature {}'.format(self.temperature_g.get()))
                f2.add_subplot(111)
                structures = 'structures' + str(self.temperature_g.get())
                for energy_type in Default.ENERGY_TYPES:
                    plt.plot([float(item[energy_type]) for item in mongo[name][structures].find(query(every), {energy_type: 1, '_id': 0})], label = energy_type)
                plt.legend()
            plt.show()

    def produce_temperatures(self):
        temp_count, temp_exp = int(self.count_t.get()), float(self.exp_t.get())
        min_temp, max_temp = self.temps.get_from(as_float = True), self.temps.get_to(as_float = True)
        threshold = float(self.threshold_t.get())
        w = temp_exp
        k = temp_count - 1
        x = np.arange(temp_count, dtype=float)
        temps = (np.exp(w*(x-k)) - np.exp(-w*k)) / (1-np.exp(-w*k)) * (max_temp - min_temp) + min_temp
        tempx = np.round(temps, 3)
        out = [tempx[0]]
        for i, n in enumerate(tempx[1:]):
            if (n-out[len(out)-1] >= threshold):
                out.append(n)
        if not self.add_temp_t.get() == '':
            for i in self.add_temp_t.get().split(','):
                out.append(float(i))
        string = str(sorted(out))
        self.text_box_t.insert("end", string + "\n")
        self.text_box_t.see("end")

    def produce_cdf(self):
        alpha, size = float(self.alpha_cdf.get()), float(self.size_cdf.get())
        if self.type_c.get() == 1:
            MIN, MAX = float(self.min_cdf.get()), float(self.max_cdf.get())
            raw = load_raw(self.z.get())
            raw[raw>MAX] = MAX
            raw[raw<MIN] = MIN
            x = np.arange(- size, size + 1)
            X, Y = np.meshgrid(x, x)
            kernel = np.exp(-alpha * (X ** 2 + Y ** 2))
            kernel /= kernel.sum()
            grad = signal.convolve2d(raw, kernel, boundary='symm', mode='same')
            phi = grad.max(1)
            psi = grad.max(0)
            dx = 2.0
            x = np.deg2rad(np.arange(-179.0, 180.0, dx))
            cdf_psi = np.cumsum(psi/np.trapz(psi,x))*np.deg2rad(dx)
            cdf_phi = np.cumsum(phi/np.trapz(phi,x))*np.deg2rad(dx)
            phi /= max(phi)
            psi /= max(psi)
        elif self.type_c.get() == 2:
            baseline, n_points, radius = float(self.baseline.get()), int(self.n_points.get()), int(self.radius.get())
            grad = np.zeros((n_points, n_points))
            x = y = np.linspace(-180, 180, n_points)
            X, Y = np.meshgrid(x, y)
            for center in Default.CENTERS:
                grad = grad + (np.exp(-alpha * ((X-center[0])**2 + (Y-center[1])**2) / radius**2))
            MIN = max(grad.flatten()) * baseline
            grad[grad<MIN] = MIN
            phi = grad.max(1)
            psi = grad.max(0)
            cdf_phi = np.cumsum(phi/np.trapz(phi,x))
            cdf_psi = np.cumsum(psi/np.trapz(psi,y))
        else:
            warning(Default.INPUT_ERROR, "No mode of CDF.H creation selected")
        #Print cdf.h as file
        if os.path.isfile(Default.CDF_H):
            result = tkMessageBox.askquestion(Default.OUTPUT_ERROR, "{} file found. Delete and continue?".format(Default.CDF_H), icon='warning')
            if result == 'no': return
        with open(Default.CDF_H, 'w') as fout:
            print >> fout, '#define CDF_LEN %d'%len(x)
            for name,array in [('psi',cdf_psi),('phi',cdf_phi), ('edges',x)]:
                print >> fout, 'const static double cdf_%s[] = {'%name
                for n,point in enumerate(array):
                    print >> fout, '%12.5f,'%point,
                    if (n%10) == 9:
                        print >> fout,''
                print >> fout, '};'
        if self.visualize.get_status(): self.visualize_cdf(grad, phi, psi, cdf_phi, cdf_psi)

    def visualize_cdf(self, g, phi, psi, phi_c, psi_c):
        fig = plt.figure()
        ticks = [-180, -90, 0, 90, 180]
        ramachandran_projection = plt.subplot2grid((3, 3), (1, 0), colspan = 2, rowspan = 2)
        ramachandran_projection.axes.set_xlim(-180, 180)
        ramachandran_projection.axes.set_xticks(ticks)
        ramachandran_projection.axes.set_ylim(-180, 180)
        ramachandran_projection.axes.set_yticks(ticks)
        phi_projection = plt.subplot2grid((3, 3), (0, 0), colspan = 2)
        phi_projection.axes.set_xlim(-180, 180)
        phi_projection.axes.set_xticks(ticks)
        psi_projection = plt.subplot2grid((3, 3), (1, 2), rowspan = 2)
        psi_projection.axes.set_ylim(-180, 180)
        psi_projection.axes.set_yticks(ticks)
        if self.type_c.get() == 1: x = np.arange(-179, 180, 2)
        else: x = np.linspace(-180, 180, int(self.n_points.get()))
        X, Y = np.meshgrid(x, x)
        ramachandran_projection.pcolormesh(X, Y, g)
        phi_projection.plot(x, psi_c)
        phi_projection.plot(x, psi)
        psi_projection.plot(phi_c, x)
        psi_projection.plot(phi, x)
        fig.tight_layout()
        plt.show()

    def validate_truncation(self, max_trajectory):
        if self.slice.get_from() > max_trajectory:
            warning(Default.INPUT_ERROR, "Minimum value to truncate file ({}) is over the trajectory length ({})".format(self.slice.get_from(), max_trajectory))
            return False
        elif self.slice.get_from() < 0:
            warning(Default.INPUT_ERROR, "Minimum value to truncate file ({}) should be an integer equal or greater than 0".format(self.slice.get_from()))
            return False
        if self.slice.get_to() > max_trajectory:
            warning(Default.INPUT_ERROR, "Maximum value to truncate file ({}) is over the trajectory length ({})".format(self.slice.get_to(), max_trajectory))
            return False
        elif self.slice.get_to() < 0:
            warning(Default.INPUT_ERROR, "Maximum value to truncate file ({}) should be an integer equal or greater than 0".format(self.slice.get_to()))
            return False
        if self.slice.get_to() < self.slice.get_from() and not self.slice.get_to() == 0:
            warning(Default.INPUT_ERROR, "Maximum value to truncate file ({}) should be higher than minimum value ({})".format(self.slice.get_to(), self.slice.get_from()))
            return False
        if not self.slice.get_to() == 0: max_trajectory = self.slice.get_to()
        resulting_trajectory = abs(max_trajectory - self.slice.get_from())
        if int(self.every_x.get()) > resulting_trajectory:
            warning(Default.INPUT_ERROR, "Print every {} exceeds resulting trajectory ({} models)".format(self.every_x.get(), resulting_trajectory))
        return True
        
    def produce_xyz(self):
        self.models = None #RESET
        name, host, port= self.name_x.get(), self.host_x.get(), int(self.port_x.get())
        temperature, output_file = int(self.temperature_x.get()), self.output_file.get()

        #Validate XYZ input
        if not validate_input_data(name, host, port, temperature): return
        if int(self.every_x.get()) <= 0:
            warning(Default.INPUT_ERROR, "Print Every value ({}) should be 1 or higher".format(self.every_x.get()))
            return
            
        #Validate truncation input
        current_models = 0
        max_trajectory = find_database_length(name, host, port, temperature)
        max_models = max_trajectory
        if self.slice.get_status():
            if not self.validate_truncation(max_trajectory): return
            #Slice: FROM
            if self.slice.get_from == 0: current_models = 0
            else: current_models = self.slice.get_from()
            print "from_", current_models
            #Slice: TO
            if self.slice.get_to() == 0 and self.slice.get_status(): max_models = max_trajectory
            else: max_models = self.slice.get_to()
        self.pb_x['value'] = current_models
        self.pb_x['maximum'] = max_models

        if self.type.get() == 1: output_file += '.xyz'
        elif self.type.get() == 2: output_file += '.xtc'
        else:
            residue_list = read_original_pdb(self.pdb_x.get())
            output_file += '.pdb'
        if not validate_file(output_file): return

        run_loop = 0
        while current_models < (max_models/float(self.every_x.get())):
            run_loop += 1
            objective_models = current_models + Default.BATCH_STEP
            if objective_models > max_models: objective_models = max_models
            
            trajectory = self.download_xyz(name, host, port, temperature,
                start = current_models, end = objective_models, every = int(self.every_x.get()))

            #Write output file
            atom_names = ['H', 'N', 'CA', 'C', 'O', 'X']
            #Type 1 - XYZ
            if self.type.get() == 1:
                index = 0
                with open(output_file, 'a') as filout:
                    for model_n, model in enumerate(trajectory):
                        self.pb_x['value'] = (current_models + model_n+1) * int(self.every_x.get())
                        root.update_idletasks()
                        filout.write("  {}\n".format(len(model)))
                        for xyz in model:
                            filout.write("{} {} {} {}\n".format(atom_names[index], xyz[0], xyz[1], xyz[2]))
                            index += 1
                            if index % len(atom_names) == 0: index = 0
            #Type 2 - XTC
            elif self.type.get() == 2:
                f = output_file.split('.')[0] + '_' + str(run_loop) + '.xtc'
                with XTCTrajectoryFile(f, 'w') as xtc:
                    for model_n, model in enumerate(trajectory):
                        self.pb_x['value'] = (current_models + model_n+1) * int(self.every_x.get())
                        xtc.write(0.1*np.array(model), step = model_n, time = model_n)
            #Type 3 - PDB
            elif self.type.get() == 3:
                with open(output_file, 'a') as filout:
                    for model_n, model in enumerate(trajectory):
                        self.pb_x['value'] = (current_models + model_n+1) * int(self.every_x.get())
                        root.update_idletasks()
                        filout.write("MODEL {}\n".format(current_models + model_n+1))
                        index = 0
                        residue_n = 0
                        for atom_n, xyz in enumerate(model):
                            filout.write("ATOM {:>6d}  {:<3s} {:<3s} A {:>3d}     {:>7.3f} {:>7.3f} {:>7.3f}  1.00  1.00 {:>11s}\n".format(
                                atom_n, atom_names[index], residue_list[residue_n], residue_n, float(xyz[0]), float(xyz[1]), float(xyz[2]), atom_names[index][0]
                            ))
                            index += 1
                            if index % len(atom_names) == 0:
                                residue_n += 1
                                index = 0
                        filout.write("TER\nENDMDL\n")
            else:
                warning(Default.INPUT_ERROR, "Please select an output file type.")
                return
            current_models = objective_models

    def download_xyz(self, name, host, port, temperature_id, start = 0, end = 0, every = 1):
        #Retrieve input data
        trajectory = []
        structures = 'structures' + str(temperature_id)
        #Return XYZ data
        if not every == 1: query = {'step': {'$mod': [every, 0]}}
        else: query = {}
        if self.models == None:
            with pymongo.MongoClient(host, port) as mongo:
                self.models = mongo[name][structures].find(query, {'xyz': 1, '_id': 0})
        models = self.models.clone()
        models = models[start:end]
        for model in models:
            trajectory.append(model['xyz'])
        return trajectory

    def read_xyz(self, input_file):
        trajectory_data = []
        model_data = []
        self.pb_r['maximum'] = find_file_length(input_file)
        with open(input_file, 'r') as f:
            for i, line in enumerate(f):
                self.pb_r['value'] = i
                root.update_idletasks()
                if i == 0: continue
                elem = line.split()
                if len(elem) > 1: model_data.append([float(elem[1]), float(elem[2]), float(elem[3])])
                else:
                    model_data = []
                    trajectory_data.append(model_data)
        return trajectory_data

    def add_to_ramachandran(self, data, current_models = 0):
        psi, phi = [], []
        for i, model in enumerate(data):
            self.pb_r['value'] = current_models + i
            root.update_idletasks()
            model_length = len(model) - 3
            j = 2
            while j < model_length:
                a1 = model[j-3]
                a2 = model[j-1]
                a3 = model[j]
                a4 = model[j+1]
                a5 = model[j+2]
                if j % 5 == 1: psi.append(np.rad2deg(calculate_dihedral(a1, a3, a4, a5)))
                else:
                    a6 = model[j+4]
                    phi.append(np.rad2deg(calculate_dihedral(a2, a3, a4, a6)))
                if j % 5 == 1: j += 1
                else: j += 4
            # for j in xrange(1, model_length-10, 5):
            #     a1 = model[j]
            #     a2 = model[j+1]
            #     a3 = model[j+2]
            #     a4 = model[j+5]
            #     psi.append(np.rad2deg(
                # (a1, a2, a3, a4)))
            # for h in xrange(3, model_length-7, 5):
            #     b1 = model[h]
            #     b2 = model[h+3]
            #     b3 = model[h+4]
            #     b4 = model[h+5]
            #     phi.append(np.rad2deg(calculate_dihedral(b1, b2, b3, b4)))
        new_H, X, Y = np.histogram2d(phi, psi, bins = Default.BINS)
        self.H += new_H
        self.X = X
        self.Y = Y

    def produce_ramachandran(self):
        name, host, port= self.name_r.get(), self.host_r.get(), int(self.port_r.get())
        temperature = int(self.temperature_r.get())
        self.H = np.zeros((180, 180))
        self.models = None
        #Obtain data
        if self.v.get() == 1: #FROM FILE
            data = self.read_xyz(self.input_file.get() + '.xyz')
            self.add_to_ramachandran(data)
        elif self.v.get() == 2: #FROM DATABASE
            max_models = find_database_length(name, host, port, temperature)
            print "MAX models:", max_models
            self.pb_r['maximum'] = max_models
            root.update_idletasks()
            current_models = 0
            while current_models < max_models:
                objective_models = current_models + Default.BATCH_STEP
                if objective_models > max_models: objective_models = max_models
                print current_models, "->", objective_models

                data = self.download_xyz(name, host, port, temperature, start = current_models, end = objective_models)
                self.add_to_ramachandran(data, current_models)
                current_models = objective_models
        else:
            warning(Default.INPUT_ERROR, "No data source selected.")
            return
        
        X, Y = np.meshgrid(self.X, self.Y)
        plt.pcolormesh(X, Y, self.H)
        plt.colorbar()
        #phi, psi = self.add_to_ramachandran(data)
        # plt.hist2d(phi, psi, bins = 180)
        # plt.xlabel("PHI")
        # plt.ylabel("PSI")
        plt.show()

root = Tk()
style = ttk.Style()
style.theme_create( "MyStyle", parent="alt", settings={"TNotebook.Tab": {"configure": {"width": 1000},}})
style.theme_use("default")
utils_gui = Utils(root)
root.mainloop()