# __ Application information
# TITLE ........... MD Maker (No-GUI)
# AUTHOR .......... Jose Pereira
# VERSION ......... 1.3
# LAST MODIFIED ... May 2020
#
# Version history:
# V1.0 - Initial version, performs MD simulations
# V1.1 - Added Step class, allowing the manipulation of several variables
# V1.2 - Added support for simulations with ligands (with an extra .itp file)
# V1.3 - Added support for PMF simulations

import os
import sys
sys.path.append('/home/jpereira/Desktop/scripts/ze_utils')
from ze_utils.molecule_manipulation import *

class Step:

    def __init__(self, title, time_step = 0.0, n_steps = 0, print_every = 0, overwrite = False):
        self.title = title
        self.mdp   = "/home/jpereira/Desktop/utility_scripts/mdps/" + title + '.mdp'
        self.time_step    = time_step
        self.n_steps      = n_steps
        self.print_every  = print_every
        self.overwrite    = overwrite
        if self.title != 'solvate':
            self.config_mdp()


    def config_mdp(self):
        lines = []
        with open(self.mdp, 'r') as file_in:
            for line in file_in:
                if line.startswith("dt"):
                    line = line[:-2] + " " + str(self.time_step) + "\n"
                if line.startswith("nsteps"):
                    line = line[:-2] + " " + str(self.n_steps) + "\n"
                elif line.startswith("nstxout") or line.startswith("nstenergy") or line.startswith("nstlog") or line.startswith("nstxout-compressed"):
                    line = line[:-2] + " " + str(self.print_every) + "\n"
                lines.append(line)

        with open(os.path.split(self.mdp)[-1], 'w') as file_out:
            for line in lines:
                file_out.write(line)


    """
    Run GROMPP/MDRUN to perform simulation. Automatically extract trajectory and
    last frames to specific files. Returns the last frame file name.
    """
    def run(self, input_file, active_site = None):

        if os.path.isdir(os.getcwd() + '/' + self.title) and not self.overwrite:
            print("BYPASSING STEP", self.title)
            return self.title + '/' + self.title + '.gro'

        # 1. Config MDP
        self.config_mdp()

        # 2.
        print("Pre-processing topology with GROMPP ...")
        pmf_entry = ""
        if self.title.endswith("_pmf"):
            os.system("echo q | gmx make_ndx -f {i}".format(i = input_file))
            for line in active_site.split("\n"):
                os.system("echo {text} >> index.ndx".format(text = line))
            pmf_entry = "-n index.ndx -r {i} -t {cpt_file}"\
                .format(i = input_file, cpt_file = input_file[:-4] + ".cpt")
        os.system("gmx grompp -f {mdp} -c {i} -p topol.top -o {t}.tpr -maxwarn 1 {_pmf}"\
            .format(mdp = os.path.split(self.mdp)[-1], i = input_file,
                    t = self.title, _pmf = pmf_entry))

        # 3.
        print("Running {title} ...".format(title = self.title))
        pmf_entry = ""
        if self.title.endswith("_pmf"):
            pmf_entry = "-pf {t}_F.xvg -px {t}_X.xvg".format(t = self.title)
        os.system("gmx mdrun -deffnm {t} -nt 1 -s {t}.tpr -v {_pmf}"\
            .format(t = self.title, _pmf = pmf_entry))

        # 4.
        print("Extracting trajectory's last frame to {title}.pdb ..."\
            .format(title = self.title))
        last_frame = self.verify_output()
        os.system("echo '0 \n 0' | gmx trjconv -f {title}.trr -o {title}.pdb -s {title}.tpr -pbc mol -conect -center -b {frame}"\
            .format(title = self.title, frame = last_frame))

        # 5.
        print("Extracting trajectory to {title}_traj.pdb ...".format(title = self.title))
        os.system("echo '0 \n 0' | gmx trjconv -f {title}.trr -o {title}_traj.pdb -s {title}.tpr -pbc mol -conect -center"\
            .format(title = self.title))

        # 6.
        print("Storing data in {folder} ...".format(folder = '/' + self.title))
        self.store_data()

        return self.title + '/' + self.title + '.gro'


    """
    1) Verifies if the minimization steps converged in the input number of steps;
    2) Returns the frame number/time of the last strucuture in the trajectory.
    """
    def verify_output(self):

        with open(self.title + '.log', 'r') as file_in:
            for line in file_in:
                if self.title.startswith('minimization'):
                    if line.startswith("Steepest Descents did not converge"):
                        print("Convergence Error\n", line)
                        self.store_data()
                        self.clean()
                        exit(1)
                    elif line.startswith("Steepest Descents converged"):
                        elems = line.split()
                        if elems[4] == "machine":
                            return elems[7]
                        else:
                            return elems[8]
                elif self.title.startswith('steered_md_pmf'):
                    gather = False
                    if gather:
                        gather = False
                        frame = float(self.time_step) * float(line.split()[1])
                    if line.startswith("Step"):
                        gather = True
                    if line.startswith("Fatar error:"):
                        return frame

            return float(self.n_steps) * float(self.time_step)


    """
    Creates a new folder and moves all necessary files for correct data keeping.
    """
    def store_data(self, folder = None):

        if folder == None:
            folder = self.title

        os.system("mkdir {folder}".format(folder = folder))
        os.system("mv {title}* mdout.mdp {folder}"\
            .format(title = self.title, folder = folder))


    """
    Delete unnecessary files.
    """
    def clean(self):

        print("Cleaning unnecessary files ...")
        os.system("rm -rf *#")
        print("All tasks successefully performed.")


class PMFStep(Step):

    def __init__(self, title, time_step = 0.0, n_steps = 0, print_every = 0, d_step_size = 0, n_d_steps = 0, overwrite = False):
        self.title = title
        self.mdps   = ["/home/jpereira/Desktop/utility_scripts/mdps/equilibration_npt_pmf.mdp", "/home/jpereira/Desktop/utility_scripts/mdps/collection_pmf.mdp"]
        self.time_step    = time_step
        self.n_steps      = n_steps
        self.print_every  = print_every
        self.d_step_size  = d_step_size
        self.n_d_steps    = n_d_steps
        self.tolerance    = 1e-2
        self.overwrite    = overwrite

    """
    """
    def config_mdps(self):
        for mdp in self.mdps:
            lines = []
            with open(mdp, 'r') as file_in:
                for line in file_in:
                    if line.startswith("dt"):
                        line = line[:-2] + " " + str(self.time_step) + "\n"
                    if line.startswith("nsteps"):
                        line = line[:-2] + " " + str(self.n_steps) + "\n"
                    elif line.startswith("nstxout") or line.startswith("nstenergy") or line.startswith("nstlog") or line.startswith("nstxout-compressed"):
                        line = line[:-2] + " " + str(self.print_every) + "\n"
                    lines.append(line)

            print("Writting to ", os.path.split(mdp)[-1])
            with open(os.path.split(mdp)[-1], 'w') as file_out:
                for line in lines:
                    file_out.write(line)


    def read_X_file(self, input_file):
        # FIX FOR THE CASE WHERE THE FILE Has a SINGLE EXTRA X (DIDNT FINISH ON THE END)
        self.x, self.y = [], []
        with open(input_file, "r") as xvg:
            for line in xvg:
                if line.startswith(("@", "#")):
                    continue
                elem = line.split()
                if len(elem) < 2:
                    continue
                self.x.append(float(elem[0]))
                self.y.append(float(elem[1]))


    """
    Run GROMPP/MDRUN to perform simulation. Automatically extract trajectory and
    last frames to specific files. Returns the last frame file name.
    """
    def run(self, input_file, active_site = None):

        # 1. Config MDP
        self.config_mdps()

        if os.path.exists("tpr-files.dat") and self.overwrite:
            os.system("rm -r tpr-files.dat")
            os.system("touch tpr-files.dat")
        if os.path.exists("pullf-files.dat") and self.overwrite:
            os.system("rm -r pullf-files.dat")
            os.system("touch pullf-files.dat")

        # 2. Create the correct input file, at the desired distance.
        self.read_X_file(input_file[:-4] + "_X.xvg")
        for index in range(self.n_d_steps):
            distance = self.y[0] + index * self.d_step_size
            print("Looking for ", distance)
            found = False
            for x, y in zip(self.x, self.y):
                if y > distance - self.tolerance and y < distance + self.tolerance:
                    print("FOUND CORRECT DISTANCE {d} at time {x}".format(d = distance, x = x))
                    found = True
                    time = x
                    break
            if not found:
                print("Distance {d} was not found with the given steered MD trajectory (with tolerance = {t})."\
                    .format(d = distance, t = self.tolerance))

            _title = "{title}_{distance:.3f}"\
                .format(title = self.title, distance = distance)

            if os.path.isdir(_title):
                if self.overwrite:
                    os.system("rm -rf {title}".format(title = _title))
                else:
                    continue

            print("Extracting trajectory's frame at time {time} to {title}_{distance:.2f}.pdb ..."\
                .format(time = time, title = self.title, distance = distance))
            os.system("echo '1 \n 0' | gmx trjconv -f {i}.trr -o {title}.pdb -s {i}.tpr -pbc mol -conect -center -dump {time}"\
                .format(i = input_file[:-4], title = _title, distance = distance, time = time))


            # 3. NPT Equilibration
            print("Pre-processing topology with GROMPP (FOR NPT EQUILIBRATION) ...")
            os.system("echo q | gmx make_ndx -f {i}.pdb".format(i = _title))
            for line in active_site.split("\n"):
                os.system("echo {text} >> index.ndx".format(text = line))
            pmf_entry = "-n index.ndx -r {t}.pdb".format(t = _title)
            os.system("gmx grompp -f {mdp} -c {t}.pdb -p topol.top -o {t}_NPT.tpr -maxwarn 1 {_pmf}"\
                .format(mdp = os.path.split(self.mdps[0])[-1], t = _title,
                    _pmf = pmf_entry))

            print("Running {title}_NPT ...".format(title = _title))
            pmf_entry = "-pf {t}_NPT_F.xvg -px {t}_NPT_X.xvg".format(t = _title)
            os.system("gmx mdrun -deffnm {t}_NPT -v {_pmf}"\
                .format(t = _title, _pmf = pmf_entry))

            # 4. Collection
            print("Pre-processing topology with GROMPP (FOR COLLECTION) ...")
            pmf_entry = "-n index.ndx -r {t}_NPT.gro -t {t}_NPT.cpt"\
                .format(t = _title)
            os.system("gmx grompp -f {mdp} -c {t}_NPT.gro -p topol.top -o {t}_collection.tpr -maxwarn 1 {_pmf}"\
                .format(mdp = os.path.split(self.mdps[1])[-1], t = _title, _pmf = pmf_entry))

            print("Running {title}_collection ...".format(title = _title))
            pmf_entry = "-pf {t}_collection_F.xvg -px {t}_collection_X.xvg"\
                .format(t = _title)
            os.system("gmx mdrun -deffnm {t}_collection -v {_pmf}"\
                .format(t = _title, _pmf = pmf_entry))

            print("Extracting trajectory to {title}_traj.pdb ...".format(title = _title))
            os.system("echo '1 \n 17' | gmx trjconv -f {title}_collection.trr -o {title}_collection_traj.pdb -s {title}_collection.tpr -pbc mol -conect -center"\
                .format(title = _title))

            print("Storing data in {folder} ...".format(folder = '/' + _title))
            self.store_data(_title)

            with open("tpr-files.dat", "a") as tpr_files:
                tpr_files.write("{title}/{title}_collection.tpr\n".format(title = _title))
            with open("pullf-files.dat", "a") as tpr_files:
                tpr_files.write("{title}/{title}_collection_F.xvg\n".format(title = _title))

        # 5. Analysing data with WHAM
        os.system("gmx wham -it tpr-files.dat -if pullf-files.dat -o -hist -unit kCal")
        os.system("mkdir pmf_results")
        os.system("cp profile.xvg histo.xvg pmf_results")

    """
    Creates a new folder and moves all necessary files for correct data keeping.
    """
    def store_data(self, folder = None):

        if folder == None:
            folder = self.title

        os.system("mkdir {folder}".format(folder = folder))
        os.system("mv {folder}* mdout.mdp {folder}"\
            .format(title = self.title, folder = folder))


# ------------------------------------------------------------------------------


class Default():
    """
    Holds default values for various variables. This defines how the process
    will work. When defining a Step to use PMF MDP's and settings, add the
    "_pmf" suffix to end of the Step title.
    """

    # Define inputs
    MDP_DIR        = "/home/jpereira/Desktop/utility_scripts/mdps"
    INPUT_FILE     = "md.pdb"
    INPUT_ITP      = "ligand.itp"

    # Define intermediary files
    ALL_H_FILE     = INPUT_FILE[:-4] + '_all_h.pdb'
    BOXED_FILE     = INPUT_FILE[:-4] + '_boxed.pdb'
    SOLVATED_FILE  = INPUT_FILE[:-4] + '_solvt.pdb'
    IONIZED_TPR    = INPUT_FILE[:-4] + '_ion.tpr'
    IONIZED_FILE   = INPUT_FILE[:-4] + '_ion.pdb'

    # Define simulation variables
    BOX_SHAPE      = "cubic"
    BOX_SIZE       = 8 # in nm
    WATER_TYPE     = "spc216"
    LIGAND_ONLY    = False # Used in topology creation when molecule being simulated is NOT a protein

    # Define steps
    #                           Name                       Time-step   N-steps   Print-every      Overwrite
    # STEPS          = [
    #                        Step('minimization_1',          0.002,      50000,    1000,            False),
    #                        Step('solvate'),
    #                        Step('minimization_2',          0.002,      50000,    1000,            False),
    #                        Step('equilibration_nvt',       0.002,      250000,   1000,            False),
    #                        Step('equilibration_npt',       0.002,      250000,   1000,            False),
    #                        Step('distance_lock_pmf',       0.002,      20000,    1000,            False),
    #                        Step('steered_md_pmf',          0.002,      500000,   10000,           True),
    #                     PMFStep('umbrella_sampling_pmf',   0.002,      500000,   10000, 0.05, 30, True)
    #                 ]
    STEPS          = [
                           Step('minimization_1',          0.002,      50000,    10000,            False),
                           Step('solvate'),
                           Step('minimization_2',          0.002,      50000,    10000,            False),
                           Step('equilibration_nvt',       0.002,      250000,   10000,            False),
                           Step('equilibration_npt',       0.002,      250000,   10000,            False),
                           Step('collection',              0.002,      5000000,  10000,            False)
                    ]

    # Define PMF-specific variables
    CUTOFF     = 8.0 # Cut off used when defining the active site (if ligand is present)


# ------------------------------------------------------------------------------


class MDMaker:

    def __init__(self, default = None):
        if default == None:
            self.default = Default()
            if os.path.exists(self.default.INPUT_ITP):
                ligand_name      = self.extract_ligand_name()
                self.active_site = self.get_active_site(ligand_name)
        else:
            self.default = default

        self.process()


    """
    """
    def extract_ligand_name(self):
        with open(self.default.INPUT_ITP, "r") as ligand_itp:
            listening = False
            index     = 0
            for line in ligand_itp:
                if listening:
                    if index == 1:
                        return line.split()[0]
                    else:
                        index += 1
                if line.startswith("[ moleculetype ]"):
                    listening = True

    """
    """
    def get_active_site(self, ligand_name):
        print("Defining active site (this may take a minute or two ...)")
        m = Molecule(self.default.INPUT_FILE, read_hetatms = True)
        ligand_atoms = m.get_atoms_from_residue_names([ligand_name])[1]
        neighbours = m.get_atoms_within_of(self.default.CUTOFF, ligand_atoms)

        residues = []
        for i in neighbours:
            residues.append(m.atoms[i].res_index)
        residues_list = list(np.unique(residues))

        t = "[ ACTIVESITE ]\n"
        for i in residues_list:
            t += "%d " % (i)
        return t[:-1]


    """
    Process control. Generate data.
    """
    def process(self):

        # Verify that all necessary files exist
        self.verify_input()

        # Correctly compile all necessary topology files
        self.create_topology_files(self.default.INPUT_FILE)

        # Add a bounding box to the input file using gmx editconf
        self.add_bounding_box(self.default.BOX_SIZE)

        previous_step_file = self.default.BOXED_FILE
        for step in self.default.STEPS:
            if step.title == 'solvate':
                previous_step_file = self.solvate(previous_step_file)
            else:
                previous_step_file = step.run(previous_step_file,
                    active_site = self.active_site)

        self.clean()


    """
    1) Verifies the existance of the necessary .mdp files;
    2) Verifies the existance of previous MD data in the current work directory.
    """
    def verify_input(self):

        for step in self.default.STEPS:
            if step.title.startswith('solvate') or type(step) == PMFStep:
                continue
            file_name = self.default.MDP_DIR + '/' + step.title + '.mdp'
            if not os.path.isfile(file_name):
                print("Missing File", "{file} not found".format(file = file_name))
                exit(1)

        for step in self.default.STEPS:
            if (os.path.isdir(os.getcwd() + '/' + step.title) and step.overwrite) == True:
                os.system("rm -rf {title}".format(title = os.getcwd() + '/' + step.title))


    def create_topology_files(self, input_file, output_file = None):
        """
        1) Creates a topology file using GMX PDB2GMX;
        2) Edits the produced topology file to include solvent and ions forcefield.
        """

        if output_file == None:
            output_file = self.default.ALL_H_FILE

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

        print("Creating topology files ...")
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
        if os.path.exists(self.default.INPUT_ITP) and not self.default.LIGAND_ONLY:
            itp_exists = True
            itp = open(self.default.INPUT_ITP, "r")
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
            with open(self.default.INPUT_ITP[:-4] + "_auto.itp", "w") as new_itp_fh:
                for line in new_itp:
                    new_itp_fh.write(line)

            # 4) Create temporary PDB of just protein for pdb2gmx topology
            input_name = "%s_auto_no_ligand.pdb" % (input_file[:-4])
            with open(input_name, "w") as tmp:
                for line in protein_lines:
                    tmp.write(line)

            # 5) pdb2gmx on protein only
            output_name = "%s_auto_no_ligand_corrected.pdb" % (input_file[:-4])
            os.system("echo '6 \n 7 \n 0' | gmx pdb2gmx -f {input_file} -o {output_file} -his -ignh"\
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
            with open(output_name, "r") as output_name_in:
                for line in output_name_in:
                    if line == "TER\n":
                        for ligand_line in ligand_lines:
                            joined_pdb.append(ligand_line)
                    joined_pdb.append(line)

            with open(output_file, "w") as joined_out:
                for line in joined_pdb:
                    joined_out.write(line)

        elif not self.default.LIGAND_ONLY:
            os.system("echo '6 \n 7 \n 0' | gmx pdb2gmx -f {input_file} -o {output_file} -his"\
                .format(input_file = input_file, output_file = output_file))

        if self.default.LIGAND_ONLY:
            os.system("cp %s topol.top" % (self.itp()))
            os.system("cp %s %s" % (input_file, output_file))

        print("Adding solvent and ligand data to topology file ...")
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
                    if self.default.LIGAND_ONLY:
                        lines.append('; Include forcefield parameters\n#include "amber99sb-ildn.ff/forcefield.itp"\n\n'.encode())
                    lines.append(line)
                elif len(elem) == 2 and elem[1][-15:-1] == "forcefield.itp":
                    lines.append(line)
                    if itp_exists:
                        lines.append('#include "atomtypes.itp"\n'.encode())
                elif len(elem) == 3 and elem[1] == 'system':
                    lines.append('; Include water topology\n#include "amber99sb-ildn.ff/tip3p.itp"\n\n; Include topology for ions\n#include "amber99sb-ildn.ff/ions.itp"'.encode())
                    if itp_exists:
                        lines.append('\n\n; Include ligand topology\n#include "%s"'.encode() % (self.default.INPUT_ITP[:-4] + "_auto.itp").encode())
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

        # if itp_exists:
        #     os.system("cp %s %s" % (input_file, output_file))
        #     os.system("rm -rf %s %s" % (input_name, output_name))


    """
    Adds a bounding box to the input structure, of the requested size, using GMX EDITCONF
    """
    def add_bounding_box(self, size, input_file = None, output_file = None,
        shape = None):

        if input_file == None:
            input_file = self.default.ALL_H_FILE
        if output_file == None:
            output_file = self.default.BOXED_FILE
        if shape == None:
            shape = self.default.BOX_SHAPE

        print("Boxing protein with {shape} box of size {size} ..."\
            .format(shape = shape, size = size))
        os.system("gmx editconf -f {input_file} -o {output_file} -c -box {size} -bt {shape}"\
            .format(input_file = input_file, output_file = output_file, size = size, shape = shape))


    """
    Adds solvent molecules to the simulation box, correctly labeling the added
    molecules in the topology file.
    """
    def solvate(self, input_file, water_type = None, solvated_file = None,
        tpr_file = None, output_file = None):

        if water_type == None:
            water_type = self.default.WATER_TYPE
        if solvated_file == None:
            solvated_file = self.default.SOLVATED_FILE
        if tpr_file == None:
            tpr_file = self.default.IONIZED_TPR
        if output_file == None:
            output_file = self.default.IONIZED_FILE

        print("Adding {water_type} water molecules to simulation ..."\
            .format(water_type = water_type))
        os.system("gmx solvate -cp {input_file} -cs {water_type}.gro -o {output_file} -p topol.top"\
                .format(input_file = input_file, water_type = water_type, output_file = solvated_file))

        print("Adding ions to system ...")
        os.system("gmx grompp -f {mdp} -c {input_file} -p topol.top -o {output_file}"\
            .format(mdp = self.default.MDP_DIR + '/ions.mdp', input_file = solvated_file, output_file = tpr_file))

        if self.default.LIGAND_ONLY:
            os.system("echo '4' | gmx genion -s {input_file} -p topol.top -o {output_file} -pname NA -nname CL -neutral"\
                .format(input_file = tpr_file, output_file = output_file))
        elif os.path.exists(self.default.INPUT_ITP):
            os.system("echo '15' | gmx genion -s {input_file} -p topol.top -o {output_file} -pname NA -nname CL -neutral"\
                .format(input_file = tpr_file, output_file = output_file))
        else:
            os.system("echo '13' | gmx genion -s {input_file} -p topol.top -o {output_file} -pname NA -nname CL -neutral"\
                .format(input_file = tpr_file, output_file = output_file))

        return output_file


    """
    Delete unnecessary files.
    """
    def clean(self):

        print("Cleaning unnecessary files ...")
        os.system("rm -rf *#")
        print("All tasks successefully performed.")


if __name__ == "__main__":
    MDMaker()