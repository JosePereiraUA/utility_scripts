import os

def zalign(target, state="default", export=0, unit="nm"):
    """
DESCRIPTION

    "zalign" performs a sequence aligment followed by structural superposition, using PyMol "align" function.
    Jose Pereira, Fev 2019 @ Universidade de Aveiro

USAGE
    
    zalign target_file [, state [, export [, unit ]]]

ARGUMENTS
    
    target_file = string: file containing the target molecule to be loaded

    state = string/int: state selected for aligment. "default" selects last frame. "all" shows RMSD values for all frames. n[int] selects frame n for aligment. {default: "default"}
    
    export = bool: if state is set to "all", export the results to a file ("rmsd_zalign.xvg"). {default: 0}
    
    unit = string: units to be used when displaying RMSD values. Supported units: "nm", "ang" {default: nm}
    
EXAMPLE

    zalign ../data/1i2t.gro, all, export

SEE ALSO

    align, intra_fit
    """

    #Set units
    if unit == "ang":
        d_unit = 1
    elif unit == "nm":
        d_unit = 10
    else:
        print "ERROR: \"%s\" unit unsuported." % unit
        return

    #Load export file
    if export and state == "all":
        data_file = open("rmsd_zalign.xvg", "w")

    #Load target structure
    if os.path.isfile(target) == False:
        print "ERROR: \"%s\" file not found." % target
        return
    cmd.load(target, "ref")

    #Define mobile structure as the first object loaded
    objects = cmd.get_object_list("(all)")

    if state == "default":
        target_state = cmd.count_states("(all)")
    elif state == "all":
        for target_state in range(1, cmd.count_states("(all)")+1):
            rmsd = float(cmd.align("%s and name CA" % (objects[0]), "ref and name CA", 999, mobile_state = target_state)[0])
            cmd.set("state", target_state)
            cmd.show_as("cartoon")
            print "%3d%9.4f" % (target_state, rmsd/d_unit)
            if export:
                data_file.write("%3d%9.4f\n" % (target_state, rmsd/d_unit))
        if export:
            print "Exported RMSD values to \"rmsd_zalign.xvg\""
            data_file.close()
        return 
    else:
        target_state = state
    rmsd = float(cmd.align("%s and name CA" % (objects[0]), "ref and name CA", 999, mobile_state = target_state)[0])
    cmd.set("state", target_state)
    cmd.show_as("cartoon")
    print "RMSD=%f" % (rmsd/d_unit)
    return 

cmd.extend("zalign", zalign)