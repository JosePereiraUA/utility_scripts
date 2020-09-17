# -*- coding: utf-8 -*-
"""
     Created   José Pereira     April      2019
Last updated   José Pereira     April      2019
"""
import csv
import json
import argparse

def main(csv_file, out_file):

    # Extract data from CSV file
    data = []
    with open(csv_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            row = filter(lambda a: a != '', [a.decode('utf-8').encode('ascii', 'ignore') for a in row if not row[0] == ''])
            if len(row) > 0:
                data.append(row)

    # Format data in JSON
    rot_lib = {}
    cur_aa = None
    n_chis = 0
    for line in data:
        if len(line) == 2:
            cur_aa = line[0]
            n_chis = int(line[1])
            rot_lib[cur_aa] = {'w': [], 'rot': [], 'range': []}
        else:
            rot_lib[cur_aa]['w'].append(float(line[1]))
            rot_lib[cur_aa]['rot'].append([float(chi) for chi in line[2:2+n_chis]])
            rot_lib[cur_aa]['range'].append([float(r) for r in line[2+n_chis:]])

    # Normalize weights
    for aa in rot_lib.keys():
        s = sum(rot_lib[aa]['w'])
        rot_lib[aa]['w'] = [w/s for w in rot_lib[aa]['w']]

    # Export JSON file
    with open(out_file, 'w') as fout:
        json.dump(rot_lib, fout, sort_keys=True, indent=2, separators=(',', ': '))
        
    print "All done! Rotamer library exported to %s" % (out_file)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Reads a Rotamer Library in CSV format and exports it in JSON format')
    parser.add_argument('-f', '--file', type=str, help='CSV rotamer data file (Required)', required = True, metavar = '')
    parser.add_argument('-o', '--out', type=str, help='Output JSON file name (Default: rot_lib.json)', default = 'rot_lib.json', metavar = '')
    args = parser.parse_args()

    main(args.file, args.out)