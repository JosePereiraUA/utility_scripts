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
        for row in csv_reader:
            row = filter(lambda a: a != '', [a.decode('utf-8').encode('ascii', 'ignore') for a in row if not row[0] == ''])
            if len(row) > 0:
                data.append(row)

    # Format data in JSON
    sc_hb_lib = {}
    cur_aa = None
    for (aa, ad, base, charged) in data:
        print aa, ad, base, charged
        if aa != cur_aa:
            sc_hb_lib[aa] = {'A': [], 'D': []}
        sc_hb_lib[aa][ad].append((base, charged))
        cur_aa = aa
    
    # Export JSON file
    with open(out_file, 'w') as fout:
        json.dump(sc_hb_lib, fout, sort_keys=True, indent=2, separators=(',', ': '))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Reads a Rotamer Library in CSV format and exports it in JSON format')
    parser.add_argument('-f', '--file', type=str, help='CSV rotamer data file (Required)', required = True, metavar = '')
    parser.add_argument('-o', '--out', type=str, help='Output JSON file name (Default: sc_hb.json)', default = 'sc_hb.json', metavar = '')
    args = parser.parse_args()

    main(args.file, args.out)