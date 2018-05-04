#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pymongo
from matplotlib import pyplot as plt


if __name__ == '__main__':

    db_host = '192.168.179.114'
    db_port = 27017
    db_name = '1ctf_19_4_2018'
    db_temp = 6
    db_from = 30000
    alpha   = 0.1
    beta    = 0.015

    with pymongo.MongoClient(db_host, db_port) as server:
        db = server[db_name]['structures{temp}'.format(temp = db_temp)]
        energies = db.find({}, {'eBond': 1, 'eH': 1, 'energy': 1, 'eSol': 1, 'eIJ': 1, '_id': 0})\
            .skip(db_from)

    for energy in energies:
        energy['eBond'] = energy['eBond'] * 50
        energy['eSol']  = energy['eSol']  * alpha * beta
        energy['eIJ']   = energy['eIJ']   * alpha
        energy['eH']    = energy['eH']

    for key in energies.distinct():
        print key
