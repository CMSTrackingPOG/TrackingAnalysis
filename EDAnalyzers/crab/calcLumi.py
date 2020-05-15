#!/bin/env python

import os

ppath = [os.path.expanduser('~')+'/.local/bin', '/cvmfs/cms-bril.cern.ch/brilconda/bin'] 
os.environ["PATH"] += os.pathsep + os.pathsep.join(ppath)

fjson = 'JSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'

os.system('brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb --byls -i '+fjson+' -o lumi.csv')
