#!/usr/bin/env python

# Run on lxplus

import os, sys
from os.path import expanduser

os.environ["PATH"] += os.pathsep + expanduser("~") + '/.local/bin'
os.environ["PATH"] += os.pathsep + '/cvmfs/cms-bril.cern.ch/brilconda/bin'

dlumi = '../'
plumi = '/results/processedLumis.json'
json = ['crab_JetHT_Run2017B__1', 'crab_JetHT_Run2017D__1', 'crab_JetHT_Run2017E__1', 'crab_JetHT_Run2017F__1', \
'crab_ZeroBias_Run2017B__1', 'crab_ZeroBias_Run2017C__1', 'crab_ZeroBias_Run2017D__1', 'crab_ZeroBias_Run2017E__1', 'crab_ZeroBias_Run2017F__1']

xsec = '69200' # Run 2

for j in json:
    
    print j
    
    os.system('brilcalc lumi --xing --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i '+dlumi+j+plumi+' -o '+j+'.csv')
    os.system('makePileupJSON.py --csvInput '+j+'.csv '+j+'.txt')
    os.system('pileupCalc.py -i '+dlumi+j+plumi+' --inputLumiJSON '+j+'.txt'+' --calcMode true --minBiasXsec '+xsec+' --maxPileupBin 100 --numPileupBins 100 '+j+'.root')
