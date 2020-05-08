#!/usr/bin/env python

# Run on lxplus under CMSSW_11_0_X

import os, sys
from os.path import expanduser

os.environ["PATH"] += os.pathsep + expanduser("~") + '/.local/bin'
os.environ["PATH"] += os.pathsep + '/cvmfs/cms-bril.cern.ch/brilconda/bin'

dlumi = 'processedLumi/'
json = ['JetHT_Run2017B', 'JetHT_Run2017D', 'JetHT_Run2017E', 'JetHT_Run2017F', \
'ZeroBias_Run2017B', 'ZeroBias_Run2017C', 'ZeroBias_Run2017D', 'ZeroBias_Run2017E', 'ZeroBias_Run2017F']

xsec = '69200' # Run 2

tmpdir = '/tmp/kskovpen/lumi/'
os.system('rm -rf '+tmpdir)
os.system('mkdir '+tmpdir)

for j in json:
    
    print j
    
    os.system('brilcalc lumi --xing --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i '+dlumi+j+'.json -o '+tmpdir+j+'.csv')
    os.system('makePileupJSON.py --csvInput '+tmpdir+j+'.csv '+j+'.txt')
    os.system('pileupCalc.py -i '+dlumi+j+plumi+' --inputLumiJSON '+j+'.txt'+' --calcMode true --minBiasXsec '+xsec+' --maxPileupBin 100 --numPileupBins 100 '+j+'.root')
    
os.system('rm -rf '+tmpdir)
