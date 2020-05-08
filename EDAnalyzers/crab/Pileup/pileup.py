#!/usr/bin/env python

# Run on lxplus under CMSSW_11_0_X

import os, sys
from os.path import expanduser

os.environ["PATH"] += os.pathsep + expanduser("~") + '/.local/bin'
os.environ["PATH"] += os.pathsep + '/cvmfs/cms-bril.cern.ch/brilconda/bin'

dlumi = 'processedLumi/'
json = {}
json['JetHT'] = ['JetHT_Run2017B', 'JetHT_Run2017D', 'JetHT_Run2017E', 'JetHT_Run2017F']
json['ZeroBias'] = ['ZeroBias_Run2017B', 'ZeroBias_Run2017C', 'ZeroBias_Run2017D', 'ZeroBias_Run2017E', 'ZeroBias_Run2017F']

xsec = '69200' # Run 2

tmpdir = '/tmp/kskovpen/lumi/'
os.system('rm -rf '+tmpdir)
os.system('mkdir '+tmpdir)

for k, v in json.iteritems():
    
    print k
    
    jfiles = ' '.join([dlumi+s+'.json' for s in v])

    os.system('rm -rf '+dlumi+k+'.json')
    os.system('mergeJSON.py '+jfiles+' --output '+dlumi+k+'.json')
    
    os.system('brilcalc lumi --xing --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -i '+dlumi+k+'.json -o '+tmpdir+k+'.csv')
    os.system('makePileupJSON.py '+tmpdir+k+'.csv '+k+'.txt')
    os.system('pileupCalc.py -i '+dlumi+k+'.json --inputLumiJSON '+k+'.txt'+' --calcMode true --minBiasXsec '+xsec+' --maxPileupBin 100 --numPileupBins 100 '+k+'.root')
    
os.system('rm -rf '+tmpdir)
