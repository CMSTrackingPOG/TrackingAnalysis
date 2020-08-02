#!/bin/env python

import os, sys, json

dsaod = '/SingleNeutrino/RunIISummer19UL17RECO-FEVTDEBUG_106X_mc2017_realistic_v6-v2/AODSIM'
dsraw = '/SingleNeutrino/RunIISummer19UL17DIGI-FEVTDEBUG_106X_mc2017_realistic_v6-v2/GEN-SIM-DIGI-RAW'

filesaod = os.popen('dasgoclient --query="file dataset='+dsaod+'"').read().splitlines()
filesraw = os.popen('dasgoclient --query="file dataset='+dsraw+'"').read().splitlines()

nmax = -1

print 'Read AOD DBS data -----'

lumiaod = {}
siteaod = {}
aod = {}
for i, f in enumerate(filesaod):
    if i >= nmax and nmax >= 0: break
    print f
    lumiaod[f] = os.popen('dasgoclient --query="lumi file='+f+'"').read().splitlines()
    siteaod[f] = os.popen('dasgoclient --query="site file='+f+'"').read().splitlines()
    aod[f] = [lumiaod[f], siteaod[f]]

print 'Read RAW DBS data -----'
    
lumiraw = {}
siteraw = {}
raw = {}
for i, f in enumerate(filesraw):
    if i >= nmax and nmax >= 0: break
    print f
    lumiraw[f] = os.popen('dasgoclient --query="lumi file='+f+'"').read().splitlines()
    siteraw[f] = os.popen('dasgoclient --query="site file='+f+'"').read().splitlines()
    raw[f] = [lumiraw[f], siteraw[f]]

with open('aod.json', 'w') as write_file:
    json.dump(aod, write_file, indent=2)
with open('raw.json', 'w') as write_file:
    json.dump(raw, write_file, indent=2)

print 'Find related AOD and RAW files -----'

with open('match.json', 'w') as write_file:

    match = {}
    
    for ak, av in lumiaod.iteritems():
        for rk, rv in lumiraw.iteritems():
            for iav in av:
                for irv in rv:
                    for iaf in iav.replace('[','').replace(']','').split(','):
                        for irf in irv.replace('[','').replace(']','').split(','):
                            if iaf == irf:
                                match[iaf] = [ak, siteaod[ak], rk, siteraw[rk]]
                                
    json.dump(match, write_file, indent=2)
