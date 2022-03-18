#!/bin/env python

import os, sys

official = False

if official:
    
    pileup = 'Collisions17_13TeV_PileUp_pileup_latest.txt'
    json = '../../../crab/JSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'

    os.system('pileupCalc.py --verbose -i '+json+' --inputLumiJSON '+pileup+' --calcMode true --minBiasXsec 66000 --maxPileupBin 100 --numPileupBins 100 puDown.root')
    os.system('pileupCalc.py --verbose -i '+json+' --inputLumiJSON '+pileup+' --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 puCentral.root')
    os.system('pileupCalc.py --verbose -i '+json+' --inputLumiJSON '+pileup+' --calcMode true --minBiasXsec 72400 --maxPileupBin 100 --numPileupBins 100 puUp.root')

else:
        
    s = ['JetHT']
        
    pileup = ['../../../crab/Pileup/'+s[0]+'.txt']
    json = ['../../../crab/Pileup/processedLumi/'+s[0]+'.json']

    for i, p in enumerate(pileup):
        
        os.system('pileupCalc.py --verbose -i '+json[i]+' --inputLumiJSON '+p+' --calcMode true --minBiasXsec 66000 --maxPileupBin 100 --numPileupBins 100 '+s[i]+'_puDown.root')
        os.system('pileupCalc.py --verbose -i '+json[i]+' --inputLumiJSON '+p+' --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 '+s[i]+'_puCentral.root')
        os.system('pileupCalc.py --verbose -i '+json[i]+' --inputLumiJSON '+p+' --calcMode true --minBiasXsec 72400 --maxPileupBin 100 --numPileupBins 100 '+s[i]+'_puUp.root')
