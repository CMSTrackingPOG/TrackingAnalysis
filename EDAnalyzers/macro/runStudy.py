#!/usr/bin/env python

import os, sys

from optparse import OptionParser

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]
        
    usage = "usage: %prog [options]\n Run the analysis of histograms and produce plots"
    
    parser = OptionParser(usage)
    parser.add_option("-p","--pv",action='store_true',help="Run PV study [default: %default]")
    parser.add_option("--ippv",action='store_true',help="Run IPPV measurement[default: %default]")
    parser.add_option("--ipbs",action='store_true',help="Run IPBS measurement [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options                                        

if __name__ == '__main__':
    
    options = main()

    jobDir = 'jobs/'

    inputMC = ['MinBiasTune4C8TeVpythia8_Summer12DR53XPUS10START53V7Av1AODSIM']
    inputData = ['MinimumBias_Run2012A22Jan2013v1AOD','MinimumBias_Run2012B22Jan2013v1AOD',\
    'MinimumBias_Run2012C22Jan2013v1AOD','MinimumBias_Run2012D22Jan2013v1AOD']

    for f in range(len(inputMC)): inputMC[f] = jobDir+inputMC[f]+'.root'
    for f in range(len(inputData)): inputData[f] = jobDir+inputData[f]+'.root'

    fMC = jobDir+'mc.root'
    fData = jobDir+'data.root'

    if not os.path.isfile(fMC):
        filesMC = " ".join(inputMC)
        os.system('hadd -f '+fMC+' '+filesMC)

    if not os.path.isfile(fData):    
        filesData = " ".join(inputData)
        os.system('hadd -f '+fData+' '+filesData)

    picdir = 'pics'
    if os.path.isdir(picdir):
        os.system("rm -rf "+picdir)
    os.system("mkdir "+picdir)

    if options.pv:
        print "Run PV study"
        os.system('python pvStudy.py --data='+fData+' --mc='+fMC)

    if options.ippv:
        print "Run IPPV study"
        os.system('python ipStudy.py --data='+fData+' --mc='+fMC+' --type=pv')

    if options.ipbs:
        print "Run IPBS study"
        os.system('python ipStudy.py --data='+fData+' --mc='+fMC+' --type=bs')
        
