#!/usr/bin/env python

import os, sys
import ROOT

from optparse import OptionParser

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]
        
    usage = "usage: %prog [options]\n Run the analysis of histograms and produce plots"
    
    parser = OptionParser(usage)
    parser.add_option("--pv",action='store_true',help="Run PV study [default: %default]")
    parser.add_option("--ippv",action='store_true',help="Run IPPV measurement[default: %default]")
    parser.add_option("--ipbs",action='store_true',help="Run IPBS measurement [default: %default]")
    parser.add_option("--withbs",action='store_true',help="Use WithBS samples [default: %default]")
    parser.add_option("--qcd",action='store_true',help="Use QCD events [default: %default]")
    parser.add_option("--reweight",action='store_true',help="Produce the reweighting file [default: %default]")
    parser.add_option("--jobs",default="jobs_dps/",help="Directory with input files [default: %default]")
    parser.add_option("--clean",action='store_true',help="Clean merged files [default: %default]")
    parser.add_option("--fit",action='store_true',help="Extract IP resolution from fits [default: %default]")
    parser.add_option("--method", default="fwhm", help="Method to extract IP resolution (fwhm or fit) [default: %default]")
    parser.add_option("--selection", default="", help="Selection bin for IP resolution [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options                                        

if __name__ == '__main__':
    
    options = main()

    jobDir = options.jobs

    inputMC = ['SingleNeutrino_RunIISummer19UL17RECO']
    inputData = ['ZeroBias_Run2017B','ZeroBias_Run2017C',\
    'ZeroBias_Run2017D','ZeroBias_Run2017E','ZeroBias_Run2017F']
    if options.qcd:
        inputMC = ['QCDPt15to7000TuneCP5Flat13TeVpythia8_RunIISummer19UL17RECO']
        inputData = ['JetHT_Run2017B','JetHT_Run2017D',\
        'JetHT_Run2017E','JetHT_Run2017F']
        
    if options.withbs:
        inputMC = ['SingleNeutrino_RunIISummer19UL17RECOwithBS']
        inputData = ['ZeroBias_Run2017CwithBS',\
        'ZeroBias_Run2017DwithBS','ZeroBias_Run2017FwithBS']

    for f in range(len(inputMC)): inputMC[f] = jobDir+inputMC[f]+'.root'
    for f in range(len(inputData)): inputData[f] = jobDir+inputData[f]+'.root'

    fMC = jobDir+'mc.root'
    fData = jobDir+'data.root'

    if not os.path.isfile(fMC) or options.clean:
        filesMC = " ".join(inputMC)
        os.system('hadd -ff '+fMC+' '+filesMC)

    if not os.path.isfile(fData) or options.clean:
        filesData = " ".join(inputData)
        os.system('hadd -ff '+fData+' '+filesData)

    if options.reweight:
        
        for r in ['jetPtMax','jetHT']:
            
            frw = ROOT.TFile.Open('data/reweight/'+r+'.root','RECREATE')
            hname = 'h_'+r
        
            fmc = ROOT.TFile.Open(fMC,'READ')        
            if (fmc.GetListOfKeys().Contains(hname)):
                frw.cd()
                h_mc = fmc.Get(hname).Clone('h_mc')
                fmc.Close()
                
            fdata = ROOT.TFile.Open(fData,'READ')
            if (fdata.GetListOfKeys().Contains(hname)):
                frw.cd()
                h_data = fdata.Get(hname).Clone('h_data')
                fdata.Close()        
                
            frw.Write()
            frw.Close()
        
    picdir = 'pics'
    if os.path.isdir(picdir):
        os.system("rm -rf "+picdir)
    os.system("mkdir "+picdir)

    isQCD = '--qcd' if options.qcd else ''
    doFit = '--fit' if options.fit else ''
    
    if options.pv:
        print "Run PV study"
        os.system('python pvStudy.py --data='+fData+' --mc='+fMC+' '+isQCD)

    if options.ippv:
        print "Run IPPV study"
        os.system('python ipStudy.py --data='+fData+' --mc='+fMC+' --type=pv'+' '+isQCD+' '+doFit+' --method='+options.method+' --selection='+options.selection)

    if options.ipbs:
        print "Run IPBS study"
        os.system('python ipStudy.py --data='+fData+' --mc='+fMC+' --type=bs'+' '+isQCD+' '+doFit+' --method='+options.method+' --selection='+options.selection)
        
