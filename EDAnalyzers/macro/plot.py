import os
import sys
import subprocess
import common as c
import functions as fun
import utils
import tree
import datetime as dt
from subprocess import call
import xml.etree.ElementTree as ET
from array import array
import math
import json
import ROOT

ROOT.PyConfig.IgnoreCommandLineOptions = True
from optparse import OptionParser

def main(argv = None):

    if argv == None:
        argv = sys.argv[1:]

    usage = "usage: %prog [options]\n Analysis script to create histograms"

    parser = OptionParser(usage)
    parser.add_option("--input", default="list.json", help="input file list [default: %default]")
    parser.add_option("--output", default="output.root", help="output file name [default: %default]")
    parser.add_option("--param", default="nTracks", help="list of parameterisations for PV resolution measurement [default: %default]")
    parser.add_option("--nmax", type=int, default=-1, help="number of events per job [default: %default]")
    parser.add_option("--etamax", type=float, default=3.0, help="max track eta [default: %default]")
    parser.add_option("--ptmin", type=float, default=0.1, help="min track pT in GeV [default: %default]")
    parser.add_option("--pileup", default="/user/kskovpen/analysis/Track/CMSSW_10_5_0_pre2/src/TrackingAnalysis/EDAnalyzers/macro/data/pileup/", help="path to pileup data files [default: %default]")
#    parser.add_option("--pileup", default="", help="path to pileup data files [default: %default]")
    parser.add_option("--reweight", default="/user/kskovpen/analysis/Track/CMSSW_10_5_0_pre2/src/TrackingAnalysis/EDAnalyzers/macro/data/reweight/", help="path to reweight data files [default: %default]")
#    parser.add_option("--reweight", default="", help="path to reweight data files [default: %default]")
    parser.add_option("--reweightvar", default="jetHT", help="variable to reweight [default: %default]")
    parser.add_option("--time", action='store_true', help="Print out run time information [default: %default]")

    (options, args) = parser.parse_args(sys.argv[1:])

    return options

if __name__ == '__main__':

    options = main()
    
    storeHist = True
    storeTree = False

    ROOT.gROOT.SetBatch()

    ts = {}
    
    if options.time: ts['init'] = dt.datetime.now()
    
    tr = ROOT.TChain('residuals/tree')

    with open(options.input, "r") as read_file:
        flist = json.load(read_file)
        
    isDataZeroBias = False
    isDataJetHT = False
    isMCZeroBias = False
    isMCJetHT = False
    for k,v in flist.items():
        if 'ZeroBias' in k: isDataZeroBias = True
        elif 'JetHT' in k: isDataJetHT = True
        elif 'SingleNeutrino' in k: isMCZeroBias = True
        elif 'QCD' in k: isMCJetHT = True
        for f in flist[k]: tr.Add(f)
    
    isData = (isDataZeroBias or isDataJetHT)
    isMC = (isMCZeroBias or isMCJetHT)
    
    evt = 'zb' if (isDataZeroBias or isMCZeroBias) else 'qcd'

    ppath = '/user/kskovpen/analysis/Track/CMSSW_10_5_0_pre2/src/TrackingAnalysis/EDAnalyzers/macro/data/bins/'
    param = {}
    param['bs'] = fun.param(ppath+'zb_bs.json') if evt == 'zb' else fun.param(ppath+'qcd_bs.json')
    param['bsw'] = fun.param(ppath+'zb_bsw.json') if evt == 'zb' else fun.param(ppath+'qcd_bsw.json')
    param['pv'] = fun.param(ppath+'zb_pv.json') if evt == 'zb' else fun.param(ppath+'qcd_pv.json')
    
    pvParamList = options.param.split(',')
    ipParamList = ['pt', 'eta', 'phi', 'npv', 'dr']
    
    ParamList = {}
    
    for t in ['bs', 'pv']:
        ParamList[t] = {}        
        for p in pvParamList+ipParamList:
            ParamList[t][p] = param[t].get(p)
            del ParamList[t][p]['allbins']

    for t in ['bsw']:
        ParamList[t] = {}
        for p in ['runstart', 'runend', 'lumistart', 'lumiend', 'beamwidthx', 'beamwidthy']:
            ParamList[t][p] = param[t].get(p)
    
    if options.pileup != '' and isMC:
        pu = fun.pileup(options.pileup, evt)

    if options.reweight != '' and isMC:
        rw = fun.reweight(options.reweight, 'jetHT')
        
    nEvents = tr.GetEntries()
    print 'Run on ' + ('Data' if isData else 'MC')
    print 'Processed events = ' + str(nEvents)

    sel = ['']
    for ks, vs in c.sel.iteritems():
        sel += vs.keys()
    
    if storeHist:
        
        h = {}    
        hd = {}
        h2d = {}

        hd['evNpv'] = {'xtit':'Number of primary vertices','nb':70,'xmin':0.,'xmax':70.,'ytit':'Events'}

        hs = c.hJetPtMax
        hd['jetPtMax'] = {'xtit':'Max jet p_{T} [GeV]','nb':hs[evt][0],'xmin':hs[evt][1],'xmax':hs[evt][2],'ytit':'Events'}
        hs = c.hJetHT
        hd['jetHT'] = {'xtit':'H_{T} [GeV]','nb':hs[evt][0],'xmin':hs[evt][1],'xmax':hs[evt][2],'ytit':'Events'}
        
        hs = c.hPVnTrks
        hd['pvNTrks'] = {'xtit':'Number of tracks','nb':hs[evt][0],'xmin':hs[evt][1],'xmax':hs[evt][2],'ytit':'Events'}
        hs = c.hPVsumTrackPt
        hd['pvSumTrackPt'] = {'xtit':'#sum p_{T} [GeV]','nb':hs[evt][0],'xmin':hs[evt][1],'xmax':hs[evt][2],'ytit':'Events'}
        hs = c.hPVsumTrackPt2
        hd['pvSumTrackPt2'] = {'xtit':'#sqrt{#sum p_{T}^{2}} [GeV]','nb':hs[evt][0],'xmin':hs[evt][1],'xmax':hs[evt][2],'ytit':'Events'}
        
        hd['pvXError'] = {'xtit':'Primary vertex fit uncertainty in x [#mum]','nb':50,'xmin':0.,'xmax':20.,'ytit':'Events'}
        hd['pvYError'] = {'xtit':'Primary vertex fit uncertainty in y [#mum]','nb':50,'xmin':0.,'xmax':20.,'ytit':'Events'}
        hd['pvZError'] = {'xtit':'Primary vertex fit uncertainty in z [#mum]','nb':50,'xmin':0.,'xmax':50.,'ytit':'Events'}
        
        hd['pv1XError'] = {'xtit':'Primary vertex (v1) fit uncertainty in x [#mum]','nb':50,'xmin':0.,'xmax':20.,'ytit':'Events'}
        hd['pv1YError'] = {'xtit':'Primary vertex (v1) fit uncertainty in y [#mum]','nb':50,'xmin':0.,'xmax':20.,'ytit':'Events'}
        hd['pv1ZError'] = {'xtit':'Primary vertex (v1) fit uncertainty in z [#mum]','nb':50,'xmin':0.,'xmax':50.,'ytit':'Events'}
    
        hd['pv2XError'] = {'xtit':'Primary vertex (v2) fit uncertainty in x [#mum]','nb':50,'xmin':0.,'xmax':20.,'ytit':'Events'}
        hd['pv2YError'] = {'xtit':'Primary vertex (v2) fit uncertainty in y [#mum]','nb':50,'xmin':0.,'xmax':20.,'ytit':'Events'}
        hd['pv2ZError'] = {'xtit':'Primary vertex (v2) fit uncertainty in z [#mum]','nb':50,'xmin':0.,'xmax':50.,'ytit':'Events'}
    
        hd['ipPt'] = {'xtit':'Track p_{T} [GeV]','nb':100,'xmin':0.,'xmax':20.,'ytit':'Events'}
        hd['ipEta'] = {'xtit':'Track #eta','nb':100,'xmin':-options.etamax,'xmax':options.etamax,'ytit':'Events'}
        hd['ipPhi'] = {'xtit':'Track #phi','nb':100,'xmin':-3.15,'xmax':3.15,'ytit':'Events'}
        hd['ipDrTrkJet'] = {'xtit':'#DeltaR(track,jet axis)','nb':100,'xmin':0.,'xmax':0.2,'ytit':'Events'}
        hd['ipNTrkJet'] = {'xtit':'Number of tracks in a jet','nb':15,'xmin':0.,'xmax':15,'ytit':'Events'}
        hd['ipD0'] = {'xtit':'Track d_{xy} [mm]','nb':100,'xmin':-1.5,'xmax':1.5,'ytit':'Events'}
        hd['ipDz'] = {'xtit':'Track d_{z} [mm]','nb':100,'xmin':-300.,'xmax':300.,'ytit':'Events'}
        hd['ipSD0'] = {'xtit':'Track d_{xy} significance','nb':100,'xmin':-100.,'xmax':100.,'ytit':'Events'}
        hd['ipSDz'] = {'xtit':'Track d_{z} significance','nb':100,'xmin':-300.,'xmax':300.,'ytit':'Events'}
        hd['ippvD0'] = {'xtit':'Track d_{xy}(PV) [mm]','nb':100,'xmin':-1.5,'xmax':1.5,'ytit':'Events'}
        hd['ippvDz'] = {'xtit':'Track d_{z}(PV) [mm]','nb':100,'xmin':-4.,'xmax':4.,'ytit':'Events'}
        hd['ippvSD0'] = {'xtit':'Track d_{xy}(PV) significance','nb':100,'xmin':-6.,'xmax':6.,'ytit':'Events'}
        hd['ippvSDz'] = {'xtit':'Track d_{z}(PV) significance','nb':100,'xmin':-6.,'xmax':6.,'ytit':'Events'}
        hd['ipbsD0'] = {'xtit':'Track d_{xy}(BS) [mm]','nb':100,'xmin':-1.5,'xmax':1.5,'ytit':'Events'}
        hd['ipbsDz'] = {'xtit':'Track d_{z}(BS) [mm]','nb':100,'xmin':-300.,'xmax':300.,'ytit':'Events'}
        hd['ipbsSD0'] = {'xtit':'Track d_{xy}(BS) significance','nb':100,'xmin':-6.,'xmax':6.,'ytit':'Events'}
        hd['ipbsSDz'] = {'xtit':'Track d_{z}(BS) significance','nb':100,'xmin':-300.,'xmax':300.,'ytit':'Events'}
        hd['ipbszpcaD0'] = {'xtit':'Track d_{xy}(BS,z=PCA) [mm]','nb':100,'xmin':-1.5,'xmax':1.5,'ytit':'Events'}
        hd['ipbszpcaSD0'] = {'xtit':'Track d_{xy}(BS,z=PCA) significance [mm]','nb':100,'xmin':-6.,'xmax':6.,'ytit':'Events'}
        hd['ipbszpvD0'] = {'xtit':'Track d_{xy}(BS,z=PV) [mm]','nb':100,'xmin':-1.5,'xmax':1.5,'ytit':'Events'}
        hd['ipbszpvSD0'] = {'xtit':'Track d_{xy}(BS,z=PV) significance [mm]','nb':100,'xmin':-6.,'xmax':6.,'ytit':'Events'}
        
        for ksel in sel:
        
            for p in ipParamList:
            
                bins = ParamList['pv'][p]
            
                for kk, vv in bins.iteritems():
                    
                    if kk in ['']: continue
                    
                    IP = bins[kk]
                    
                    for pvp in pvParamList:
                        
                        pvbins = ParamList['pv'][pvp]
                        
                        for kpv, vpv in pvbins.iteritems():
                            
                            hd['ippvd0'+kpv+kk+ksel] = {'xtit':'d_{xy}(PV) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
                            hd['ippvdz'+kpv+kk+ksel] = {'xtit':'d_{z}(PV) [#mum]','nb':IP['dz'][0],'xmin':IP['dz'][1],'xmax':IP['dz'][2],'ytit':'Events'}

                bins = ParamList['bs'][p]
            
                for kk, vv in bins.iteritems():

                    if kk in ['']: continue
                
                    IP = bins[kk]
                
                    if isData:
                        bsbins = ParamList['bsw']['beamwidthx']
                
                        for ibs in range(len(bsbins)):
                        
                            hd['ipbsd0'+kk+'_'+str(ibs)+ksel] = {'xtit':'d_{xy}(BS) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
                            hd['ipbsd0zpv'+kk+'_'+str(ibs)+ksel] = {'xtit':'d_{xy}(BS) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
                            hd['ipbsdz'+kk+'_'+str(ibs)+ksel] = {'xtit':'d_{z}(BS) [#mum]','nb':IP['dz'][0],'xmin':IP['dz'][1],'xmax':IP['dz'][2],'ytit':'Events'}
                    else:
                        hd['ipbsd0'+kk+ksel] = {'xtit':'d_{xy}(BS) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
                        hd['ipbsd0zpv'+kk+ksel] = {'xtit':'d_{xy}(BS) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
                        hd['ipbsdz'+kk+ksel] = {'xtit':'d_{z}(BS) [#mum]','nb':IP['dz'][0],'xmin':IP['dz'][1],'xmax':IP['dz'][2],'ytit':'Events'}

        if isData:
            PV = c.hPVData
            hd['pvx'] = {'xtit':'x [mm]','nb':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'ytit':'Events'}
            hd['pvy'] = {'xtit':'y [mm]','nb':PV['y'][0],'xmin':PV['y'][1],'xmax':PV['y'][2],'ytit':'Events'}
            hd['pvz'] = {'xtit':'z [cm]','nb':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'ytit':'Events'}
            
            h2d['pvx_y'] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            h2d['pvx_z'] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['x'][0],'ymin':PV['x'][1],'ymax':PV['x'][2]}
            h2d['pvy_z'] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            
        else:
            PV = c.hPVMC
            hd['pvx'] = {'xtit':'x [mm]','nb':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'ytit':'Events'}
            hd['pvy'] = {'xtit':'y [mm]','nb':PV['y'][0],'xmin':PV['y'][1],'xmax':PV['y'][2],'ytit':'Events'}
            hd['pvz'] = {'xtit':'z [cm]','nb':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'ytit':'Events'}
            
            h2d['pvx_y'] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            h2d['pvx_z'] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['x'][0],'ymin':PV['x'][1],'ymax':PV['x'][2]}
            h2d['pvy_z'] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            
        for pvp in pvParamList:
            
            pvbins = ParamList['pv'][pvp]
            
            for k, v in pvbins.iteritems():
            
                hd['pvdx12'+k] = {'xtit':'Primary vertex resolution in x [#mum]','nb':v['resox'][0],'xmin':v['resox'][1],'xmax':v['resox'][2],'ytit':'Events'}
                hd['pvdy12'+k] = {'xtit':'Primary vertex resolution in y [#mum]','nb':v['resoy'][0],'xmin':v['resoy'][1],'xmax':v['resoy'][2],'ytit':'Events'}
                hd['pvdz12'+k] = {'xtit':'Primary vertex resolution in z [#mum]','nb':v['resoz'][0],'xmin':v['resoz'][1],'xmax':v['resoz'][2],'ytit':'Events'}
                
                hd['pvdxPull12'+k] = {'xtit':'Primary vertex pull in x','nb':v['pullx'][0],'xmin':v['pullx'][1],'xmax':v['pullx'][2],'ytit':'Events'}
                hd['pvdyPull12'+k] = {'xtit':'Primary vertex pull in y','nb':v['pully'][0],'xmin':v['pully'][1],'xmax':v['pully'][2],'ytit':'Events'}
                hd['pvdzPull12'+k] = {'xtit':'Primary vertex pull in z','nb':v['pullz'][0],'xmin':v['pullz'][1],'xmax':v['pullz'][2],'ytit':'Events'}

        if isData:
            PV = c.hPVData
            hd['bsx0'] = {'xtit':'x0 [mm]','nb':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'ytit':'Events'}
            hd['bsy0'] = {'xtit':'y0 [mm]','nb':PV['y'][0],'xmin':PV['y'][1],'xmax':PV['y'][2],'ytit':'Events'}
            hd['bsz0'] = {'xtit':'z0 [cm]','nb':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'ytit':'Events'}
            
            hd['bsx'] = {'xtit':'x [mm]','nb':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'ytit':'Events'}
            hd['bsy'] = {'xtit':'y [mm]','nb':PV['y'][0],'xmin':PV['y'][1],'xmax':PV['y'][2],'ytit':'Events'}
            
            h2d['bsx0_y0'] = {'xtit':'x0 [mm]','ytit':'y0 [mm]','nbx':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            h2d['bsx0_z0'] = {'xtit':'z0 [cm]','ytit':'x0 [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['x'][0],'ymin':PV['x'][1],'ymax':PV['x'][2]}
            h2d['bsy0_z0'] = {'xtit':'z0 [cm]','ytit':'y0 [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            
            h2d['bsx_y'] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            h2d['bsx_z'] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['x'][0],'ymin':PV['x'][1],'ymax':PV['x'][2]}
            h2d['bsy_z'] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            
            BS = c.hBSData
            hd['bsBeamWidthX'] = {'xtit':'Beam width (x) [#mum]','nb':BS['bwx'][0],'xmin':BS['bwx'][1],'xmax':BS['bwx'][2],'ytit':'Events'}
            hd['bsBeamWidthY'] = {'xtit':'Beam width (y) [#mum]','nb':BS['bwy'][0],'xmin':BS['bwy'][1],'xmax':BS['bwy'][2],'ytit':'Events'}
            hd['bsSigmaZ'] = {'xtit':'Beam sigma z [mm]','nb':BS['bwz'][0],'xmin':BS['bwz'][1],'xmax':BS['bwz'][2],'ytit':'Events'}
            
        else:
            PV = c.hPVMC
            hd['bsx0'] = {'xtit':'x0 [mm]','nb':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'ytit':'Events'}
            hd['bsy0'] = {'xtit':'y0 [mm]','nb':PV['y'][0],'xmin':PV['y'][1],'xmax':PV['y'][2],'ytit':'Events'}
            hd['bsz0'] = {'xtit':'z0 [cm]','nb':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'ytit':'Events'}
            
            hd['bsx'] = {'xtit':'x [mm]','nb':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'ytit':'Events'}
            hd['bsy'] = {'xtit':'y [mm]','nb':PV['y'][0],'xmin':PV['y'][1],'xmax':PV['y'][2],'ytit':'Events'}
            
            h2d['bsx0_y0'] = {'xtit':'x0 [mm]','ytit':'y0 [mm]','nbx':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            h2d['bsx0_z0'] = {'xtit':'z0 [cm]','ytit':'x0 [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['x'][0],'ymin':PV['x'][1],'ymax':PV['x'][2]}
            h2d['bsy0_z0'] = {'xtit':'z0 [cm]','ytit':'y0 [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            
            h2d['bsx_y'] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            h2d['bsx_z'] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['x'][0],'ymin':PV['x'][1],'ymax':PV['x'][2]}
            h2d['bsy_z'] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            
            BS = c.hBSMC
            hd['bsBeamWidthX'] = {'xtit':'Beam width (x) [#mum]','nb':BS['bwx'][0],'xmin':BS['bwx'][1],'xmax':BS['bwx'][2],'ytit':'Events'}
            hd['bsBeamWidthY'] = {'xtit':'Beam width (y) [#mum]','nb':BS['bwy'][0],'xmin':BS['bwy'][1],'xmax':BS['bwy'][2],'ytit':'Events'}
            hd['bsSigmaZ'] = {'xtit':'Beam sigma z [mm]','nb':BS['bwz'][0],'xmin':BS['bwz'][1],'xmax':BS['bwz'][2],'ytit':'Events'}

    if options.time: ts['hist'] = dt.datetime.now()
    
    outFile = ROOT.TFile.Open(options.output,"RECREATE")
    
    if storeTree: 
        outFile.SetCompressionAlgorithm(ROOT.ROOT.kLZ4)
        outFile.SetCompressionLevel(3)
        trkTree = tree.trackTree()

    if storeHist:
        
        h = {}
        for k, v in hd.iteritems():
            hname = 'h_'+k
            if not ROOT.gDirectory.FindObject(hname):
                h[hname] = ROOT.TH1F(hname,hname,v['nb'],v['xmin'],v['xmax'])
                h[hname].GetXaxis().SetTitle(v['xtit'])
                h[hname].GetYaxis().SetTitle(v['ytit'])
                h[hname].Sumw2()

        h2 = {}
        for k, v in h2d.iteritems():
            hname = 'h2_'+k
            if not ROOT.gDirectory.FindObject(hname):
                h2[hname] = ROOT.TH2F(hname,hname,v['nbx'],v['xmin'],v['xmax'],v['nby'],v['ymin'],v['ymax'])
                h2[hname].GetXaxis().SetTitle(v['xtit'])
                h2[hname].GetYaxis().SetTitle(v['ytit'])
                h2[hname].Sumw2()

    if options.time: 
        ts['main'] = dt.datetime.now()
                
    # Fill histograms
    for i in range(nEvents):
        
        ts['main_checkpoints'] = []
        
        if i > options.nmax and options.nmax >= 0: break

        if options.time: ts['main_checkpoints'].append(dt.datetime.now())
        
        tr.GetEntry(i)
        
        if options.time: ts['main_checkpoints'].append(dt.datetime.now())

        # Event selection
        
        jetPtMax = 0.
        jetHT = 0.
        for j in range(tr.pfjet_n):
            
            jpt = tr.pfjet_pt[j]
            jeta = tr.pfjet_eta[j]
            
            if jpt > jetPtMax: jetPtMax = jpt
            
            if jpt < 50 or math.fabs(jeta) > 2.5: continue
            jetHT += tr.pfjet_pt[j]
        
        if isDataZeroBias or isMCZeroBias:                
            if not tr.trig_ZeroBias_pass: continue

        if isDataJetHT or isMCJetHT:
            if jetPtMax < 50. or jetHT < 200.: continue
            if not (tr.trig_PFHT180_pass or tr.trig_PFHT250_pass or tr.trig_PFHT370_pass or\
            tr.trig_PFHT430_pass or tr.trig_PFHT510_pass or tr.trig_PFHT590_pass or\
            tr.trig_PFHT680_pass or tr.trig_PFHT780_pass or tr.trig_PFHT890_pass or tr.trig_PFHT1050_pass): continue
            
        isValid = tr.pv_IsValid
        isFake = tr.pv_IsFake
            
        if not isValid or isFake: continue

        nTracks = tr.pv_NTracks
        sumTrackPt = tr.pv_SumTrackPt
        sumTrackPtSq = math.sqrt(tr.pv_SumTrackPt2)
        
        if nTracks < 3: continue

        pvParamListVal = []
        for pvp in pvParamList:
            pvParamListVal.append(eval(pvp))
        
        npv = tr.ev_nPV
        
        we = 1.
        if options.pileup != '' and isMC:
            we = we*pu.getWeight(tr.mc_pu_trueNumInt)

        if options.reweight != '' and isMC and evt == 'qcd':
            we = we*rw.getWeight(eval(options.reweightvar))
            
        run = tr.ev_run
        lumi = tr.ev_lumi

        if options.time: ts['main_checkpoints'].append(dt.datetime.now())
        
        # PV/BS study
        
        bs_x0 = tr.bs_x0
        bs_y0 = tr.bs_y0
        bs_z0 = tr.bs_z0
        
        bs_x = tr.bs_x_zpv
        bs_y = tr.bs_y_zpv
        
        bs_beamWidthX = tr.bs_BeamWidthX
        bs_beamWidthY = tr.bs_BeamWidthY
        bs_sigmaZ = tr.bs_sigmaZ*10 # mm
        bs_beamWidthXError = tr.bs_BeamWidthXError
        bs_beamWidthYError = tr.bs_BeamWidthYError
        bs_sigmaZError = tr.bs_sigmaZ0Error*10 # mm
        bs_emittanceX = tr.bs_emittanceX
        bs_emittanceY = tr.bs_emittanceY
        bs_betaStar = tr.bs_betaStar

        # from BS monitoring
        if (bs_beamWidthXError >= 0.1): continue

        pv_x = tr.pv_x
        pv_y = tr.pv_y
        pv_z = tr.pv_z
        pv_xError = tr.pv_xError
        pv_yError = tr.pv_yError
        pv_zError = tr.pv_zError

        pv_x1 = tr.pv_x_p1
        pv_y1 = tr.pv_y_p1
        pv_z1 = tr.pv_z_p1
        pv_xError1 = tr.pv_xError_p1
        pv_yError1 = tr.pv_yError_p1
        pv_zError1 = tr.pv_zError_p1
        
        pv_x2 = tr.pv_x_p2
        pv_y2 = tr.pv_y_p2
        pv_z2 = tr.pv_z_p2
        pv_xError2 = tr.pv_xError_p2
        pv_yError2 = tr.pv_yError_p2
        pv_zError2 = tr.pv_zError_p2
        
        if pv_x1 == -777 or pv_x2 == -777: continue
            
        pv_dx12 = (pv_x1-pv_x2)/math.sqrt(2)
        pv_dy12 = (pv_y1-pv_y2)/math.sqrt(2)
        pv_dz12 = (pv_z1-pv_z2)/math.sqrt(2)
        
        pv_dxPull12 = (pv_x1-pv_x2)/math.sqrt(pv_xError1*pv_xError1+pv_xError2*pv_xError2)
        pv_dyPull12 = (pv_y1-pv_y2)/math.sqrt(pv_yError1*pv_yError1+pv_yError2*pv_yError2)
        pv_dzPull12 = (pv_z1-pv_z2)/math.sqrt(pv_zError1*pv_zError1+pv_zError2*pv_zError2)

        pv_cm = 1/10000.
        pv_mm = 1/1000.

        bs_cm = 1.
        bs_mm = 10.
        bs_micron = 10000.
        
        if bs_beamWidthXError*bs_micron > 0.1: continue

        if storeHist:
            
            h['h_evNpv'].Fill(npv, we)            
            h['h_jetPtMax'].Fill(jetPtMax, we)
            h['h_jetHT'].Fill(jetHT, we)
            
            h['h_pvNTrks'].Fill(nTracks, we)
            h['h_pvSumTrackPt'].Fill(sumTrackPt, we)
            h['h_pvSumTrackPt2'].Fill(sumTrackPtSq, we)
            
            h['h_pvXError'].Fill(pv_xError, we)
            h['h_pvYError'].Fill(pv_yError, we)
            h['h_pvZError'].Fill(pv_zError, we)
            
            h['h_pv1XError'].Fill(pv_xError1, we)
            h['h_pv1YError'].Fill(pv_yError1, we)
            h['h_pv1ZError'].Fill(pv_zError1, we)
            
            h['h_pv2XError'].Fill(pv_xError2, we)
            h['h_pv2YError'].Fill(pv_yError2, we)
            h['h_pv2ZError'].Fill(pv_zError2, we)
            
            h['h_bsx0'].Fill(bs_x0*bs_mm, we)
            h['h_bsy0'].Fill(bs_y0*bs_mm, we)
            h['h_bsz0'].Fill(bs_z0*bs_cm, we)
            
            h['h_bsx'].Fill(bs_x*bs_mm, we)
            h['h_bsy'].Fill(bs_y*bs_mm, we)
            
            h2['h2_bsx0_y0'].Fill(bs_x0*bs_mm, bs_y0*bs_mm, we)
            h2['h2_bsx0_z0'].Fill(bs_z0*bs_cm, bs_x0*bs_mm, we)
            h2['h2_bsy0_z0'].Fill(bs_z0*bs_cm, bs_y0*bs_mm, we)
            
            h2['h2_bsx_y'].Fill(bs_x*bs_mm, bs_y*bs_mm, we)
            h2['h2_bsx_z'].Fill(bs_z0*bs_cm, bs_x*bs_mm, we)
            h2['h2_bsy_z'].Fill(bs_z0*bs_cm, bs_y*bs_mm, we)
            
            h['h_bsBeamWidthX'].Fill(bs_beamWidthX*bs_micron, we)
            h['h_bsBeamWidthY'].Fill(bs_beamWidthY*bs_micron, we)
            h['h_bsSigmaZ'].Fill(bs_sigmaZ*bs_cm, we)

            h['h_pvx'].Fill(pv_x*pv_mm, we)
            h['h_pvy'].Fill(pv_y*pv_mm, we)
            h['h_pvz'].Fill(pv_z*pv_cm, we)
            
            h2['h2_pvx_y'].Fill(pv_x*pv_mm, pv_y*pv_mm, we)
            h2['h2_pvx_z'].Fill(pv_z*pv_cm, pv_x*pv_mm, we)
            h2['h2_pvy_z'].Fill(pv_z*pv_cm, pv_y*pv_mm, we)
            
            for pvp in pvParamList:
                
                param = eval(pvp)
                pvbins = ParamList['pv'][pvp]
                
                for k, v in pvbins.iteritems():
                    
                    paramMin = v['bins'][1]
                    paramMax = v['bins'][2]
                    
                    if param >= paramMin and param < paramMax:
                        
                        h['h_pvdx12'+k].Fill(pv_dx12, we)
                        h['h_pvdy12'+k].Fill(pv_dy12, we)
                        h['h_pvdz12'+k].Fill(pv_dz12, we)
                        
                        h['h_pvdxPull12'+k].Fill(pv_dxPull12, we)
                        h['h_pvdyPull12'+k].Fill(pv_dyPull12, we)
                        h['h_pvdzPull12'+k].Fill(pv_dzPull12, we)

        if storeTree:
            
            trkTree.clear()
            
            trkTree.run[0] = run
            trkTree.lumi[0] = lumi
            trkTree.beamWidthX[0] = bs_beamWidthX
            trkTree.beamWidthY[0] = bs_beamWidthY
            trkTree.beamSigmaZ[0] = bs_sigmaZ
            trkTree.sumTrackPtSq[0] = sumTrackPtSq
            trkTree.beamWidthXError[0] = bs_beamWidthXError
            trkTree.beamWidthYError[0] = bs_beamWidthYError
            trkTree.beamSigmaZError[0] = bs_sigmaZError
            trkTree.beamEmittanceX[0] = bs_emittanceX
            trkTree.beamEmittanceY[0] = bs_emittanceY
            trkTree.beamBetaStar[0] = bs_betaStar

        if options.time: ts['main_checkpoints'].append(dt.datetime.now())
        
        # IP study
        
        # Fast access to vector branches
        
        trk_pt = tr.trk_pt
        trk_eta = tr.trk_eta
        trk_phi = tr.trk_phi
        
        trk_hasPixelBarrelLayer1 = tr.trk_hasPixelBarrelLayer1
        trk_hasPixelEndcapLayer1 = tr.trk_hasPixelEndcapLayer1
        trk_hasPixelBarrelLayer2 = tr.trk_hasPixelBarrelLayer2
        trk_hasPixelEndcapLayer2 = tr.trk_hasPixelEndcapLayer2
        trk_hasPixelBarrelLayer3 = tr.trk_hasPixelBarrelLayer3
        trk_hasPixelEndcapLayer3 = tr.trk_hasPixelEndcapLayer3
        trk_hasPixelBarrelLayer4 = tr.trk_hasPixelBarrelLayer4
        trk_hasPixelEndcapLayer4 = tr.trk_hasPixelEndcapLayer4
        
        trk_jet_found = tr.trk_jet_found
        trk_jet_eta = tr.trk_jet_eta
        trk_jet_phi = tr.trk_jet_phi
        trk_jet_nTracks = tr.trk_jet_nTracks
        
        trk_d0 = tr.trk_d0
        trk_dz = tr.trk_dz
        
        trk_d0_pv = tr.trk_d0_pv
        trk_dz_pv = tr.trk_dz_pv
        trk_d0_bs = tr.trk_d0_bs
        trk_dz_bs = tr.trk_dz_bs
        
        trk_d0_bs_zpv = tr.trk_d0_bs_zpv
        trk_d0_bs_zpca = tr.trk_d0_bs_zpca
        
        trk_d0Err = tr.trk_d0Err
        trk_dzErr = tr.trk_dzErr
        
        nTracks = trk_pt.size()

        if options.time: print 'nTracks =', nTracks
        
        for t in range(nTracks):

            pt = trk_pt[t]
            eta = trk_eta[t]
            phi = trk_phi[t]
            
            hasPXL1 = (trk_hasPixelBarrelLayer1[t] or trk_hasPixelEndcapLayer1[t])
            hasPXL2 = (trk_hasPixelBarrelLayer2[t] or trk_hasPixelEndcapLayer2[t])
            hasPXL3 = (trk_hasPixelBarrelLayer3[t] or trk_hasPixelEndcapLayer3[t])
            hasPXL4 = (trk_hasPixelBarrelLayer4[t] or trk_hasPixelEndcapLayer4[t])
            
            trkSelPXL1 = bool(hasPXL1)
            trkSelPXL2 = bool(hasPXL1 and hasPXL2)
            trkSelPXL3 = bool(hasPXL1 and hasPXL2 and hasPXL3)
            trkSelPXL4 = bool(hasPXL1 and hasPXL2 and hasPXL3 and hasPXL4)
            
##            if pt > 10.: continue
##            if pt < 1.: continue

            isJet = trk_jet_found[t]
            jetEta = trk_jet_eta[t]
            jetPhi = trk_jet_phi[t]
            nTrkJet = trk_jet_nTracks[t]

            if isJet:
                drTrkJet = utils.deltaR2(eta, phi, jetEta, jetPhi)
            
            d0 = trk_d0[t]
            dz = trk_dz[t]

            d0_pv = trk_d0_pv[t]
            dz_pv = trk_dz_pv[t]
            d0_bs = trk_d0_bs[t]
            dz_bs = trk_dz_bs[t]
            d0_bs_zpv = trk_d0_bs_zpv[t]
            d0_bs_zpca = trk_d0_bs_zpca[t]
            
            d0Err = trk_d0Err[t]
            dzErr = trk_dzErr[t]

            if pt < options.ptmin: continue
            if math.fabs(eta) > options.etamax: continue

            if storeHist:
                
                h['h_ipPt'].Fill(pt, we)
                h['h_ipEta'].Fill(eta, we)
                h['h_ipPhi'].Fill(phi, we)
                if isJet: h['h_ipDrTrkJet'].Fill(drTrkJet, we)
                if isJet: h['h_ipNTrkJet'].Fill(nTrkJet, we)

            sd0 = d0/d0Err if d0Err > 0 else -777
            sdz = dz/dzErr if dzErr > 0 else -777
            sd0_pv = d0_pv/d0Err if d0Err > 0 else -777
            sdz_pv = dz_pv/dzErr if dzErr > 0 else -777
            sd0_bs = d0_bs/d0Err if d0Err > 0 else -777
            sdz_bs = dz_bs/dzErr if dzErr > 0 else -777
            sd0_bs_zpv = d0_bs_zpv/d0Err if d0Err > 0 else -777
            sd0_bs_zpca = d0_bs_zpca/d0Err if d0Err > 0 else -777
            
            mm = 1E-3

            if storeTree:
            
                trkTree.pt.push_back(pt)
                trkTree.eta.push_back(eta)
                trkTree.phi.push_back(phi)
                trkTree.npv.push_back(npv)
                if isJet: trkTree.dr.push_back(drTrkJet)
            
            if storeHist:
                
                h['h_ipD0'].Fill(d0*mm, we)
                h['h_ipDz'].Fill(dz*mm, we)
                h['h_ipSD0'].Fill(sd0, we)
                h['h_ipSDz'].Fill(sdz, we)
                
                h['h_ippvD0'].Fill(d0_pv*mm, we)
                h['h_ippvDz'].Fill(dz_pv*mm, we)
                h['h_ippvSD0'].Fill(sd0_pv, we)
                h['h_ippvSDz'].Fill(sdz_pv, we)
                
                h['h_ipbsD0'].Fill(d0_bs*mm, we)
                h['h_ipbsDz'].Fill(dz_bs*mm, we)
                h['h_ipbsSD0'].Fill(sd0_bs, we)
                h['h_ipbsSDz'].Fill(sdz_bs, we)
                
                h['h_ipbszpvD0'].Fill(d0_bs_zpv*mm, we)
                h['h_ipbszpcaD0'].Fill(d0_bs_zpca*mm, we)
                h['h_ipbszpvSD0'].Fill(sd0_bs_zpv, we)
                h['h_ipbszpcaSD0'].Fill(sd0_bs_zpca, we)

                sellist = ['']
                
                for ksel, vsel in c.sel.iteritems():
                    
                    if ksel == 'pt': sv = pt
                    elif ksel == 'eta': sv = math.fabs(eta)
                    else:
                        print 'Uknown selection:', ksel
                        sys.exit()
                    
                    for kksel, vvsel in vsel.iteritems():
                        
                        vselMin = vvsel[0]
                        vselMax = vvsel[1]
                        if (sv >= vselMin and sv < vselMax):
                            sellist.append(kksel)
                
                for p in ipParamList:

                    varp = pt
                    if p == 'eta': varp = eta
                    elif p == 'phi': varp = phi
                    elif p == 'npv': varp = npv
                    elif p == 'dr':
                        if not isJet: continue
                        else: varp = drTrkJet
                        
                    blistpv = ParamList['pv'][p]
                    
                    for kp, vp in blistpv.iteritems():

                        if kp in ['']: continue
                        
                        IP = blistpv[kp]['bins']

                        pMin = IP[1]
                        pMax = IP[2]
                        
                        if not (varp >= pMin and varp < pMax): continue

                        for ipvp, pvp in enumerate(pvParamList):
 
                            param = pvParamListVal[ipvp]
                            pvbins = ParamList['pv'][pvp]
                            
                            for kpv, vpv in pvbins.iteritems():

                                if kpv in ['']: continue
                                
                                paramEdge = vpv['bins']
                                
                                paramMin = paramEdge[1]
                                paramMax = paramEdge[2]
                
                                if not (param >= paramMin and param < paramMax): continue
                        
                                for kselip in sellist:

#                                    if kpv+kp+kselip == '_sumTrackPtSq10p587to10p868_pt0p8207831493914127to0p931683250322938':
#                                        print d0_pv, dz_pv
                                    
                                    h['h_ippvd0'+kpv+kp+kselip].Fill(d0_pv, we)
                                    h['h_ippvdz'+kpv+kp+kselip].Fill(dz_pv, we)

                    blistbs = ParamList['bs'][p]
                                
                    for kp, vp in blistbs.iteritems():

                        if kp in ['']: continue
                        
                        IP = blistbs[kp]['bins']

                        pMin = IP[1]
                        pMax = IP[2]
                        
                        if not (varp >= pMin and varp < pMax): continue

                        if isData:
                            
                            bsbins = ParamList['bsw']
                            
                            for ibs in range(len(bsbins['runstart'])):
                                
                                runMin = bsbins['runstart'][ibs]
                                runMax = bsbins['runend'][ibs]
                                lumiMin = bsbins['lumistart'][ibs]
                                lumiMax = bsbins['lumiend'][ibs]
                                
                                if not (run >= runMin and run < runMax): continue
                                if (run == runMin) and (lumi < lumiMin): continue
                                if (run == runMax) and (lumi >= lumiMax): continue

                                for kselip in sellist:
                                
                                    h['h_ipbsd0zpv'+kp+'_'+str(ibs)+kselip].Fill(d0_bs_zpv, we)
                                    h['h_ipbsd0'+kp+'_'+str(ibs)+kselip].Fill(d0_bs, we)
                                    h['h_ipbsdz'+kp+'_'+str(ibs)+kselip].Fill(dz_bs, we)
                        else:
                            
                            for kselip in sellist:
                                
                                h['h_ipbsd0zpv'+kp+kselip].Fill(d0_bs_zpv, we)
                                h['h_ipbsd0'+kp+kselip].Fill(d0_bs, we)
                                h['h_ipbsdz'+kp+kselip].Fill(dz_bs, we)
                                
        if storeTree: trkTree.fill()
        
        if options.time: ts['main_checkpoints'].append(dt.datetime.now())

        if options.time: print '     -> ', \
        'total: ', (ts['main_checkpoints'][-1]-ts['main_checkpoints'][0]).total_seconds(), 's', \
        'read: ', (ts['main_checkpoints'][1]-ts['main_checkpoints'][0]).total_seconds(), 's', \
        'selection: ', (ts['main_checkpoints'][2]-ts['main_checkpoints'][1]).total_seconds(), 's', \
        'pv: ', (ts['main_checkpoints'][3]-ts['main_checkpoints'][2]).total_seconds(), 's', \
        'ip: ', (ts['main_checkpoints'][4]-ts['main_checkpoints'][3]).total_seconds(), 's'

    if options.time: ts['end'] = dt.datetime.now()
    
    if options.time:
        
        print '----> Runtime stats ----------->'
        print 'Initialisation     =', (ts['hist']-ts['init']).total_seconds(), 's'
        print 'Histogram booking  =', (ts['main']-ts['hist']).total_seconds(), 's'
        print 'Event loop         =', (ts['end']-ts['main']).total_seconds(), 's'
        print '---->-------------------------->'
    
    print '\033[1;32mdone\033[1;m'

    outFile.Write()
    outFile.Close()
