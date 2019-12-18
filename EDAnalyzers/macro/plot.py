import os
import sys
import subprocess
import common as c
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
    parser.add_option("-i","--input",default="list.json",help="input file list [default: %default]")
    parser.add_option("-o","--output",default="output.root",help="output file name [default: %default]")

    (options, args) = parser.parse_args(sys.argv[1:])

    return options

if __name__ == '__main__':

    options = main()

    ROOT.gROOT.SetBatch()

    trIP_data = ROOT.TChain('residuals/treeIP')
    trIP_mc = ROOT.TChain('residuals/treeIP')

    trPV_data = ROOT.TChain('residuals/treePV')
    trPV_mc = ROOT.TChain('residuals/treePV')

    with open(options.input, "r") as read_file:
        flist = json.load(read_file)
    
    for f in flist['ZeroBias']:
        trIP_data.Add(f)
        trPV_data.Add(f)

    for f in flist['RelValMinBias13']:
        trIP_mc.Add(f)
        trPV_mc.Add(f)
        
    nDataIP = trIP_data.GetEntries()
    nDataPV = trPV_data.GetEntries()
    print 'Processed events in data: IP(' + str(nDataIP) + ') PV(' + str(nDataPV) + ')'

    nMCIP = trIP_mc.GetEntries()
    nMCPV = trPV_mc.GetEntries()    
    print 'Processed events in MC: IP(' + str(nMCIP) + ') PV(' + str(nMCPV) + ')'
    
    h = {}
    hd = {}
    h2d = {}

    hd['pvNTrks'] = {'xtit':'Number of tracks','nb':30,'xmin':0.,'xmax':150.,'ytit':'Events'}
    
    hd['ipPt'] = {'xtit':'Track p_{T} [GeV]','nb':100,'xmin':0.,'xmax':5.,'ytit':'Events'}
    hd['ipEta'] = {'xtit':'Track #eta','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
    hd['ipPhi'] = {'xtit':'Track #phi','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}

    for k, v in c.PVnTracks.iteritems():
            
        for kk, vv in c.IPpt.iteritems():
    
            nBins = 100
            xMin = -1500
            xMax = 1500
            if v[1] > 1: nBins = 150
            
            hd['ipd0'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            hd['ipdz'+k+kk] = {'xtit':'d_{z} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            hd['ipd0Err'+k+kk] = {'xtit':'#sigma(d_{xy}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            hd['ipdzErr'+k+kk] = {'xtit':'#sigma(d_{z}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            
            hd['ipd0NoRefit'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            hd['ipdzNoRefit'+k+kk] = {'xtit':'d_{z} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            hd['ipd0ErrNoRefit'+k+kk] = {'xtit':'#sigma(d_{xy}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            hd['ipdzErrNoRefit'+k+kk] = {'xtit':'#sigma(d_{z}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}

    for k, v in c.PVnTracks.iteritems():
            
        for kk, vv in c.IPeta.iteritems():
        
            hd['ipd0'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipdz'+k+kk] = {'xtit':'d_{z} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipd0Err'+k+kk] = {'xtit':'#sigma(d_{xy}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            hd['ipdzErr'+k+kk] = {'xtit':'#sigma(d_{z}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            
            hd['ipd0NoRefit'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipdzNoRefit'+k+kk] = {'xtit':'d_{z} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipd0ErrNoRefit'+k+kk] = {'xtit':'#sigma(d_{xy}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            hd['ipdzErrNoRefit'+k+kk] = {'xtit':'#sigma(d_{z}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
        
    for k, v in c.PVnTracks.iteritems():
        
        hd['pvx'+k+'__data'] = {'xtit':'x [mm]','nb':30,'xmin':0.4,'xmax':0.9,'ytit':'Events'}
        hd['pvy'+k+'__data'] = {'xtit':'y [mm]','nb':30,'xmin':0.8,'xmax':1.2,'ytit':'Events'}
        hd['pvz'+k+'__data'] = {'xtit':'z [cm]','nb':30,'xmin':-20.,'xmax':20.,'ytit':'Events'}
        
        h2d['pvx_y'+k+'__data'] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':30,'xmin':0.4,'xmax':0.9,'nby':30,'ymin':0.8,'ymax':1.2}
        h2d['pvx_z'+k+'__data'] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':30,'xmin':-20.,'xmax':20.,'nby':30,'ymin':0.4,'ymax':0.9}
        h2d['pvy_z'+k+'__data'] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':30,'xmin':-20.,'xmax':20.,'nby':30,'ymin':0.8,'ymax':1.2}

        hd['pvx'+k+'__mc'] = {'xtit':'x [mm]','nb':30,'xmin':0.8,'xmax':1.3,'ytit':'Events'}
        hd['pvy'+k+'__mc'] = {'xtit':'y [mm]','nb':30,'xmin':1.4,'xmax':2.,'ytit':'Events'}
        hd['pvz'+k+'__mc'] = {'xtit':'z [cm]','nb':30,'xmin':-20.,'xmax':20.,'ytit':'Events'}

        h2d['pvx_y'+k+'__mc'] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':30,'xmin':0.8,'xmax':1.3,'nby':30,'ymin':1.4,'ymax':2.}
        h2d['pvx_z'+k+'__mc'] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':30,'xmin':-20.,'xmax':20.,'nby':30,'ymin':0.8,'ymax':1.3}
        h2d['pvy_z'+k+'__mc'] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':30,'xmin':-20.,'xmax':20.,'nby':30,'ymin':1.4,'ymax':2.}
        
        hd['pvdx12'+k] = {'xtit':'Primary vertex resolution in x [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}
        hd['pvdy12'+k] = {'xtit':'Primary vertex resolution in y [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}
        hd['pvdz12'+k] = {'xtit':'Primary vertex resolution in z [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}

        hd['pvdxPull12'+k] = {'xtit':'Primary vertex pull in x','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}
        hd['pvdyPull12'+k] = {'xtit':'Primary vertex pull in y','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}
        hd['pvdzPull12'+k] = {'xtit':'Primary vertex pull in z','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}

    outFile = ROOT.TFile.Open(options.output,"RECREATE")
        
    h = {}
    for t in ['data','mc']:
        for k, v in hd.iteritems():
            if any(map(lambda w: w in k, ('data','mc'))):
                hname = 'h_'+k
            else:
                hname = 'h_'+k+'__'+t
            if not ROOT.gDirectory.FindObject(hname):
                h[hname] = ROOT.TH1F(hname,hname,v['nb'],v['xmin'],v['xmax'])
                h[hname].GetXaxis().SetTitle(v['xtit'])
                h[hname].GetYaxis().SetTitle(v['ytit'])
                h[hname].Sumw2()

    h2 = {}
    for t in ['data','mc']:
        for k, v in h2d.iteritems():
            if any(map(lambda w: w in k, ('data','mc'))):
                hname = 'h2_'+k
            else:
                hname = 'h2_'+k+'__'+t
            if not ROOT.gDirectory.FindObject(hname):
                h2[hname] = ROOT.TH2F(hname,hname,v['nbx'],v['xmin'],v['xmax'],v['nby'],v['ymin'],v['ymax'])
                h2[hname].GetXaxis().SetTitle(v['xtit'])
                h2[hname].GetYaxis().SetTitle(v['ytit'])
                h2[hname].Sumw2()

    for t in ['data','mc']:

        # IP study
        nStat = nDataIP
        if t == 'mc': nStat = nMCIP
        for i in range(nStat):

            trName = 'trIP_'+t
            eval(trName+'.GetEntry('+str(i)+')')

            isValid = eval(trName+'.pvIsValid')
            isFake = eval(trName+'.pvIsFake')
            
            if not isValid or isFake: continue

            nTracks = eval(trName+'.pvNTracks')
            
            hasPXL = eval(trName+'.hasPXL')
            quality = eval(trName+'.quality')
            
            pt = eval(trName+'.pt')
            eta = eval(trName+'.eta')
            phi = eval(trName+'.phi')

            d0 = eval(trName+'.d0')
            dz = eval(trName+'.dz')
            d0Err = eval(trName+'.d0Err')
            dzErr = eval(trName+'.dzErr')
            
            d0NoRefit = eval(trName+'.d0NoRefit')
            dzNoRefit = eval(trName+'.dzNoRefit')
            d0ErrNoRefit = eval(trName+'.d0ErrNoRefit')
            dzErrNoRefit = eval(trName+'.dzErrNoRefit')

            pt = eval(trName+'.pt')
            eta = eval(trName+'.eta')
            phi = eval(trName+'.phi')
            
            if math.fabs(d0) > 2000 or math.fabs(dz) > 2000: continue

            h['h_ipPt__'+t].Fill(pt)
            h['h_ipEta__'+t].Fill(eta)
            h['h_ipPhi__'+t].Fill(phi)

            for ktrk, vtrk in c.PVnTracks.iteritems():
                
                nTrkMin = vtrk[1]
                nTrkMax = vtrk[2]
                
                if not (nTracks >= nTrkMin and nTracks < nTrkMax): continue
            
                for keta, veta in c.IPeta.iteritems():
                        
                    etaMin = veta[1]
                    etaMax = veta[2]
                
                    if not (eta >= etaMin and eta < etaMax): continue

                    h['h_ipd0'+ktrk+keta+'__'+t].Fill(d0)
                    h['h_ipdz'+ktrk+keta+'__'+t].Fill(dz)
                    h['h_ipd0Err'+ktrk+keta+'__'+t].Fill(d0Err)
                    h['h_ipdzErr'+ktrk+keta+'__'+t].Fill(dzErr)
                    
                    h['h_ipd0NoRefit'+ktrk+keta+'__'+t].Fill(d0NoRefit)
                    h['h_ipdzNoRefit'+ktrk+keta+'__'+t].Fill(dzNoRefit)
                    h['h_ipd0ErrNoRefit'+ktrk+keta+'__'+t].Fill(d0ErrNoRefit)
                    h['h_ipdzErrNoRefit'+ktrk+keta+'__'+t].Fill(dzErrNoRefit)
                
                for kpt, vpt in c.IPpt.iteritems():
                    
                    ptMin = vpt[1]
                    ptMax = vpt[2]
                
                    if not (pt >= ptMin and pt < ptMax): continue
                    
                    h['h_ipd0'+ktrk+kpt+'__'+t].Fill(d0)
                    h['h_ipdz'+ktrk+kpt+'__'+t].Fill(dz)
                    h['h_ipd0Err'+ktrk+kpt+'__'+t].Fill(d0Err)
                    h['h_ipdzErr'+ktrk+kpt+'__'+t].Fill(dzErr)
                
                    h['h_ipd0NoRefit'+ktrk+kpt+'__'+t].Fill(d0NoRefit)
                    h['h_ipdzNoRefit'+ktrk+kpt+'__'+t].Fill(dzNoRefit)
                    h['h_ipd0ErrNoRefit'+ktrk+kpt+'__'+t].Fill(d0ErrNoRefit)
                    h['h_ipdzErrNoRefit'+ktrk+kpt+'__'+t].Fill(dzErrNoRefit)

        # PV study
        nStat = nDataPV
        if t == 'mc': nStat = nMCPV
        for i in range(nStat):

            trName = 'trPV_'+t
            eval(trName+'.GetEntry('+str(i)+')')
            
            isValid = eval(trName+'.pvIsValid')
            isFake = eval(trName+'.pvIsFake')
            
            if not isValid or isFake: continue

            x = eval(trName+'.pvx')
            y = eval(trName+'.pvy')
            z = eval(trName+'.pvz')
            xError = eval(trName+'.pvxError')
            yError = eval(trName+'.pvyError')
            zError = eval(trName+'.pvzError')
            nTracks = eval(trName+'.pvNTracks')
            
            x1 = eval(trName+'.pv1x')
            y1 = eval(trName+'.pv1y')
            z1 = eval(trName+'.pv1z')
            xError1 = eval(trName+'.pv1xError')
            yError1 = eval(trName+'.pv1yError')
            zError1 = eval(trName+'.pv1zError')
            
            x2 = eval(trName+'.pv2x')
            y2 = eval(trName+'.pv2y')
            z2 = eval(trName+'.pv2z')
            xError2 = eval(trName+'.pv2xError')
            yError2 = eval(trName+'.pv2yError')
            zError2 = eval(trName+'.pv2zError')
            
            dx12 = (x1-x2)/math.sqrt(2)
            dy12 = (y1-y2)/math.sqrt(2)
            dz12 = (z1-z2)/math.sqrt(2)
            
            dxPull12 = (x1-x2)/math.sqrt(xError1*xError1+xError2*xError2)
            dyPull12 = (y1-y2)/math.sqrt(yError1*yError1+yError2*yError2)
            dzPull12 = (z1-z2)/math.sqrt(zError1*zError1+zError2*zError2)

            cm = 1/10000.
            mm = 1/1000.

            h['h_pvNTrks__'+t].Fill(nTracks)
            
            for k, v in c.PVnTracks.iteritems():
                
                nTrkMin = v[1]
                nTrkMax = v[2]
                
                if nTracks >= nTrkMin and nTracks < nTrkMax:
                    
                    h['h_pvx'+k+'__'+t].Fill(x*mm)
                    h['h_pvy'+k+'__'+t].Fill(y*mm)
                    h['h_pvz'+k+'__'+t].Fill(z*cm)
                    
                    h2['h2_pvx_y'+k+'__'+t].Fill(x*mm,y*mm)
                    h2['h2_pvx_z'+k+'__'+t].Fill(z*cm,x*mm)
                    h2['h2_pvy_z'+k+'__'+t].Fill(z*cm,y*mm)
                    
                    h['h_pvdx12'+k+'__'+t].Fill(dx12)
                    h['h_pvdy12'+k+'__'+t].Fill(dy12)
                    h['h_pvdz12'+k+'__'+t].Fill(dz12)
                
                    h['h_pvdxPull12'+k+'__'+t].Fill(dxPull12)
                    h['h_pvdyPull12'+k+'__'+t].Fill(dyPull12)
                    h['h_pvdzPull12'+k+'__'+t].Fill(dzPull12)
                
    print '\033[1;32mdone\033[1;m'

    outFile.Write()
    outFile.Close()
