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

    tr = ROOT.TChain('residuals/tree')

    with open(options.input, "r") as read_file:
        flist = json.load(read_file)

    isData = False
    for k,v in flist.items():
        if 'MinimumBias' in k: isData = True
        for f in flist[k]: tr.Add(f)
    
    nEvents = tr.GetEntries()
    print 'Run on ' + ('Data' if isData else 'MC')
    print 'Processed events = ' + str(nEvents)

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

    for k, v in c.PVnTracks.iteritems():
            
        for kk, vv in c.IPeta.iteritems():
        
            hd['ipd0'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipdz'+k+kk] = {'xtit':'d_{z} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipd0Err'+k+kk] = {'xtit':'#sigma(d_{xy}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            hd['ipdzErr'+k+kk] = {'xtit':'#sigma(d_{z}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            
            hd['ipd0NoRefit'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipdzNoRefit'+k+kk] = {'xtit':'d_{z} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
        
    for k, v in c.PVnTracks.iteritems():

        if isData:
            hd['pvx'+k] = {'xtit':'x [mm]','nb':30,'xmin':0.4,'xmax':0.9,'ytit':'Events'}
            hd['pvy'+k] = {'xtit':'y [mm]','nb':30,'xmin':0.8,'xmax':1.2,'ytit':'Events'}
            hd['pvz'+k] = {'xtit':'z [cm]','nb':30,'xmin':-20.,'xmax':20.,'ytit':'Events'}
            
            h2d['pvx_y'+k] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':30,'xmin':0.4,'xmax':0.9,'nby':30,'ymin':0.8,'ymax':1.2}
            h2d['pvx_z'+k] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':30,'xmin':-20.,'xmax':20.,'nby':30,'ymin':0.4,'ymax':0.9}
            h2d['pvy_z'+k] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':30,'xmin':-20.,'xmax':20.,'nby':30,'ymin':0.8,'ymax':1.2}

        else:
            hd['pvx'+k] = {'xtit':'x [mm]','nb':30,'xmin':0.8,'xmax':1.3,'ytit':'Events'}
            hd['pvy'+k] = {'xtit':'y [mm]','nb':30,'xmin':1.4,'xmax':2.,'ytit':'Events'}
            hd['pvz'+k] = {'xtit':'z [cm]','nb':30,'xmin':-20.,'xmax':20.,'ytit':'Events'}
            
            h2d['pvx_y'+k] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':30,'xmin':0.8,'xmax':1.3,'nby':30,'ymin':1.4,'ymax':2.}
            h2d['pvx_z'+k] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':30,'xmin':-20.,'xmax':20.,'nby':30,'ymin':0.8,'ymax':1.3}
            h2d['pvy_z'+k] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':30,'xmin':-20.,'xmax':20.,'nby':30,'ymin':1.4,'ymax':2.}
        
        hd['pvdx12'+k] = {'xtit':'Primary vertex resolution in x [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}
        hd['pvdy12'+k] = {'xtit':'Primary vertex resolution in y [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}
        hd['pvdz12'+k] = {'xtit':'Primary vertex resolution in z [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}

        hd['pvdxPull12'+k] = {'xtit':'Primary vertex pull in x','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}
        hd['pvdyPull12'+k] = {'xtit':'Primary vertex pull in y','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}
        hd['pvdzPull12'+k] = {'xtit':'Primary vertex pull in z','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}

    outFile = ROOT.TFile.Open(options.output,"RECREATE")
        
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

    # IP study
    for i in range(nEvents):

        tr.GetEntry(i)

        isValid = tr.pv_IsValid
        isFake = tr.pv_IsFake
            
        if not isValid or isFake: continue

        nPVTracks = tr.pv_NTracks
        
        nTracks = tr.trk_pt.size()
            
        for t in range(nTracks):
        
            hasPXL = tr.trk_hasPXL[t]
            quality = tr.trk_quality[t]
            
            pt = tr.trk_pt[t]
            eta = tr.trk_eta[t]
            phi = tr.trk_phi[t]

            d0 = tr.trk_d0_pv[t]
            dz = tr.trk_dz_pv[t]
            d0Err = tr.trk_d0Err[t]
            dzErr = tr.trk_dzErr[t]

            d0NoRefit = tr.trk_d0_pv_NoRefit[t]
            dzNoRefit = tr.trk_dz_pv_NoRefit[t]
            
            if math.fabs(d0) > 2000 or math.fabs(dz) > 2000: continue

            h['h_ipPt'].Fill(pt)
            h['h_ipEta'].Fill(eta)
            h['h_ipPhi'].Fill(phi)

            for ktrk, vtrk in c.PVnTracks.iteritems():
                
                nTrkMin = vtrk[1]
                nTrkMax = vtrk[2]
                
                if not (nPVTracks >= nTrkMin and nPVTracks < nTrkMax): continue
            
                for keta, veta in c.IPeta.iteritems():
                        
                    etaMin = veta[1]
                    etaMax = veta[2]
                
                    if not (eta >= etaMin and eta < etaMax): continue

                    h['h_ipd0'+ktrk+keta].Fill(d0)
                    h['h_ipdz'+ktrk+keta].Fill(dz)
                    h['h_ipd0Err'+ktrk+keta].Fill(d0Err)
                    h['h_ipdzErr'+ktrk+keta].Fill(dzErr)
                    
                    h['h_ipd0NoRefit'+ktrk+keta].Fill(d0NoRefit)
                    h['h_ipdzNoRefit'+ktrk+keta].Fill(dzNoRefit)
                
                for kpt, vpt in c.IPpt.iteritems():
                    
                    ptMin = vpt[1]
                    ptMax = vpt[2]
                
                    if not (pt >= ptMin and pt < ptMax): continue
                    
                    h['h_ipd0'+ktrk+kpt].Fill(d0)
                    h['h_ipdz'+ktrk+kpt].Fill(dz)
                    h['h_ipd0Err'+ktrk+kpt].Fill(d0Err)
                    h['h_ipdzErr'+ktrk+kpt].Fill(dzErr)
                
                    h['h_ipd0NoRefit'+ktrk+kpt].Fill(d0NoRefit)
                    h['h_ipdzNoRefit'+ktrk+kpt].Fill(dzNoRefit)

    # PV study
    for i in range(nEvents):
        
        tr.GetEntry(i)
        
        isValid = tr.pv_IsValid
        isFake = tr.pv_IsFake
            
        if not isValid or isFake: continue

        x = tr.pv_x
        y = tr.pv_y
        z = tr.pv_z
        xError = tr.pv_xError
        yError = tr.pv_yError
        zError = tr.pv_zError
        nTracks = tr.pv_NTracks

        x1 = tr.pv_x_p1
        y1 = tr.pv_y_p1
        z1 = tr.pv_z_p1
        xError1 = tr.pv_xError_p1
        yError1 = tr.pv_yError_p1
        zError1 = tr.pv_zError_p1
            
        x2 = tr.pv_x_p2
        y2 = tr.pv_y_p2
        z2 = tr.pv_z_p2
        xError2 = tr.pv_xError_p2
        yError2 = tr.pv_yError_p2
        zError2 = tr.pv_zError_p2
            
        dx12 = (x1-x2)/math.sqrt(2)
        dy12 = (y1-y2)/math.sqrt(2)
        dz12 = (z1-z2)/math.sqrt(2)
        
        dxPull12 = (x1-x2)/math.sqrt(xError1*xError1+xError2*xError2)
        dyPull12 = (y1-y2)/math.sqrt(yError1*yError1+yError2*yError2)
        dzPull12 = (z1-z2)/math.sqrt(zError1*zError1+zError2*zError2)

        cm = 1/10000.
        mm = 1/1000.

        h['h_pvNTrks'].Fill(nTracks)
            
        for k, v in c.PVnTracks.iteritems():
                
            nTrkMin = v[1]
            nTrkMax = v[2]
                
            if nTracks >= nTrkMin and nTracks < nTrkMax:
                    
                h['h_pvx'+k].Fill(x*mm)
                h['h_pvy'+k].Fill(y*mm)
                h['h_pvz'+k].Fill(z*cm)
                    
                h2['h2_pvx_y'+k].Fill(x*mm,y*mm)
                h2['h2_pvx_z'+k].Fill(z*cm,x*mm)
                h2['h2_pvy_z'+k].Fill(z*cm,y*mm)
                    
                h['h_pvdx12'+k].Fill(dx12)
                h['h_pvdy12'+k].Fill(dy12)
                h['h_pvdz12'+k].Fill(dz12)
                
                h['h_pvdxPull12'+k].Fill(dxPull12)
                h['h_pvdyPull12'+k].Fill(dyPull12)
                h['h_pvdzPull12'+k].Fill(dzPull12)
                
    print '\033[1;32mdone\033[1;m'

    outFile.Write()
    outFile.Close()
