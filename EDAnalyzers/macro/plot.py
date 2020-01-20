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
    parser.add_option("-p","--param",default="PVnTracks",help="parameterisation for PV resolution measurement [default: %default]")

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
    hd['pvSumTrackPt'] = {'xtit':'Sum of track p_{T}^2','nb':30,'xmin':0.,'xmax':150.,'ytit':'Events'}
    
    hd['ipPt'] = {'xtit':'Track p_{T} [GeV]','nb':100,'xmin':0.,'xmax':5.,'ytit':'Events'}
    hd['ipEta'] = {'xtit':'Track #eta','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
    hd['ipPhi'] = {'xtit':'Track #phi','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}

    PVparam = c.PVnTracks
    if 'PVsumTrackPt' in options.param: PVparam = c.PVsumTrackPt
    
    for k, v in PVparam.iteritems():
            
        for kk, vv in c.IPpt.iteritems():
    
            nBins = 100
            xMin = -1500
            xMax = 1500
            if v[1] > 1: nBins = 150
            
            hd['ippvd0'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            hd['ippvdz'+k+kk] = {'xtit':'d_{z} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            hd['ipbsd0'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            hd['ipbsdz'+k+kk] = {'xtit':'d_{z} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            
            hd['ipd0Err'+k+kk] = {'xtit':'#sigma(d_{xy}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            hd['ipdzErr'+k+kk] = {'xtit':'#sigma(d_{z}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            
            hd['ippvd0NoRefit'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}
            hd['ippvdzNoRefit'+k+kk] = {'xtit':'d_{z} [#mum]','nb':nBins,'xmin':xMin,'xmax':xMax,'ytit':'Events'}

    for k, v in PVparam.iteritems():
            
        for kk, vv in c.IPeta.iteritems():
        
            hd['ippvd0'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ippvdz'+k+kk] = {'xtit':'d_{z} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipbsd0'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ipbsdz'+k+kk] = {'xtit':'d_{z} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            
            hd['ipd0Err'+k+kk] = {'xtit':'#sigma(d_{xy}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            hd['ipdzErr'+k+kk] = {'xtit':'#sigma(d_{z}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
            
            hd['ippvd0NoRefit'+k+kk] = {'xtit':'d_{xy} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}
            hd['ippvdzNoRefit'+k+kk] = {'xtit':'d_{z} [#mum]','nb':100,'xmin':-3.,'xmax':3.,'ytit':'Events'}

    for k, v in PVparam.iteritems():

        if isData:
            hd['pvx'+k] = {'xtit':'x [mm]','nb':30,'xmin':0.5,'xmax':1.0,'ytit':'Events'}
            hd['pvy'+k] = {'xtit':'y [mm]','nb':30,'xmin':0.3,'xmax':0.9,'ytit':'Events'}
            hd['pvz'+k] = {'xtit':'z [cm]','nb':30,'xmin':-20.,'xmax':20.,'ytit':'Events'}
            
            h2d['pvx_y'+k] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':60,'xmin':0.5,'xmax':1.0,'nby':60,'ymin':0.3,'ymax':0.9}
            h2d['pvx_z'+k] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':60,'xmin':-20.,'xmax':20.,'nby':60,'ymin':0.5,'ymax':1.0}
            h2d['pvy_z'+k] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':60,'xmin':-20.,'xmax':20.,'nby':60,'ymin':0.3,'ymax':0.9}

        else:
            hd['pvx'+k] = {'xtit':'x [mm]','nb':30,'xmin':2.1,'xmax':2.8,'ytit':'Events'}
            hd['pvy'+k] = {'xtit':'y [mm]','nb':30,'xmin':3.6,'xmax':4.3,'ytit':'Events'}
            hd['pvz'+k] = {'xtit':'z [cm]','nb':30,'xmin':-25.,'xmax':25.,'ytit':'Events'}
            
            h2d['pvx_y'+k] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':60,'xmin':2.1,'xmax':2.8,'nby':60,'ymin':3.6,'ymax':4.3}
            h2d['pvx_z'+k] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':60,'xmin':-25.,'xmax':25.,'nby':60,'ymin':2.1,'ymax':2.8}
            h2d['pvy_z'+k] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':60,'xmin':-25.,'xmax':25.,'nby':60,'ymin':3.6,'ymax':4.3}
        
        hd['pvdx12'+k] = {'xtit':'Primary vertex resolution in x [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}
        hd['pvdy12'+k] = {'xtit':'Primary vertex resolution in y [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}
        hd['pvdz12'+k] = {'xtit':'Primary vertex resolution in z [#mum]','nb':30,'xmin':-200.,'xmax':200.,'ytit':'Events'}

        hd['pvdxPull12'+k] = {'xtit':'Primary vertex pull in x','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}
        hd['pvdyPull12'+k] = {'xtit':'Primary vertex pull in y','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}
        hd['pvdzPull12'+k] = {'xtit':'Primary vertex pull in z','nb':30,'xmin':-4,'xmax':4,'ytit':'Events'}

    for k, v in PVparam.iteritems():

        if isData:
            hd['bsx0'+k] = {'xtit':'x0 [mm]','nb':30,'xmin':0.5,'xmax':1.0,'ytit':'Events'}
            hd['bsy0'+k] = {'xtit':'y0 [mm]','nb':30,'xmin':0.3,'xmax':0.9,'ytit':'Events'}
            hd['bsz0'+k] = {'xtit':'z0 [cm]','nb':30,'xmin':-20.,'xmax':20.,'ytit':'Events'}

            hd['bsx'+k] = {'xtit':'x [mm]','nb':30,'xmin':0.5,'xmax':1.0,'ytit':'Events'}
            hd['bsy'+k] = {'xtit':'y [mm]','nb':30,'xmin':0.3,'xmax':0.9,'ytit':'Events'}
            
            h2d['bsx0_y0'+k] = {'xtit':'x0 [mm]','ytit':'y0 [mm]','nbx':60,'xmin':0.5,'xmax':1.0,'nby':60,'ymin':0.3,'ymax':0.9}
            h2d['bsx0_z0'+k] = {'xtit':'z0 [cm]','ytit':'x0 [mm]','nbx':60,'xmin':-20.,'xmax':20.,'nby':60,'ymin':0.5,'ymax':1.0}
            h2d['bsy0_z0'+k] = {'xtit':'z0 [cm]','ytit':'y0 [mm]','nbx':60,'xmin':-20.,'xmax':20.,'nby':60,'ymin':0.3,'ymax':0.9}

            h2d['bsx_y'+k] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':60,'xmin':0.5,'xmax':1.0,'nby':60,'ymin':0.3,'ymax':0.9}
            h2d['bsx_z'+k] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':60,'xmin':-20.,'xmax':20.,'nby':60,'ymin':0.5,'ymax':1.0}
            h2d['bsy_z'+k] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':60,'xmin':-20.,'xmax':20.,'nby':60,'ymin':0.3,'ymax':0.9}
            
        else:
            hd['bsx0'+k] = {'xtit':'x0 [mm]','nb':30,'xmin':2.1,'xmax':2.8,'ytit':'Events'}
            hd['bsy0'+k] = {'xtit':'y0 [mm]','nb':30,'xmin':3.6,'xmax':4.3,'ytit':'Events'}
            hd['bsz0'+k] = {'xtit':'z0 [cm]','nb':30,'xmin':-25.,'xmax':25.,'ytit':'Events'}

            hd['bsx'+k] = {'xtit':'x [mm]','nb':30,'xmin':2.1,'xmax':2.8,'ytit':'Events'}
            hd['bsy'+k] = {'xtit':'y [mm]','nb':30,'xmin':3.6,'xmax':4.3,'ytit':'Events'}
            
            h2d['bsx0_y0'+k] = {'xtit':'x0 [mm]','ytit':'y0 [mm]','nbx':60,'xmin':2.1,'xmax':2.8,'nby':60,'ymin':3.6,'ymax':4.3}
            h2d['bsx0_z0'+k] = {'xtit':'z0 [cm]','ytit':'x0 [mm]','nbx':60,'xmin':-25.,'xmax':25.,'nby':60,'ymin':2.1,'ymax':2.8}
            h2d['bsy0_z0'+k] = {'xtit':'z0 [cm]','ytit':'y0 [mm]','nbx':60,'xmin':-25.,'xmax':25.,'nby':60,'ymin':3.6,'ymax':4.3}

            h2d['bsx_y'+k] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':60,'xmin':2.1,'xmax':2.8,'nby':60,'ymin':3.6,'ymax':4.3}
            h2d['bsx_z'+k] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':60,'xmin':-25.,'xmax':25.,'nby':60,'ymin':2.1,'ymax':2.8}
            h2d['bsy_z'+k] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':60,'xmin':-25.,'xmax':25.,'nby':60,'ymin':3.6,'ymax':4.3}
            
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

        if not tr.trig_ZeroBias_pass and not tr.trig_ZeroBiasPixel_DoubleTrack_pass: continue
        
        isValid = tr.pv_IsValid
        isFake = tr.pv_IsFake
            
        if not isValid or isFake: continue

        nPVTracks = tr.pv_NTracks
        sumPtTrackPV = tr.pv_SumTrackPt

        param = nPVTracks
        if 'PVsumTrackPt' in options.param: param = sumPtTrackPV

        nTracks = tr.trk_pt.size()
            
        for t in range(nTracks):
        
            hasPXL = tr.trk_hasPXL[t]
            quality = tr.trk_quality[t]
            
            pt = tr.trk_pt[t]
            eta = tr.trk_eta[t]
            phi = tr.trk_phi[t]

            d0_pv = tr.trk_d0_pv[t]
            dz_pv = tr.trk_dz_pv[t]
            d0_bs = tr.trk_d0_bs[t]
            dz_bs = tr.trk_dz_bs[t]
            
            d0Err = tr.trk_d0Err[t]
            dzErr = tr.trk_dzErr[t]

            d0NoRefit = tr.trk_d0_pv_NoRefit[t]
            dzNoRefit = tr.trk_dz_pv_NoRefit[t]
            
            if math.fabs(d0_pv) > 2000 or math.fabs(dz_pv) > 2000: continue

            h['h_ipPt'].Fill(pt)
            h['h_ipEta'].Fill(eta)
            h['h_ipPhi'].Fill(phi)

            for kparam, vparam in PVparam.iteritems():
                
                paramMin = vparam[1]
                paramMax = vparam[2]
                
                if not (param >= paramMin and param < paramMax): continue
            
                for keta, veta in c.IPeta.iteritems():
                        
                    etaMin = veta[1]
                    etaMax = veta[2]
                
                    if not (eta >= etaMin and eta < etaMax): continue

                    h['h_ippvd0'+kparam+keta].Fill(d0_pv)
                    h['h_ippvdz'+kparam+keta].Fill(dz_pv)
                    h['h_ipbsd0'+kparam+keta].Fill(d0_bs)
                    h['h_ipbsdz'+kparam+keta].Fill(dz_bs)
                    
                    h['h_ipd0Err'+kparam+keta].Fill(d0Err)
                    h['h_ipdzErr'+kparam+keta].Fill(dzErr)
                    
                    h['h_ippvd0NoRefit'+kparam+keta].Fill(d0NoRefit)
                    h['h_ippvdzNoRefit'+kparam+keta].Fill(dzNoRefit)
                
                for kpt, vpt in c.IPpt.iteritems():
                    
                    ptMin = vpt[1]
                    ptMax = vpt[2]
                
                    if not (pt >= ptMin and pt < ptMax): continue
                    
                    h['h_ippvd0'+kparam+kpt].Fill(d0_pv)
                    h['h_ippvdz'+kparam+kpt].Fill(dz_pv)
                    h['h_ipbsd0'+kparam+kpt].Fill(d0_bs)
                    h['h_ipbsdz'+kparam+kpt].Fill(dz_bs)
                    
                    h['h_ipd0Err'+kparam+kpt].Fill(d0Err)
                    h['h_ipdzErr'+kparam+kpt].Fill(dzErr)
                
                    h['h_ippvd0NoRefit'+kparam+kpt].Fill(d0NoRefit)
                    h['h_ippvdzNoRefit'+kparam+kpt].Fill(dzNoRefit)

    # PV/BS study
    for i in range(nEvents):
        
        tr.GetEntry(i)
        
        if not tr.trig_ZeroBias_pass and not tr.trig_ZeroBiasPixel_DoubleTrack_pass: continue
        
        isValid = tr.pv_IsValid
        isFake = tr.pv_IsFake
            
        if not isValid or isFake: continue

        nPVTracks = tr.pv_NTracks
        sumPtTrackPV = tr.pv_SumTrackPt

        param = nPVTracks
        if 'PVsumTrackPt' in options.param: param = sumPtTrackPV
        
        bs_x0 = tr.bs_x0
        bs_y0 = tr.bs_y0
        bs_z0 = tr.bs_z0

        bs_x = tr.bs_x_zpv
        bs_y = tr.bs_y_zpv
        
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
        
        h['h_pvNTrks'].Fill(nPVTracks)
        h['h_pvSumTrackPt'].Fill(sumPtTrackPV)
            
        for k, v in PVparam.iteritems():
                
            paramMin = v[1]
            paramMax = v[2]
                
            if param >= paramMin and param < paramMax:
                    
                h['h_pvx'+k].Fill(pv_x*pv_mm)
                h['h_pvy'+k].Fill(pv_y*pv_mm)
                h['h_pvz'+k].Fill(pv_z*pv_cm)

                h2['h2_pvx_y'+k].Fill(pv_x*pv_mm,pv_y*pv_mm)
                h2['h2_pvx_z'+k].Fill(pv_z*pv_cm,pv_x*pv_mm)
                h2['h2_pvy_z'+k].Fill(pv_z*pv_cm,pv_y*pv_mm)
                    
                h['h_pvdx12'+k].Fill(pv_dx12*pv_mm)
                h['h_pvdy12'+k].Fill(pv_dy12*pv_mm)
                h['h_pvdz12'+k].Fill(pv_dz12*pv_cm)
                
                h['h_pvdxPull12'+k].Fill(pv_dxPull12)
                h['h_pvdyPull12'+k].Fill(pv_dyPull12)
                h['h_pvdzPull12'+k].Fill(pv_dzPull12)

                h['h_bsx0'+k].Fill(bs_x0*bs_mm)
                h['h_bsy0'+k].Fill(bs_y0*bs_mm)
                h['h_bsz0'+k].Fill(bs_z0*bs_cm)

                h['h_bsx'+k].Fill(bs_x*bs_mm)
                h['h_bsy'+k].Fill(bs_y*bs_mm)
                
                h2['h2_bsx0_y0'+k].Fill(bs_x0*bs_mm,bs_y0*bs_mm)
                h2['h2_bsx0_z0'+k].Fill(bs_z0*bs_cm,bs_x0*bs_mm)
                h2['h2_bsy0_z0'+k].Fill(bs_z0*bs_cm,bs_y0*bs_mm)

                h2['h2_bsx_y'+k].Fill(bs_x*bs_mm,bs_y*bs_mm)
                h2['h2_bsx_z'+k].Fill(bs_z0*bs_cm,bs_x*bs_mm)
                h2['h2_bsy_z'+k].Fill(bs_z0*bs_cm,bs_y*bs_mm)
                
    print '\033[1;32mdone\033[1;m'

    outFile.Write()
    outFile.Close()
