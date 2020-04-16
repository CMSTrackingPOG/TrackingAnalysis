import os
import sys
import subprocess
import common as c
import utils
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

    ParamList = [c.IPpt,c.IPeta,c.IPphi,c.IPnpv,c.IPdr]
        
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
        
    nEvents = tr.GetEntries()
    print 'Run on ' + ('Data' if isData else 'MC')
    print 'Processed events = ' + str(nEvents)

    h = {}
    hd = {}
    h2d = {}

    hd['evNpv'] = {'xtit':'Number of primary vertices','nb':70,'xmin':0.,'xmax':70.,'ytit':'Events'}
    
    hd['pvNTrks'] = {'xtit':'Number of tracks','nb':50,'xmin':0.,'xmax':150.,'ytit':'Events'}
    hd['pvSumTrackPt'] = {'xtit':'Sum of track p_{T} [GeV]','nb':60,'xmin':0.,'xmax':150.,'ytit':'Events'}
    
    hd['ipPt'] = {'xtit':'Track p_{T} [GeV]','nb':100,'xmin':0.,'xmax':20.,'ytit':'Events'}
    hd['ipEta'] = {'xtit':'Track #eta','nb':100,'xmin':-2.5,'xmax':2.5,'ytit':'Events'}
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

    PVparam = c.PVnTracks
    if 'PVsumTrackPt' in options.param: PVparam = c.PVsumTrackPt

    for p in ParamList:

        for kk, vv in p.iteritems():

            IP = p[kk]
        
            for k, v in PVparam.iteritems():

                hd['ippvd0'+k+kk] = {'xtit':'d_{xy}(PV) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
                hd['ippvdz'+k+kk] = {'xtit':'d_{z}(PV) [#mum]','nb':IP['dz'][0],'xmin':IP['dz'][1],'xmax':IP['dz'][2],'ytit':'Events'}
                
                hd['ipd0Err'+k+kk] = {'xtit':'#sigma(d_{xy}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
                hd['ipdzErr'+k+kk] = {'xtit':'#sigma(d_{z}) [#mum]','nb':100,'xmin':0.,'xmax':1000.,'ytit':'Events'}
                
                hd['ippvd0NoRefit'+k+kk] = {'xtit':'d_{xy}(PV) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
                hd['ippvdzNoRefit'+k+kk] = {'xtit':'d_{z}(PV) [#mum]','nb':IP['dz'][0],'xmin':IP['dz'][1],'xmax':IP['dz'][2],'ytit':'Events'}
                
            hd['ipbsd0'+kk] = {'xtit':'d_{xy}(BS) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
            hd['ipbsd0zpv'+kk] = {'xtit':'d_{xy}(BS) [#mum]','nb':IP['d0'][0],'xmin':IP['d0'][1],'xmax':IP['d0'][2],'ytit':'Events'}
            hd['ipbsdz'+kk] = {'xtit':'d_{z}(BS) [#mum]','nb':IP['dz'][0],'xmin':IP['dz'][1],'xmax':IP['dz'][2],'ytit':'Events'}
        
    for k, v in PVparam.iteritems():

        if isData:
            PV = c.PVData
            hd['pvx'+k] = {'xtit':'x [mm]','nb':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'ytit':'Events'}
            hd['pvy'+k] = {'xtit':'y [mm]','nb':PV['y'][0],'xmin':PV['y'][1],'xmax':PV['y'][2],'ytit':'Events'}
            hd['pvz'+k] = {'xtit':'z [cm]','nb':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'ytit':'Events'}
            
            h2d['pvx_y'+k] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            h2d['pvx_z'+k] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['x'][0],'ymin':PV['x'][1],'ymax':PV['x'][2]}
            h2d['pvy_z'+k] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}

        else:
            PV = c.PVMC
            hd['pvx'+k] = {'xtit':'x [mm]','nb':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'ytit':'Events'}
            hd['pvy'+k] = {'xtit':'y [mm]','nb':PV['y'][0],'xmin':PV['y'][1],'xmax':PV['y'][2],'ytit':'Events'}
            hd['pvz'+k] = {'xtit':'z [cm]','nb':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'ytit':'Events'}
            
            h2d['pvx_y'+k] = {'xtit':'x [mm]','ytit':'y [mm]','nbx':PV['x'][0],'xmin':PV['x'][1],'xmax':PV['x'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}
            h2d['pvx_z'+k] = {'xtit':'z [cm]','ytit':'x [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['x'][0],'ymin':PV['x'][1],'ymax':PV['x'][2]}
            h2d['pvy_z'+k] = {'xtit':'z [cm]','ytit':'y [mm]','nbx':PV['z'][0],'xmin':PV['z'][1],'xmax':PV['z'][2],'nby':PV['y'][0],'ymin':PV['y'][1],'ymax':PV['y'][2]}

        hd['pvdx12'+k] = {'xtit':'Primary vertex resolution in x [#mum]','nb':v['resox'][0],'xmin':v['resox'][1],'xmax':v['resox'][2],'ytit':'Events'}
        hd['pvdy12'+k] = {'xtit':'Primary vertex resolution in y [#mum]','nb':v['resoy'][0],'xmin':v['resoy'][1],'xmax':v['resoy'][2],'ytit':'Events'}
        hd['pvdz12'+k] = {'xtit':'Primary vertex resolution in z [#mum]','nb':v['resoz'][0],'xmin':v['resoz'][1],'xmax':v['resoz'][2],'ytit':'Events'}

        hd['pvdxPull12'+k] = {'xtit':'Primary vertex pull in x','nb':v['pullx'][0],'xmin':v['pullx'][1],'xmax':v['pullx'][2],'ytit':'Events'}
        hd['pvdyPull12'+k] = {'xtit':'Primary vertex pull in y','nb':v['pully'][0],'xmin':v['pully'][1],'xmax':v['pully'][2],'ytit':'Events'}
        hd['pvdzPull12'+k] = {'xtit':'Primary vertex pull in z','nb':v['pullz'][0],'xmin':v['pullz'][1],'xmax':v['pullz'][2],'ytit':'Events'}

    if isData:
        PV = c.PVData
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

        BS = c.BSData
        hd['bsBeamWidthX'] = {'xtit':'Beam width (x) [#mum]','nb':BS['bwx'][0],'xmin':BS['bwx'][1],'xmax':BS['bwx'][2],'ytit':'Events'}
        hd['bsBeamWidthY'] = {'xtit':'Beam width (y) [#mum]','nb':BS['bwy'][0],'xmin':BS['bwy'][1],'xmax':BS['bwy'][2],'ytit':'Events'}
        hd['bsSigmaZ'] = {'xtit':'Beam sigma z [mm]','nb':BS['bwz'][0],'xmin':BS['bwz'][1],'xmax':BS['bwz'][2],'ytit':'Events'}
        
    else:
        PV = c.PVMC
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

        BS = c.BSMC
        hd['bsBeamWidthX'] = {'xtit':'Beam width (x) [#mum]','nb':BS['bwx'][0],'xmin':BS['bwx'][1],'xmax':BS['bwx'][2],'ytit':'Events'}
        hd['bsBeamWidthY'] = {'xtit':'Beam width (y) [#mum]','nb':BS['bwy'][0],'xmin':BS['bwy'][1],'xmax':BS['bwy'][2],'ytit':'Events'}
        hd['bsSigmaZ'] = {'xtit':'Beam sigma z [mm]','nb':BS['bwz'][0],'xmin':BS['bwz'][1],'xmax':BS['bwz'][2],'ytit':'Events'}
        
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

        if isDataZeroBias or isMCZeroBias:                
            if not tr.trig_ZeroBias_pass: continue

        if isDataJetHT or isMCJetHT:
            if not (tr.trig_PFJet40_pass or tr.trig_PFJet60_pass or \
            tr.trig_PFJet80_pass or tr.trig_PFJet140_pass or \
            tr.trig_PFJet200_pass or tr.trig_PFJet260_pass or \
            tr.trig_PFJet320_pass or tr.trig_PFJet400_pass or \
            tr.trig_PFJet450_pass or tr.trig_PFJet500_pass or tr.trig_PFJet550_pass): continue
            
        isValid = tr.pv_IsValid
        isFake = tr.pv_IsFake
            
        if not isValid or isFake: continue

        nPVTracks = tr.pv_NTracks
        sumPtTrackPV = tr.pv_SumTrackPt
        
        if nPVTracks < 3: continue
        
        param = nPVTracks
        if 'PVsumTrackPt' in options.param: param = sumPtTrackPV

        npv = tr.ev_nPV
        h['h_evNpv'].Fill(npv)
        
        bs_sigmaZ = tr.bs_sigmaZ*10 # mm
        
        nTracks = tr.trk_pt.size()

        for t in range(nTracks):
        
            hasPXL1 = (tr.trk_hasPixelBarrelLayer1[t] or tr.trk_hasPixelEndcapLayer1[t])
            hasPXL2 = (tr.trk_hasPixelBarrelLayer2[t] or tr.trk_hasPixelEndcapLayer2[t])
            hasPXL3 = (tr.trk_hasPixelBarrelLayer3[t] or tr.trk_hasPixelEndcapLayer3[t])
            hasPXL4 = (tr.trk_hasPixelBarrelLayer4[t] or tr.trk_hasPixelEndcapLayer4[t])
            
            trkSelPXL1 = bool(hasPXL1)
            trkSelPXL2 = bool(hasPXL1 and hasPXL2)
            trkSelPXL3 = bool(hasPXL1 and hasPXL2 and hasPXL3)
            trkSelPXL4 = bool(hasPXL1 and hasPXL2 and hasPXL3 and hasPXL4)
            
            quality = tr.trk_quality[t]
            
            pt = tr.trk_pt[t]
            eta = tr.trk_eta[t]
            phi = tr.trk_phi[t]

            isJet = tr.trk_jet_found[t]
            jetEta = tr.trk_jet_eta[t]
            jetPhi = tr.trk_jet_phi[t]
            nTrkJet = tr.trk_jet_nTracks[t]

            if isJet:
                drTrkJet = utils.deltaR2(eta,phi,jetEta,jetPhi)
            
            d0 = tr.trk_d0[t]
            dz = tr.trk_dz[t]

            d0_pv = tr.trk_d0_pv[t]
            dz_pv = tr.trk_dz_pv[t]
            d0_bs = tr.trk_d0_bs[t]
            dz_bs = tr.trk_dz_bs[t]

            d0_bs_zpv = tr.trk_d0_bs_zpv[t]
            d0_bs_zpca = tr.trk_d0_bs_zpca[t]
            
            d0Err = tr.trk_d0Err[t]
            dzErr = tr.trk_dzErr[t]

            d0NoRefit = tr.trk_d0_pv_NoRefit[t]
            dzNoRefit = tr.trk_dz_pv_NoRefit[t]

            if pt < 0.05: continue
            if math.fabs(eta) > 2.5: continue

            h['h_ipPt'].Fill(pt)
            h['h_ipEta'].Fill(eta)
            h['h_ipPhi'].Fill(phi)
            if isJet: h['h_ipDrTrkJet'].Fill(drTrkJet)
            if isJet: h['h_ipNTrkJet'].Fill(nTrkJet)

            sd0 = d0/d0Err if d0Err > 0 else -777
            sdz = dz/dzErr if dzErr > 0 else -777
            sd0_pv = d0_pv/d0Err if d0Err > 0 else -777
            sdz_pv = dz_pv/dzErr if dzErr > 0 else -777
            sd0_bs = d0_bs/d0Err if d0Err > 0 else -777
            sdz_bs = dz_bs/dzErr if dzErr > 0 else -777
            sd0_bs_zpv = d0_bs_zpv/d0Err if d0Err > 0 else -777
            sd0_bs_zpca = d0_bs_zpca/d0Err if d0Err > 0 else -777
            
            mm = 1E-3
            h['h_ipD0'].Fill(d0*mm)
            h['h_ipDz'].Fill(dz*mm)
            h['h_ipSD0'].Fill(sd0)
            h['h_ipSDz'].Fill(sdz)
            
            h['h_ippvD0'].Fill(d0_pv*mm)
            h['h_ippvDz'].Fill(dz_pv*mm)
            h['h_ippvSD0'].Fill(sd0_pv)
            h['h_ippvSDz'].Fill(sdz_pv)
            
            h['h_ipbsD0'].Fill(d0_bs*mm)
            h['h_ipbsDz'].Fill(dz_bs*mm)
            h['h_ipbsSD0'].Fill(sd0_bs)
            h['h_ipbsSDz'].Fill(sdz_bs)
            
            h['h_ipbszpvD0'].Fill(d0_bs_zpv*mm)
            h['h_ipbszpcaD0'].Fill(d0_bs_zpca*mm)
            h['h_ipbszpvSD0'].Fill(sd0_bs_zpv)
            h['h_ipbszpcaSD0'].Fill(sd0_bs_zpca)

            for pi, p in enumerate(ParamList):
                
                for kp, vp in p.iteritems():

                    pMin = vp['bins'][1]
                    pMax = vp['bins'][2]
                
                    varp = pt
                    if pi == 1: varp = eta
                    elif pi == 2: varp = phi
                    elif pi == 3: varp = npv
                    elif pi == 4:
                        if not isJet: continue
                        else: varp = drTrkJet
                    elif pi == 5: varp = bs_sigmaZ
                        
                    if not (varp >= pMin and varp < pMax): continue

                    h['h_ipbsd0zpv'+kp].Fill(d0_bs_zpv)
                    h['h_ipbsd0'+kp].Fill(d0_bs)
                    h['h_ipbsdz'+kp].Fill(dz_bs)
                
                    for kparam, vparam in PVparam.iteritems():
                
                        paramMin = vparam['bins'][1]
                        paramMax = vparam['bins'][2]
                
                        if not (param >= paramMin and param < paramMax): continue
                        
                        h['h_ippvd0'+kparam+kp].Fill(d0_pv)
                        h['h_ippvdz'+kparam+kp].Fill(dz_pv)
                        
                        h['h_ipd0Err'+kparam+kp].Fill(d0Err)
                        h['h_ipdzErr'+kparam+kp].Fill(dzErr)
                        
                        h['h_ippvd0NoRefit'+kparam+kp].Fill(d0NoRefit)
                        h['h_ippvdzNoRefit'+kparam+kp].Fill(dzNoRefit)
                    
    # PV/BS study
    for i in range(nEvents):
        
        tr.GetEntry(i)
        
        if not tr.trig_ZeroBias_pass: continue
        
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

#        print 'x0='+str(tr.bs_x0), ' x0Error='+str(tr.bs_x0Error), ' BeamWidthX='+str(tr.bs_BeamWidthX), ' z0='+str(tr.bs_z0), ' z0Error='+str(tr.bs_z0Error), ' sigmaZ='+str(tr.bs_sigmaZ)
        
        bs_x = tr.bs_x_zpv
        bs_y = tr.bs_y_zpv
        
        bs_beamWidthX = tr.bs_BeamWidthX
        bs_beamWidthY = tr.bs_BeamWidthY
        bs_sigmaZ = tr.bs_sigmaZ*10 # mm

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
        
        h['h_pvNTrks'].Fill(nPVTracks)
        h['h_pvSumTrackPt'].Fill(sumPtTrackPV)

        h['h_bsx0'].Fill(bs_x0*bs_mm)
        h['h_bsy0'].Fill(bs_y0*bs_mm)
        h['h_bsz0'].Fill(bs_z0*bs_cm)
        
        h['h_bsx'].Fill(bs_x*bs_mm)
        h['h_bsy'].Fill(bs_y*bs_mm)
        
        h2['h2_bsx0_y0'].Fill(bs_x0*bs_mm,bs_y0*bs_mm)
        h2['h2_bsx0_z0'].Fill(bs_z0*bs_cm,bs_x0*bs_mm)
        h2['h2_bsy0_z0'].Fill(bs_z0*bs_cm,bs_y0*bs_mm)
        
        h2['h2_bsx_y'].Fill(bs_x*bs_mm,bs_y*bs_mm)
        h2['h2_bsx_z'].Fill(bs_z0*bs_cm,bs_x*bs_mm)
        h2['h2_bsy_z'].Fill(bs_z0*bs_cm,bs_y*bs_mm)

        h['h_bsBeamWidthX'].Fill(bs_beamWidthX*bs_micron)
        h['h_bsBeamWidthY'].Fill(bs_beamWidthY*bs_micron)
        h['h_bsSigmaZ'].Fill(bs_sigmaZ*bs_cm)
        
        for k, v in PVparam.iteritems():
                
            paramMin = v['bins'][1]
            paramMax = v['bins'][2]
                
            if param >= paramMin and param < paramMax:
                    
                h['h_pvx'+k].Fill(pv_x*pv_mm)
                h['h_pvy'+k].Fill(pv_y*pv_mm)
                h['h_pvz'+k].Fill(pv_z*pv_cm)

                h2['h2_pvx_y'+k].Fill(pv_x*pv_mm,pv_y*pv_mm)
                h2['h2_pvx_z'+k].Fill(pv_z*pv_cm,pv_x*pv_mm)
                h2['h2_pvy_z'+k].Fill(pv_z*pv_cm,pv_y*pv_mm)
                    
                h['h_pvdx12'+k].Fill(pv_dx12)
                h['h_pvdy12'+k].Fill(pv_dy12)
                h['h_pvdz12'+k].Fill(pv_dz12)
                
                h['h_pvdxPull12'+k].Fill(pv_dxPull12)
                h['h_pvdyPull12'+k].Fill(pv_dyPull12)
                h['h_pvdzPull12'+k].Fill(pv_dzPull12)
                
    print '\033[1;32mdone\033[1;m'

    outFile.Write()
    outFile.Close()
