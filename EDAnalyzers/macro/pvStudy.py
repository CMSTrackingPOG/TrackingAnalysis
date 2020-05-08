import os
import sys
import math
import subprocess
import common as c
from subprocess import call
import style
import functions as func
import json
import ROOT

import fit

ROOT.PyConfig.IgnoreCommandLineOptions = True
from optparse import OptionParser

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]

    usage = "usage: %prog [options]\n Analysis script to study PV and BS resolution"
    
    parser = OptionParser(usage)
    parser.add_option("-d","--data",default="data.root",help="Input data file name [default: %default]")
    parser.add_option("-m","--mc",default="mc.root",help="Input mc file name [default: %default]")
    parser.add_option("-o","--output",default="pics",help="Output directory [default: %default]")
    parser.add_option("--qcd",action='store_true',help="Use QCD events [default: %default]")
    parser.add_option("-p","--param",default="PVsumTrackPtSq",help="Parameterisation for PV resolution measurement [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(1)

    if not os.path.isdir(options.output):
        os.system("mkdir "+options.output)
        
    fHistData = ROOT.TFile.Open(options.data,'read')
    fHistMC = ROOT.TFile.Open(options.mc,'read')
    
    evt = 'zb'
    if options.qcd: evt = 'qcd'

    param = c.PVnTracks[evt]
    if options.param == 'PVsumTrackPt': param = c.PVsumTrackPt[evt]
    elif options.param == 'PVsumTrackPtSq': param = c.PVsumTrackPtSq[evt]
    elif options.param != 'PVnTracks':
        print 'Unknown parameterisation requested:', options.param
        sys.exit()        
    
    hIncl = ['jetPtMax','jetHT','evNpv','pvNTrks','pvSumTrackPt','pvSumTrackPt2']

    print 'Plot cumulative distributions'
    for h in hIncl:
        
        hData = fHistData.Get('h_'+h)
        hMC = fHistMC.Get('h_'+h)
        
        c1 = ROOT.TCanvas()
        
        hMC.Draw('hist')
        hData.Draw('e1 same')
        intMC = hMC.Integral()
        intData = hData.Integral()
        hMC.Scale(intData/intMC)
        maxData = hData.GetMaximum()
        maxMC = hMC.GetMaximum()
        hMC.SetMaximum(1.2*max(maxData,maxMC))
        hMC.SetMinimum(0.)

        hData.SetMarkerSize(0.7)
        hData.SetMarkerColor(1)
        hData.SetLineColor(1)
        
        hMC.SetMarkerSize(0)
        hMC.SetMarkerColor(c.mccol)
        hMC.SetLineColor(c.mccol)
        hMC.SetFillColor(c.mccol)
        hMC.SetLineStyle(1)
        
        leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
        leg.SetFillColor(253)
        leg.SetBorderSize(0)
        leg.AddEntry(hData,"Data","p")
        leg.AddEntry(hMC,"Simulation","f")
        leg.Draw()        

        t1, t2, t3, t4 = style.cmslabel(1, c.year, evt)
        t1.Draw()
        t2.Draw()
        t3.Draw()
        t4.Draw()
        
        c1.Print(options.output+'/'+h+'.pdf')
        c1.Clear()
    
    #### Basic distributions

    print 'Plot PV and BS profiles'
    
    #### PV
    
    pv_hXData = {}; pv_hYData = {}; pv_hZData = {}
    pv_hXMC = {}; pv_hYMC = {}; pv_hZMC = {}
    pv_hXYData = {}; pv_hXZData = {}; pv_hYZData = {}
    pv_hXYMC = {}; pv_hXZMC = {}; pv_hYZMC = {}

    for k, v in param.iteritems():
        
        pv_hXData[k] = fHistData.Get('h_pvx'+k)
        pv_hYData[k] = fHistData.Get('h_pvy'+k)
        pv_hZData[k] = fHistData.Get('h_pvz'+k)

        pv_hXMC[k] = fHistMC.Get('h_pvx'+k)
        pv_hYMC[k] = fHistMC.Get('h_pvy'+k)
        pv_hZMC[k] = fHistMC.Get('h_pvz'+k)
        
        pv_hXYData[k] = fHistData.Get('h2_pvx_y'+k)
        pv_hXZData[k] = fHistData.Get('h2_pvx_z'+k)
        pv_hYZData[k] = fHistData.Get('h2_pvy_z'+k)
        
        pv_hXYMC[k] = fHistMC.Get('h2_pvx_y'+k)
        pv_hXZMC[k] = fHistMC.Get('h2_pvx_z'+k)
        pv_hYZMC[k] = fHistMC.Get('h2_pvy_z'+k)

    #### BS
        
    bs_hXData = {}; bs_hYData = {}; bs_hZData = {}
    bs_hXMC = {}; bs_hYMC = {}; bs_hZMC = {}
    bs_hXYData = {}; bs_hXZData = {}; bs_hYZData = {}
    bs_hXYMC = {}; bs_hXZMC = {}; bs_hYZMC = {}
    bs_hBeamWidthXData = {}; bs_hBeamWidthYData = {}; bs_hSigmaZData = {}
    bs_hBeamWidthXMC = {}; bs_hBeamWidthYMC = {}; bs_hSigmaZMC = {}
        
    bs_hXData[''] = fHistData.Get('h_bsx0')
    bs_hYData[''] = fHistData.Get('h_bsy0')
    bs_hZData[''] = fHistData.Get('h_bsz0')
    
    bs_hXMC[''] = fHistMC.Get('h_bsx0')
    bs_hYMC[''] = fHistMC.Get('h_bsy0')
    bs_hZMC[''] = fHistMC.Get('h_bsz0')
    
    bs_hXYData[''] = fHistData.Get('h2_bsx0_y0')
    bs_hXZData[''] = fHistData.Get('h2_bsx0_z0')
    bs_hYZData[''] = fHistData.Get('h2_bsy0_z0')
    
    bs_hXYMC[''] = fHistMC.Get('h2_bsx0_y0')
    bs_hXZMC[''] = fHistMC.Get('h2_bsx0_z0')
    bs_hYZMC[''] = fHistMC.Get('h2_bsy0_z0')
    
    bs_hBeamWidthXData[''] = fHistData.Get('h_bsBeamWidthX')
    bs_hBeamWidthYData[''] = fHistData.Get('h_bsBeamWidthY')
    bs_hSigmaZData[''] = fHistData.Get('h_bsSigmaZ')

    bs_hBeamWidthXMC[''] = fHistMC.Get('h_bsBeamWidthX')
    bs_hBeamWidthYMC[''] = fHistMC.Get('h_bsBeamWidthY')
    bs_hSigmaZMC[''] = fHistMC.Get('h_bsSigmaZ')

    for idx, h in enumerate([pv_hXData,pv_hYData,pv_hZData,pv_hXMC,pv_hYMC,pv_hZMC]):
        for k, hh in h.iteritems():
        
            c1 = ROOT.TCanvas()

            hh.Draw('e1')

            proj = 'X'
            if 'y' in hh.GetName(): proj = 'Y'
            if 'z' in hh.GetName(): proj = 'Z'
            samp = 'data'
            if idx > 2: samp = 'mc'
            
            c1.Print(options.output+'/pv'+proj+k+'_'+samp+'.pdf')
            c1.Clear()

    for idx, h2 in enumerate([pv_hXYData,pv_hXZData,pv_hYZData,pv_hXYMC,pv_hXZMC,pv_hYZMC]):
        for k, hh in h2.iteritems():
        
            c1 = ROOT.TCanvas()
            hh.Draw('COLZ')
            hh.SetMinimum(3)
            
            proj = 'XY'
            if 'x_z' in hh.GetName(): proj = 'XZ'
            if 'y_z' in hh.GetName(): proj = 'YZ'
            samp = 'data'
            if idx > 2: samp = 'mc'

            if 'x_y' in hh.GetName() and k == '':
                if samp == 'data':
                    hpvplot = pv_hXYData[k]
                    hbsplot = bs_hXYData[k]
                    hbsplot.Draw('BOX SAME')
                    hbsplot.SetLineColor(ROOT.kRed)
                else:
                    hpvplot = pv_hXYMC[k]
                    hbsplot = bs_hXYMC[k]
                    hbsplot.Draw('BOX SAME')
                    hbsplot.SetLineColor(ROOT.kRed)
                lpvrms1 = ROOT.TLatex(0.2,0.3,"RMS (PV)"); lpvrms1.SetNDC(); lpvrms1.SetTextFont(43); lpvrms1.SetTextSize(20); lpvrms1.Draw()
                lpvrms2 = ROOT.TLatex(0.2,0.25,"x = %.1f" % (hpvplot.GetRMS(1)*math.pow(10,3))+" #mum"); lpvrms2.SetNDC(); lpvrms2.SetTextFont(43); lpvrms2.SetTextSize(20); lpvrms2.Draw()
                lpvrms3 = ROOT.TLatex(0.2,0.2,"y = %.1f" % (hpvplot.GetRMS(2)*math.pow(10,3))+" #mum"); lpvrms3.SetNDC(); lpvrms3.SetTextFont(43); lpvrms3.SetTextSize(20); lpvrms3.Draw()
                if samp == 'data':
                    lbsrms1 = ROOT.TLatex(0.65,0.3,"RMS (BS)"); lbsrms1.SetNDC(); lbsrms1.SetTextFont(43); lbsrms1.SetTextSize(20); lbsrms1.Draw()
                    lbsrms2 = ROOT.TLatex(0.65,0.25,"x = %.1f" % (hbsplot.GetRMS(1)*math.pow(10,3))+" #mum"); lbsrms2.SetNDC(); lbsrms2.SetTextFont(43); lbsrms2.SetTextSize(20); lbsrms2.Draw()
                    lbsrms3 = ROOT.TLatex(0.65,0.2,"y = %.1f" % (hbsplot.GetRMS(2)*math.pow(10,3))+" #mum"); lbsrms3.SetNDC(); lbsrms3.SetTextFont(43); lbsrms3.SetTextSize(20); lbsrms3.Draw()
                    
            if 'x_z' in hh.GetName() and k == '': 
                if samp == 'data':
                    hpvplot = pv_hXZData[k]
                    hbsplot = bs_hXZData[k]
                    hbsplot.Draw('BOX SAME')
                    hbsplot.SetLineColor(ROOT.kRed)
                else:
                    hpvplot = pv_hXZMC[k]
                    hbsplot = bs_hXZMC[k]
                    hbsplot.Draw('BOX SAME')
                    hbsplot.SetLineColor(ROOT.kRed)
                lpvrms1 = ROOT.TLatex(0.2,0.3,"RMS (PV)"); lpvrms1.SetNDC(); lpvrms1.SetTextFont(43); lpvrms1.SetTextSize(20); lpvrms1.Draw()
                lpvrms2 = ROOT.TLatex(0.2,0.25,"x = %.1f" % (hpvplot.GetRMS(2)*math.pow(10,3))+" #mum"); lpvrms2.SetNDC(); lpvrms2.SetTextFont(43); lpvrms2.SetTextSize(20); lpvrms2.Draw()
                lpvrms3 = ROOT.TLatex(0.2,0.2,"z = %.1f" % (hpvplot.GetRMS(1))+" mm"); lpvrms3.SetNDC(); lpvrms3.SetTextFont(43); lpvrms3.SetTextSize(20); lpvrms3.Draw()
                if samp == 'data':
                    lbsrms1 = ROOT.TLatex(0.65,0.3,"RMS (BS)"); lbsrms1.SetNDC(); lbsrms1.SetTextFont(43); lbsrms1.SetTextSize(20); lbsrms1.Draw()
                    lbsrms2 = ROOT.TLatex(0.65,0.25,"x = %.1f" % (hbsplot.GetRMS(2)*math.pow(10,3))+" #mum"); lbsrms2.SetNDC(); lbsrms2.SetTextFont(43); lbsrms2.SetTextSize(20); lbsrms2.Draw()
                    lbsrms3 = ROOT.TLatex(0.65,0.2,"z = %.1f" % (hbsplot.GetRMS(1))+" mm"); lbsrms3.SetNDC(); lbsrms3.SetTextFont(43); lbsrms3.SetTextSize(20); lbsrms3.Draw()
                    
            if 'y_z' in hh.GetName() and k == '':
                if samp == 'data':
                    hpvplot = pv_hYZData[k]
                    hbsplot = bs_hYZData[k]
                    hbsplot.Draw('BOX SAME')
                    hbsplot.SetLineColor(ROOT.kRed)
                else:
                    hpvplot = pv_hYZMC[k]
                    hbsplot = bs_hYZMC[k]
                    hbsplot.Draw('BOX SAME')
                    hbsplot.SetLineColor(ROOT.kRed)
                lpvrms1 = ROOT.TLatex(0.2,0.3,"RMS (PV)"); lpvrms1.SetNDC(); lpvrms1.SetTextFont(43); lpvrms1.SetTextSize(20); lpvrms1.Draw()
                lpvrms2 = ROOT.TLatex(0.2,0.25,"y = %.1f" % (hpvplot.GetRMS(2)*math.pow(10,3))+" #mum"); lpvrms2.SetNDC(); lpvrms2.SetTextFont(43); lpvrms2.SetTextSize(20); lpvrms2.Draw()
                lpvrms3 = ROOT.TLatex(0.2,0.2,"z = %.1f" % (hpvplot.GetRMS(1))+" mm"); lpvrms3.SetNDC(); lpvrms3.SetTextFont(43); lpvrms3.SetTextSize(20); lpvrms3.Draw()
                if samp == 'data':
                    lbsrms1 = ROOT.TLatex(0.65,0.3,"RMS (BS)"); lbsrms1.SetNDC(); lbsrms1.SetTextFont(43); lbsrms1.SetTextSize(20); lbsrms1.Draw()
                    lbsrms2 = ROOT.TLatex(0.65,0.25,"y = %.1f" % (hbsplot.GetRMS(2)*math.pow(10,3))+" #mum"); lbsrms2.SetNDC(); lbsrms2.SetTextFont(43); lbsrms2.SetTextSize(20); lbsrms2.Draw()
                    lbsrms3 = ROOT.TLatex(0.65,0.2,"z = %.1f" % (hbsplot.GetRMS(1))+" mm"); lbsrms3.SetNDC(); lbsrms3.SetTextFont(43); lbsrms3.SetTextSize(20); lbsrms3.Draw()

            t1, t2, t3, t4 = style.cmslabel(1, c.year, evt)
            t1.Draw()
            t2.Draw()
            t3.Draw()
            t4.Draw()
                
            c1.Print(options.output+'/pv'+proj+k+'_'+samp+'.pdf')
            c1.Clear()

    for idx, h in enumerate([bs_hBeamWidthXData,bs_hBeamWidthYData,bs_hSigmaZData,bs_hBeamWidthXMC,bs_hBeamWidthYMC,bs_hSigmaZMC]):
        for k, hh in h.iteritems():
        
            c1 = ROOT.TCanvas()

            hh.Draw('hist')

            proj = 'X'
            if 'Y' in hh.GetName(): proj = 'Y'
            if 'Z' in hh.GetName(): proj = 'Z'
            samp = 'data'
            if idx > 2: samp = 'mc'

            t1, t2, t3, t4 = style.cmslabel(1, c.year, evt)
            t1.Draw()
            t2.Draw()
            t3.Draw()
            t4.Draw()
            
            c1.Print(options.output+'/bs'+proj+k+'_'+samp+'.pdf')
            c1.Clear()

    beamWidthXData = bs_hBeamWidthXData[''].GetMean()
    beamWidthYData = bs_hBeamWidthYData[''].GetMean()
    sigmaZData = bs_hSigmaZData[''].GetMean()

    beamWidthXMC = bs_hBeamWidthXMC[''].GetMean()
    beamWidthYMC = bs_hBeamWidthYMC[''].GetMean()
    sigmaZMC = bs_hSigmaZMC[''].GetMean()
    
    #### Fits
    
    print 'Run PV fits'
    
    pstyle.SetOptFit(1111)
    
    rout = {}
    for l1 in ['reso','pull','widthx','widthy','sigmaz']: rout[l1] = {}
    for t in ['data','mc']:
        for l2 in ['reso','pull','widthx','widthy','sigmaz']: rout[l2][t] = {}
        for x in c.PVmeas:
            for l3 in ['reso','pull','widthx','widthy','sigmaz']: rout[l3][t][x] = {}
            for k, v in param.iteritems():
                for l4 in ['reso','pull','widthx','widthy','sigmaz']: rout[l4][t][x][k] = {}
                
    for k, v in param.iteritems():
    
        for x in c.PVmeas:
            
            hNameResoData = 'h_pvd'+x+'12'+k
            hNamePullData = 'h_pvd'+x+'Pull12'+k
            hResoData = fHistData.Get(hNameResoData)
            hPullData = fHistData.Get(hNamePullData)
            func.addbin(hResoData)
            func.addbin(hPullData)
        
            hNameResoMC = 'h_pvd'+x+'12'+k
            hNamePullMC = 'h_pvd'+x+'Pull12'+k
            hResoMC = fHistMC.Get(hNameResoMC)
            hPullMC = fHistMC.Get(hNamePullMC)
            func.addbin(hResoMC)
            func.addbin(hPullMC)
            
            if hResoData.GetEntries() < 10:
                print 'No stats in data: '+hResoData.GetName()
                sys.exit()
            if hResoMC.GetEntries() < 10:
                print 'No stats in mc: '+hResoMC.GetName()
                sys.exit()

            for h in [hResoData,hPullData]:
                h.SetMarkerSize(0.7)
                h.SetMarkerColor(1)
                h.SetLineColor(1)

            for h in [hResoMC,hPullMC]:
                h.SetMarkerSize(0)
                h.SetMarkerColor(c.mccol)
                h.SetLineColor(c.mccol)
                h.SetFillColor(c.mccol)
                h.SetLineStyle(1)
        
            intResoMC = hResoMC.Integral()
            intResoData = hResoData.Integral()
            hResoMC.Scale(intResoData/intResoMC)        

            intPullMC = hPullMC.Integral()
            intPullData = hPullData.Integral()
            hPullMC.Scale(intPullData/intPullMC)        
            
            maxResoData = hResoData.GetMaximum()
            maxResoMC = hResoMC.GetMaximum()
            hResoMC.SetMaximum(1.2*max(maxResoData,maxResoMC))
            hResoMC.SetMinimum(0.)

            maxPullData = hPullData.GetMaximum()
            maxPullMC = hPullMC.GetMaximum()
            hPullMC.SetMaximum(1.2*max(maxPullData,maxPullMC))
            hPullMC.SetMinimum(0.)

            selName = 'N_{trk}'
            units = ''
            if options.param == 'PVsumTrackPt': 
                selName = '#sump_{T}'
                units = ' GeV'
            elif options.param == 'PVsumTrackPtSq': 
                selName = '#sqrt{#sump^{2}_{T}}'
                units = ' GeV'
            
            # Resolution
        
            c1 = ROOT.TCanvas()
        
            hResoMC.Draw('hist')
            hResoData.Draw('e1 sames')
        
            fitt = '2g'
            if '0to2' in k and evt == 'zb': fitt = ''

            resResoMC, resoMC, resoErrMC, resoChi2MC = fit.doFit('mcfit',hResoMC,x,k,c.mcfit,fitt)
            resResoMC.Draw("same")
            
            resResoData, resoData, resoErrData, resoChi2Data = fit.doFit('datafit',hResoData,x,k,1,fitt)
            resResoData.Draw("same")
            resResoData.SetLineStyle(2)
            
            c1.Update()
            
            finfoData = hResoData.GetListOfFunctions().FindObject("stats")
            finfoData.__class__ = ROOT.TPaveStats
            finfoData.SetX1NDC(0.7)
            finfoData.SetX2NDC(0.95)
            finfoData.SetY1NDC(0.20)
            finfoData.SetY2NDC(0.40)
            lData = ROOT.TText(0.81,0.41,"Data (fit)")
            lData.SetTextSize(0.035)
            lData.SetNDC()
            lData.Draw()
            
            finfoMC = hResoMC.GetListOfFunctions().FindObject("stats")
            finfoMC.__class__ = ROOT.TPaveStats
            finfoMC.SetX1NDC(0.7)
            finfoMC.SetX2NDC(0.95)
            finfoMC.SetY1NDC(0.45)
            finfoMC.SetY2NDC(0.65)
            lMC = ROOT.TText(0.81,0.66,"Simulation (fit)")
            lMC.SetTextSize(0.035)
            lMC.SetNDC()
            lMC.Draw()

            lDataReso = ROOT.TLatex(0.20,0.63,"#sigma^{Data}_{"+x+"} = %.1f #mum" % (resoData))
            lDataReso.SetNDC(); lDataReso.SetTextFont(43); lDataReso.SetTextSize(20); lDataReso.Draw()

            lMCReso = ROOT.TLatex(0.20,0.70,"#sigma^{Sim.}_{"+x+"} = %.1f #mum" % (resoMC))
            lMCReso.SetNDC(); lMCReso.SetTextFont(43); lMCReso.SetTextSize(20); lMCReso.Draw()
            
            lSel = ROOT.TLatex(0.20,0.80,str(v['bins'][1])+' < '+selName+' < '+str(v['bins'][2])+units)
            lSel.SetTextSize(0.035)
            lSel.SetNDC()
            lSel.Draw()
            
            c1.Update()

            leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
            leg.SetFillColor(253)
            leg.SetBorderSize(0)
            leg.AddEntry(hResoData,"Data","p")
            leg.AddEntry(hResoMC,"Simulation","f")
            leg.AddEntry(resResoData,"Data (fit)","l")
            leg.AddEntry(resResoMC,"Simulation (fit)","l")
            leg.Draw()        
            
            t1, t2, t3, t4 = style.cmslabel(1, c.year, evt)
            t1.Draw()
            t2.Draw()
            t3.Draw()
            t4.Draw()
            
            c1.Print(options.output+'/pvReso_'+x+k+'.pdf')
            c1.Clear()

            # Pull
        
            hPullMC.Draw('hist')
            hPullData.Draw('e1 sames')
            
            resPullMC, pullMC, pullErrMC, pullChi2MC = fit.doFit('mcfit',hPullMC,x,k,c.mcfit,'2g')
            resPullMC.Draw("same")
            
            resPullData, pullData, pullErrData, pullChi2Data = fit.doFit('datafit',hPullData,x,k,1,'2g')
            resPullData.Draw("same")
            resPullData.SetLineStyle(2)
            
            c1.Update()
            
            finfoData = hPullData.GetListOfFunctions().FindObject("stats")
            finfoData.__class__ = ROOT.TPaveStats
            finfoData.SetX1NDC(0.7)
            finfoData.SetX2NDC(0.95)
            finfoData.SetY1NDC(0.20)
            finfoData.SetY2NDC(0.40)
            lData = ROOT.TText(0.81,0.41,"Data (fit)")
            lData.SetTextSize(0.035)
            lData.SetNDC()
            lData.Draw()
            
            finfoMC = hPullMC.GetListOfFunctions().FindObject("stats")
            finfoMC.__class__ = ROOT.TPaveStats
            finfoMC.SetX1NDC(0.7)
            finfoMC.SetX2NDC(0.95)
            finfoMC.SetY1NDC(0.45)
            finfoMC.SetY2NDC(0.65)
            lMC = ROOT.TText(0.81,0.66,"Simulation (fit)")
            lMC.SetTextSize(0.035)
            lMC.SetNDC()
            lMC.Draw()

            lDataPull = ROOT.TLatex(0.20,0.63,"#sigma^{Data}_{"+x+"} = %.2f" % (pullData))
            lDataPull.SetNDC(); lDataPull.SetTextFont(43); lDataPull.SetTextSize(20); lDataPull.Draw()

            lMCPull = ROOT.TLatex(0.20,0.70,"#sigma^{Sim.}_{"+x+"} = %.2f" % (pullMC))
            lMCPull.SetNDC(); lMCPull.SetTextFont(43); lMCPull.SetTextSize(20); lMCPull.Draw()

            lSel = ROOT.TLatex(0.20,0.80,str(v['bins'][1])+' < '+selName+' < '+str(v['bins'][2])+units)
            lSel.SetTextSize(0.035)
            lSel.SetNDC()
            lSel.Draw()
            
            c1.Update()
            
            leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
            leg.SetFillColor(253)
            leg.SetBorderSize(0)
            leg.AddEntry(hPullData,"Data","p")
            leg.AddEntry(hPullMC,"Simulation","f")
            leg.AddEntry(resPullData,"Data (fit)","l")
            leg.AddEntry(resPullMC,"Simulation (fit)","l")
            leg.Draw()        
            
            t1, t2, t3, t4 = style.cmslabel(1, c.year, evt)
            t1.Draw()
            t2.Draw()
            t3.Draw()
            t4.Draw()
            
            c1.Print(options.output+'/pvPull_'+x+k+'.pdf')
            c1.Clear()
            
            # collect fit results
#            if k != '':
                # append bin index to the output string
            rout['reso']['data'][x][k].update([('value',resoData), ('error',resoErrData)])
            rout['reso']['mc'][x][k].update([('value',resoMC), ('error',resoErrMC)])
            rout['pull']['data'][x][k].update([('value',pullData), ('error',pullErrData)])
            rout['pull']['mc'][x][k].update([('value',pullMC), ('error',pullErrMC)])
            rout['widthx']['data'][x][k].update([('value',beamWidthXData), ('error',0.)])
            rout['widthx']['mc'][x][k].update([('value',beamWidthXMC), ('error',0.)])
            rout['widthy']['data'][x][k].update([('value',beamWidthYData), ('error',0.)])
            rout['widthy']['mc'][x][k].update([('value',beamWidthYMC), ('error',0.)])
            rout['sigmaz']['data'][x][k].update([('value',sigmaZData), ('error',0.)])
            rout['sigmaz']['mc'][x][k].update([('value',sigmaZMC), ('error',0.)])

    if not os.path.isdir('results'): os.system('mkdir results/')
      
    foutput = 'results/pv_zb.json'
    if options.qcd: foutput = 'results/pv_qcd.json'
    with open(foutput, 'w') as write_file:
        json.dump(rout, write_file, indent=2)
