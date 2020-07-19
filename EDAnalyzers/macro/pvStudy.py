import os
import sys
import math
import subprocess
import multiprocessing
import common as c
from subprocess import call
import style
import functions as fun
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
    parser.add_option("--data", default="data.root", help="Input data file name [default: %default]")
    parser.add_option("--mc", default="mc.root", help="Input mc file name [default: %default]")
    parser.add_option("--output", default="pics", help="Output directory [default: %default]")
    parser.add_option("--qcd", action='store_true', help="Use QCD events [default: %default]")
    parser.add_option("--param", default="sumTrackPtSq", help="Parameterisation for PV resolution measurement [default: %default]")
    parser.add_option("--image", default="eps", help="Image format [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

def runFit(evt, v, x, kstr, img, hResoData, hResoMC, hPullData, hPullMC):

    ROOT.gROOT.Reset()
    
    ROOT.gROOT.SetBatch()
    
    pstyle = style.SetPlotStyle(1, 'FIT')
    pstyle.SetOptFit(1111)
    
    figs = []
            
    if hResoData.GetEntries() < 10:
        print 'No stats in data: '+hResoData.GetName()
        sys.exit()
    if hResoMC.GetEntries() < 10:
        print 'No stats in mc: '+hResoMC.GetName()
        sys.exit()
        
    for h in [hResoData,hPullData]:
        h.SetMarkerStyle(c.datamrk)
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
    
    # Rebin in case of fine-bin measurement
#    hResoData = hResoData.Rebin(2)
#    hResoMC = hResoMC.Rebin(2)
#    hPullData = hPullData.Rebin(5)
#    hPullMC = hPullMC.Rebin(5)
            
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
    if options.param == 'sumTrackPt':
        selName = '#sump_{T}'
        units = ' GeV'
    elif options.param == 'sumTrackPtSq': 
        selName = '#sqrt{#sump^{2}_{T}}'
        units = ' GeV'
            
    # Resolution
        
    c1 = ROOT.TCanvas()
        
    fun.adjust(hResoMC, hResoData, nsig=5)
            
    hResoMC.Draw('hist')
    hResoData.Draw('e1 sames')

    resResoMC, resoMC, resoErrMC, resoChi2MC = fit.doFit('mcfit_reso', hResoMC, x, kstr, c.mcfit, '2g', nsig=3, nTries=3)
    resResoMC.Draw("same")
            
    resResoData, resoData, resoErrData, resoChi2Data = fit.doFit('datafit_reso', hResoData, x, kstr, 1, '2g', nsig=3, nTries=3)
    resResoData.Draw("same")
    resResoData.SetLineStyle(2)
    
    sysResoData = hResoData.GetXaxis().GetBinWidth(2)
    sysResoMC = hResoMC.GetXaxis().GetBinWidth(2)
            
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
            
    b = fun.isbadfit(resoChi2MC, resoChi2Data)
    b.Draw()

    foutput = 'pvReso_'+x+kstr
    figs.append(foutput)
    c1.Print(options.output+'/'+foutput+'.'+img)
    c1.Clear()

    # Pull
            
    fun.adjust(hPullMC, hPullData, nsig=5)
        
    hPullMC.Draw('hist')
    hPullData.Draw('e1 sames')

    resPullMC, pullMC, pullErrMC, pullChi2MC = fit.doFit('mcfit_pull', hPullMC, x, kstr, c.mcfit, '1g', nsig=2.5, nTries=3)
    resPullMC.Draw("same")
    
    resPullData, pullData, pullErrData, pullChi2Data = fit.doFit('datafit_pull', hPullData, x, kstr, 1, '1g', nsig=2.5, nTries=3)
    resPullData.Draw("same")
    resPullData.SetLineStyle(2)

    sysPullData = hPullData.GetXaxis().GetBinWidth(2)
    sysPullMC = hPullMC.GetXaxis().GetBinWidth(2)
    
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
    
    b = fun.isbadfit(pullChi2MC, pullChi2Data)
    b.Draw()
    
    foutput = 'pvPull_'+x+kstr
    figs.append(foutput)
    c1.Print(options.output+'/'+foutput+'.'+img)
    c1.Clear()

    return figs, x, kstr, resoData, resoErrData, resoMC, resoErrMC, pullData, pullErrData, pullMC, pullErrMC, sysResoData, sysResoMC, sysPullData, sysPullMC

if __name__ == '__main__':
    
    options = main()
    
    pool = multiprocessing.Pool(c.ncores)

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(1, 'MAIN')

    if not os.path.isdir(options.output):
        os.system("mkdir "+options.output)
        
    fHistData = ROOT.TFile.Open(options.data, 'read')
    fHistMC = ROOT.TFile.Open(options.mc, 'read')
    
    figs = []
    
    evt = 'zb'
    if options.qcd: evt = 'qcd'

    ppath = 'data/bins/'
    fparam = fun.param(ppath+evt+'_pv.json')
#    fparam = fun.param(ppath+evt+'_bs.json')
    param = fparam.get(options.param)
    
    img = options.image
    
    hIncl = ['jetPtMax', 'jetHT', 'evNpv', 'pvNTrks', 'pvSumTrackPt', 'pvSumTrackPt2']

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
        
        foutput = h
        figs.append(foutput)
        c1.Print(options.output+'/'+foutput+'.'+img)
        c1.Clear()
    
    #### Basic distributions

    print 'Plot PV and BS profiles'
    
    #### PV
        
    pv_hXData = fHistData.Get('h_pvx')
    pv_hYData = fHistData.Get('h_pvy')
    pv_hZData = fHistData.Get('h_pvz')
    
    pv_hXMC = fHistMC.Get('h_pvx')
    pv_hYMC = fHistMC.Get('h_pvy')
    pv_hZMC = fHistMC.Get('h_pvz')
    
    pv_hXYData = fHistData.Get('h2_pvx_y')
    pv_hXZData = fHistData.Get('h2_pvx_z')
    pv_hYZData = fHistData.Get('h2_pvy_z')
    
    pv_hXYMC = fHistMC.Get('h2_pvx_y')
    pv_hXZMC = fHistMC.Get('h2_pvx_z')
    pv_hYZMC = fHistMC.Get('h2_pvy_z')

    #### BS
        
    bs_hXData = fHistData.Get('h_bsx0')
    bs_hYData = fHistData.Get('h_bsy0')
    bs_hZData = fHistData.Get('h_bsz0')
    
    bs_hXMC = fHistMC.Get('h_bsx0')
    bs_hYMC = fHistMC.Get('h_bsy0')
    bs_hZMC = fHistMC.Get('h_bsz0')
    
    bs_hXYData = fHistData.Get('h2_bsx0_y0')
    bs_hXZData = fHistData.Get('h2_bsx0_z0')
    bs_hYZData = fHistData.Get('h2_bsy0_z0')
    
    bs_hXYMC = fHistMC.Get('h2_bsx0_y0')
    bs_hXZMC = fHistMC.Get('h2_bsx0_z0')
    bs_hYZMC = fHistMC.Get('h2_bsy0_z0')
    
    bs_hBeamWidthXData = fHistData.Get('h_bsBeamWidthX')
    bs_hBeamWidthYData = fHistData.Get('h_bsBeamWidthY')
    bs_hSigmaZData = fHistData.Get('h_bsSigmaZ')

    bs_hBeamWidthXMC = fHistMC.Get('h_bsBeamWidthX')
    bs_hBeamWidthYMC = fHistMC.Get('h_bsBeamWidthY')
    bs_hSigmaZMC = fHistMC.Get('h_bsSigmaZ')

    for idx, h in enumerate([pv_hXData,pv_hYData,pv_hZData,pv_hXMC,pv_hYMC,pv_hZMC]):
        
        c1 = ROOT.TCanvas()
        
        h.Draw('e1')

        proj = 'X'
        if 'y' in h.GetName(): proj = 'Y'
        if 'z' in h.GetName(): proj = 'Z'
        samp = 'data'
        if idx > 2: samp = 'mc'
            
        foutput = 'pv'+proj+'_'+samp
        figs.append(foutput)
        c1.Print(options.output+'/'+foutput+'.'+img)
        c1.Clear()

    for idx, h in enumerate([pv_hXYData,pv_hXZData,pv_hYZData,pv_hXYMC,pv_hXZMC,pv_hYZMC]):
        
        c1 = ROOT.TCanvas()
        h.Draw('COLZ')
        h.SetMinimum(3)
            
        proj = 'XY'
        if 'x_z' in h.GetName(): proj = 'XZ'
        if 'y_z' in h.GetName(): proj = 'YZ'
        samp = 'data'
        if idx > 2: samp = 'mc'

        if 'x_y' in h.GetName():
            if samp == 'data':
                hpvplot = pv_hXYData
                hbsplot = bs_hXYData
#                hbsplot.Draw('BOX SAME')
#                hbsplot.SetLineColor(ROOT.kRed)
            else:
                hpvplot = pv_hXYMC
                hbsplot = bs_hXYMC
#                hbsplot.Draw('BOX SAME')
#                hbsplot.SetLineColor(ROOT.kRed)
            lpvrms1 = ROOT.TLatex(0.2,0.3,"RMS (PV)"); lpvrms1.SetNDC(); lpvrms1.SetTextFont(43); lpvrms1.SetTextSize(20); lpvrms1.Draw()
            lpvrms2 = ROOT.TLatex(0.2,0.25,"x = %.1f" % (hpvplot.GetRMS(1)*math.pow(10,3))+" #mum"); lpvrms2.SetNDC(); lpvrms2.SetTextFont(43); lpvrms2.SetTextSize(20); lpvrms2.Draw()
            lpvrms3 = ROOT.TLatex(0.2,0.2,"y = %.1f" % (hpvplot.GetRMS(2)*math.pow(10,3))+" #mum"); lpvrms3.SetNDC(); lpvrms3.SetTextFont(43); lpvrms3.SetTextSize(20); lpvrms3.Draw()
#            if samp == 'data':
#                lbsrms1 = ROOT.TLatex(0.65,0.3,"RMS (BS)"); lbsrms1.SetNDC(); lbsrms1.SetTextFont(43); lbsrms1.SetTextSize(20); lbsrms1.SetTextColor(ROOT.kRed); lbsrms1.Draw()
#                lbsrms2 = ROOT.TLatex(0.65,0.25,"x = %.1f" % (hbsplot.GetRMS(1)*math.pow(10,3))+" #mum"); lbsrms2.SetNDC(); lbsrms2.SetTextFont(43); lbsrms2.SetTextSize(20); lbsrms2.SetTextColor(ROOT.kRed); lbsrms2.Draw()
#                lbsrms3 = ROOT.TLatex(0.65,0.2,"y = %.1f" % (hbsplot.GetRMS(2)*math.pow(10,3))+" #mum"); lbsrms3.SetNDC(); lbsrms3.SetTextFont(43); lbsrms3.SetTextSize(20); lbsrms3.SetTextColor(ROOT.kRed); lbsrms3.Draw()
                    
        if 'x_z' in h.GetName():
            if samp == 'data':
                hpvplot = pv_hXZData
                hbsplot = bs_hXZData
#                hbsplot.Draw('BOX SAME')
#                hbsplot.SetLineColor(ROOT.kRed)
            else:
                hpvplot = pv_hXZMC
                hbsplot = bs_hXZMC
#                hbsplot.Draw('BOX SAME')
#                hbsplot.SetLineColor(ROOT.kRed)
            lpvrms1 = ROOT.TLatex(0.2,0.3,"RMS (PV)"); lpvrms1.SetNDC(); lpvrms1.SetTextFont(43); lpvrms1.SetTextSize(20); lpvrms1.Draw()
            lpvrms2 = ROOT.TLatex(0.2,0.25,"x = %.1f" % (hpvplot.GetRMS(2)*math.pow(10,3))+" #mum"); lpvrms2.SetNDC(); lpvrms2.SetTextFont(43); lpvrms2.SetTextSize(20); lpvrms2.Draw()
            lpvrms3 = ROOT.TLatex(0.2,0.2,"z = %.1f" % (hpvplot.GetRMS(1)*10.)+" mm"); lpvrms3.SetNDC(); lpvrms3.SetTextFont(43); lpvrms3.SetTextSize(20); lpvrms3.Draw()
#            if samp == 'data':
#                lbsrms1 = ROOT.TLatex(0.65,0.3,"RMS (BS)"); lbsrms1.SetNDC(); lbsrms1.SetTextFont(43); lbsrms1.SetTextSize(20); lbsrms1.SetTextColor(ROOT.kRed); lbsrms1.Draw()
#                lbsrms2 = ROOT.TLatex(0.65,0.25,"x = %.1f" % (hbsplot.GetRMS(2)*math.pow(10,3))+" #mum"); lbsrms2.SetNDC(); lbsrms2.SetTextFont(43); lbsrms2.SetTextSize(20); lbsrms2.SetTextColor(ROOT.kRed); lbsrms2.Draw()
#                lbsrms3 = ROOT.TLatex(0.65,0.2,"z = %.1f" % (hbsplot.GetRMS(1))+" mm"); lbsrms3.SetNDC(); lbsrms3.SetTextFont(43); lbsrms3.SetTextSize(20); lbsrms3.SetTextColor(ROOT.kRed); lbsrms3.Draw()
                    
        if 'y_z' in h.GetName():
            if samp == 'data':
                hpvplot = pv_hYZData
                hbsplot = bs_hYZData
#                hbsplot.Draw('BOX SAME')
#                hbsplot.SetLineColor(ROOT.kRed)
            else:
                hpvplot = pv_hYZMC
                hbsplot = bs_hYZMC
#                hbsplot.Draw('BOX SAME')
#                hbsplot.SetLineColor(ROOT.kRed)
            lpvrms1 = ROOT.TLatex(0.2,0.3,"RMS (PV)"); lpvrms1.SetNDC(); lpvrms1.SetTextFont(43); lpvrms1.SetTextSize(20); lpvrms1.Draw()
            lpvrms2 = ROOT.TLatex(0.2,0.25,"y = %.1f" % (hpvplot.GetRMS(2)*math.pow(10,3))+" #mum"); lpvrms2.SetNDC(); lpvrms2.SetTextFont(43); lpvrms2.SetTextSize(20); lpvrms2.Draw()
            lpvrms3 = ROOT.TLatex(0.2,0.2,"z = %.1f" % (hpvplot.GetRMS(1)*10.)+" mm"); lpvrms3.SetNDC(); lpvrms3.SetTextFont(43); lpvrms3.SetTextSize(20); lpvrms3.Draw()
#            if samp == 'data':
#                lbsrms1 = ROOT.TLatex(0.65,0.3,"RMS (BS)"); lbsrms1.SetNDC(); lbsrms1.SetTextFont(43); lbsrms1.SetTextSize(20); lbsrms1.SetTextColor(ROOT.kRed); lbsrms1.Draw()
#                lbsrms2 = ROOT.TLatex(0.65,0.25,"y = %.1f" % (hbsplot.GetRMS(2)*math.pow(10,3))+" #mum"); lbsrms2.SetNDC(); lbsrms2.SetTextFont(43); lbsrms2.SetTextSize(20); lbsrms2.SetTextColor(ROOT.kRed); lbsrms2.Draw()
#                lbsrms3 = ROOT.TLatex(0.65,0.2,"z = %.1f" % (hbsplot.GetRMS(1))+" mm"); lbsrms3.SetNDC(); lbsrms3.SetTextFont(43); lbsrms3.SetTextSize(20); lbsrms3.SetTextColor(ROOT.kRed); lbsrms3.Draw()

        t1, t2, t3, t4 = style.cmslabel(1, c.year, evt)
        t1.Draw()
        t2.Draw()
        t3.Draw()
        t4.Draw()

        foutput = 'pv'+proj+'_'+samp
        figs.append(foutput)
        for outdir in [options.output, 'pub']:
            c1.Print(outdir+'/'+foutput+'.'+img)
        c1.Print('pub_dps/'+foutput+'.root')
        c1.Clear()

    for idx, h in enumerate([bs_hBeamWidthXData,bs_hBeamWidthYData,bs_hSigmaZData,bs_hBeamWidthXMC,bs_hBeamWidthYMC,bs_hSigmaZMC]):
        
        c1 = ROOT.TCanvas()
        
        h.Draw('hist')
        
        proj = 'X'
        if 'Y' in h.GetName(): proj = 'Y'
        if 'Z' in h.GetName(): proj = 'Z'
        samp = 'data'
        if idx > 2: samp = 'mc'
        
        t1, t2, t3, t4 = style.cmslabel(1, c.year, evt)
        t1.Draw()
        t2.Draw()
        t3.Draw()
        t4.Draw()
        
        foutput = 'bs'+proj+'_'+samp
        figs.append(foutput)
        c1.Print(options.output+'/'+foutput+'.'+img)
        c1.Clear()
        
    beamWidthXData = bs_hBeamWidthXData.GetMean()
    beamWidthYData = bs_hBeamWidthYData.GetMean()
    sigmaZData = bs_hSigmaZData.GetMean()

    beamWidthXMC = bs_hBeamWidthXMC.GetMean()
    beamWidthYMC = bs_hBeamWidthYMC.GetMean()
    sigmaZMC = bs_hSigmaZMC.GetMean()
    
    pstyle.Reset()
    
    #### Fits
    
    print 'Run PV fits'
    
    meas = ['reso','pull','widthx','widthy','sigmaz']
    rout = {}
    for l1 in meas: rout[l1] = {}
    for t in ['data','mc']:
        for l2 in meas: rout[l2][t] = {}
        for x in c.PVmeas:
            for l3 in meas: rout[l3][t][x] = {}
            for k, v in param.iteritems():
                if k in ['allbins']: continue
                for l4 in meas: rout[l4][t][x][k] = {}
            
    jobs = []
    
    for k, v in param.iteritems():
        
        if k in ['allbins', '']: continue
#        if k not in ['_sumTrackPtSq8p999to9p072']: continue
    
        for x in c.PVmeas:

            kstr = str(k)

            hNameResoData = 'h_pvd'+x+'12'+kstr
            hNamePullData = 'h_pvd'+x+'Pull12'+kstr
            hResoData = fHistData.Get(hNameResoData).Clone('hResoData')
            hPullData = fHistData.Get(hNamePullData).Clone('hPullData')
            fun.addbin(hResoData)
            fun.addbin(hPullData)
        
            hNameResoMC = 'h_pvd'+x+'12'+kstr
            hNamePullMC = 'h_pvd'+x+'Pull12'+kstr
            hResoMC = fHistMC.Get(hNameResoMC).Clone('hNameResoMC')
            hPullMC = fHistMC.Get(hNamePullMC).Clone('hNamePullMC')
            fun.addbin(hResoMC)
            fun.addbin(hPullMC)
            
            jobs.append( pool.apply_async(runFit, (evt, v, x, kstr, img, hResoData, hResoMC, hPullData, hPullMC)) )
            
            rout['widthx']['data'][x][k].update([('value',beamWidthXData), ('error',0.)])
            rout['widthx']['mc'][x][k].update([('value',beamWidthXMC), ('error',0.)])
            rout['widthy']['data'][x][k].update([('value',beamWidthYData), ('error',0.)])
            rout['widthy']['mc'][x][k].update([('value',beamWidthYMC), ('error',0.)])
            rout['sigmaz']['data'][x][k].update([('value',sigmaZData), ('error',0.)])
            rout['sigmaz']['mc'][x][k].update([('value',sigmaZMC), ('error',0.)])

    for job in jobs: 
        
        result = job.get()
        
        figs += result[0]
        
        x = result[1]
        kstr = result[2]
        
        resoData = result[3]
        resoErrData = result[4]
        resoMC = result[5]
        resoErrMC = result[6]
        pullData = result[7]
        pullErrData = result[8]
        pullMC = result[9]
        pullErrMC = result[10]
        resoSysData = result[11]
        resoSysMC = result[12]
        pullSysData = result[13]
        pullSysMC = result[14]
        
        rout['reso']['data'][x][kstr].update([('value',resoData), ('error',resoErrData), ('sys',resoSysData)])
        rout['reso']['mc'][x][kstr].update([('value',resoMC), ('error',resoErrMC), ('sys',resoSysMC)])
        rout['pull']['data'][x][kstr].update([('value',pullData), ('error',pullErrData), ('sys',pullSysData)])
        rout['pull']['mc'][x][kstr].update([('value',pullMC), ('error',pullErrMC), ('sys',pullSysMC)])
        
    pool.close()

    if not os.path.isdir('results'): os.system('mkdir results/')

    print 'Store results'
    foutput = 'results/pv_zb.json'
    if options.qcd: foutput = 'results/pv_qcd.json'
    with open(foutput, 'w') as write_file:
        json.dump(rout, write_file, indent=2)
        
    print 'Generate thumbnails'
    fun.generateThumbs(figs, options.output, img)
    print 'Put plots on a web page'
    fun.createPage('pv')
