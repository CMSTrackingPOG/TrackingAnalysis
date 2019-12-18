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

    usage = "usage: %prog [options]\n Analysis script to study PV resolution"
    
    parser = OptionParser(usage)
    parser.add_option("-i","--input",default="output.root",help="input file name [default: %default]")
    parser.add_option("-o","--output",default="pics",help="output directory [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(1)

    if not os.path.isdir(options.output):
        os.system("mkdir "+options.output)

    fHist = ROOT.TFile.Open(options.input,'read')

    hIncl = ['pvNTrks']
    for h in hIncl:
        
        hData = fHist.Get('h_'+h+'__data')
        hMC = fHist.Get('h_'+h+'__mc')
        
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
        hMC.SetMarkerColor(ROOT.kBlue-10)
        hMC.SetLineColor(ROOT.kBlue-10)
        hMC.SetFillColor(ROOT.kBlue-10)
        hMC.SetLineStyle(1)
        
        leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
        leg.SetFillColor(253)
        leg.SetBorderSize(0)
        leg.AddEntry(hData,"Data","p")
        leg.AddEntry(hMC,"Simulation","f")
        leg.Draw()        

        t1, t2 = style.cmslabel(1,777)
        t1.Draw()
        
        c1.Print(options.output+'/'+h+'.eps')
        c1.Clear()
    
    #### Basic distributions
    
    hXData = {}; hYData = {}; hZData = {}
    hXMC = {}; hYMC = {}; hZMC = {}
    hXYData = {}; hXZData = {}; hYZData = {}
    hXYMC = {}; hXZMC = {}; hYZMC = {}
    
    for k, v in c.PVnTracks.iteritems():
        
        hXData[k] = fHist.Get('h_pvx'+k+'__data')
        hYData[k] = fHist.Get('h_pvy'+k+'__data')
        hZData[k] = fHist.Get('h_pvz'+k+'__data')

        hXMC[k] = fHist.Get('h_pvx'+k+'__mc')
        hYMC[k] = fHist.Get('h_pvy'+k+'__mc')
        hZMC[k] = fHist.Get('h_pvz'+k+'__mc')
        
        hXYData[k] = fHist.Get('h2_pvx_y'+k+'__data')
        hXZData[k] = fHist.Get('h2_pvx_z'+k+'__data')
        hYZData[k] = fHist.Get('h2_pvy_z'+k+'__data')
        
        hXYMC[k] = fHist.Get('h2_pvx_y'+k+'__mc')
        hXZMC[k] = fHist.Get('h2_pvx_z'+k+'__mc')
        hYZMC[k] = fHist.Get('h2_pvy_z'+k+'__mc')

    for h in [hXData,hYData,hZData,hXMC,hYMC,hZMC]:
        for k, hh in h.iteritems():
        
            c1 = ROOT.TCanvas()

            hh.Draw('e1')

            proj = 'X'
            if 'y' in hh.GetName(): proj = 'Y'
            if 'z' in hh.GetName(): proj = 'Z'
            samp = 'data'
            if 'mc' in hh.GetName(): samp = 'mc'
            
            c1.Print(options.output+'/pv'+proj+k+'_'+samp+'.eps')
            c1.Clear()

    for h2 in [hXYData,hXZData,hYZData,hXYMC,hXZMC,hYZMC]:
        for k, hh in h2.iteritems():
        
            c1 = ROOT.TCanvas()
            hh.Draw('COLZ')
            
            proj = 'XY'
            if 'x_z' in hh.GetName(): proj = 'XZ'
            if 'y_z' in hh.GetName(): proj = 'YZ'
            samp = 'data'
            if 'mc' in hh.GetName(): samp = 'mc'
            
            c1.Print(options.output+'/pv'+proj+k+'_'+samp+'.eps')
            c1.Clear()
        
    #### Fits
    
    pstyle.SetOptFit(1111)
    
    rout = {}
    rout['reso'] = {}
    rout['pull'] = {}
    for t in ['data','mc']:
        rout['reso'][t] = {}
        rout['pull'][t] = {}
        for x in c.PVmeas:
            rout['reso'][t][x] = {}
            rout['pull'][t][x] = {}
            for k, v in c.PVnTracks.iteritems():
                rout['reso'][t][x][k] = {}
                rout['pull'][t][x][k] = {}

    for k, v in c.PVnTracks.iteritems():
    
        for x in c.PVmeas:
            
            hNameResoData = 'h_pvd'+x+'12'+k+'__data'
            hNamePullData = 'h_pvd'+x+'Pull12'+k+'__data'
            hResoData = fHist.Get(hNameResoData)
            hPullData = fHist.Get(hNamePullData)
            func.addbin(hResoData)
            func.addbin(hPullData)
        
            hNameResoMC = 'h_pvd'+x+'12'+k+'__mc'
            hNamePullMC = 'h_pvd'+x+'Pull12'+k+'__mc'
            hResoMC = fHist.Get(hNameResoMC)
            hPullMC = fHist.Get(hNamePullMC)
            func.addbin(hResoMC)
            func.addbin(hPullMC)
            
            if hResoData.GetEntries() < 10:
                print 'No stats in '+hResoData.GetName()
                sys.exit()
            if hResoMC.GetEntries() < 10:
                print 'No stats in '+hResoMC.GetName()
                sys.exit()

            for h in [hResoData,hPullData]:
                h.SetMarkerSize(0.7)
                h.SetMarkerColor(1)
                h.SetLineColor(1)

            for h in [hResoMC,hPullMC]:
                h.SetMarkerSize(0)
                h.SetMarkerColor(ROOT.kBlue-10)
                h.SetLineColor(ROOT.kBlue-10)
                h.SetFillColor(ROOT.kBlue-10)
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

            # Resolution
        
            c1 = ROOT.TCanvas()
        
            hResoMC.Draw('hist')
            hResoData.Draw('e1 sames')
        
            resResoMC, resoMC, resoErrMC, resoChi2MC = fit.doFit('mcfit',hResoMC,'g1',38)
            resResoMC.Draw("same")
            
            resResoData, resoData, resoErrData, resoChi2Data = fit.doFit('datafit',hResoData,'g1',1)
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
            
            c1.Update()

            leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
            leg.SetFillColor(253)
            leg.SetBorderSize(0)
            leg.AddEntry(hResoData,"Data","p")
            leg.AddEntry(hResoMC,"Simulation","f")
            leg.AddEntry(resResoData,"Data (fit)","l")
            leg.AddEntry(resResoMC,"Simulation (fit)","l")
            leg.Draw()        
            
            t1, t2 = style.cmslabel(1,777)
            t1.Draw()
            
            c1.Print(options.output+'/pvReso_'+x+k+'.eps')
            c1.Clear()

            # Pull
        
            hPullMC.Draw('hist')
            hPullData.Draw('e1 sames')
            
            resPullMC, pullMC, pullErrMC, pullChi2MC = fit.doFit('mcfit',hPullMC,'g1',38)
            resPullMC.Draw("same")
            
            resPullData, pullData, pullErrData, pullChi2Data = fit.doFit('datafit',hPullData,'g1',1)
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
            
            c1.Update()
            
            leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
            leg.SetFillColor(253)
            leg.SetBorderSize(0)
            leg.AddEntry(hPullData,"Data","p")
            leg.AddEntry(hPullMC,"Simulation","f")
            leg.AddEntry(resPullData,"Data (fit)","l")
            leg.AddEntry(resPullMC,"Simulation (fit)","l")
            leg.Draw()        
            
            t1, t2 = style.cmslabel(1,777)
            t1.Draw()
            
            c1.Print(options.output+'/pvPull_'+x+k+'.eps')
            c1.Clear()
            
            # collect fit results
            if k != '':
                # append bin index to the output string
                rout['reso']['data'][x][k].update([('value',resoData), ('error',resoErrData)])
                rout['reso']['mc'][x][k].update([('value',resoMC), ('error',resoErrMC)])
                rout['pull']['data'][x][k].update([('value',pullData), ('error',pullErrData)])
                rout['pull']['mc'][x][k].update([('value',pullMC), ('error',pullErrMC)])
            
    with open("pv.json", "w") as write_file:
        json.dump(rout, write_file)
