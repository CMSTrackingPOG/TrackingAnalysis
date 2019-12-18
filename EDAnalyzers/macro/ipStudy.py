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

    usage = "usage: %prog [options]\n Analysis script to study IP resolution"
    
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

    hIncl = ['ipPt','ipEta','ipPhi']
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
      
    #### Fits
    
    pstyle.SetOptFit(1111)
    
    rout = {}
    rout['reso'] = {}
    for t in ['data','mc']:
        rout['reso'][t] = {}
        for x in c.IPmeas:
            rout['reso'][t][x] = {}
            for k, v in c.IPpt.iteritems():
                rout['reso'][t][x][k] = {}
                for ktrk, vtrk in c.PVnTracks.iteritems():
                    rout['reso'][t][x][k][ktrk] = {}

    for ktrk, vtrk in c.PVnTracks.iteritems():
        
        for k, v in c.IPpt.iteritems():
    
            for x in c.IPmeas:

                hNameResoData = 'h_ip'+x+ktrk+k+'__data'
                hResoData = fHist.Get(hNameResoData)
                func.addbin(hResoData)
        
                hNameResoMC = 'h_ip'+x+ktrk+k+'__mc'
                hResoMC = fHist.Get(hNameResoMC)
                func.addbin(hResoMC)
            
                if hResoData.GetEntries() < 1:
                    print 'No stats in '+hResoData.GetName()
                    sys.exit()
                if hResoMC.GetEntries() < 1:
                    print 'No stats in '+hResoMC.GetName()
                    sys.exit()

                for h in [hResoData]:
                    h.SetMarkerSize(0.7)
                    h.SetMarkerColor(1)
                    h.SetLineColor(1)

                for h in [hResoMC]:
                    h.SetMarkerSize(0)
                    h.SetMarkerColor(ROOT.kBlue-10)
                    h.SetLineColor(ROOT.kBlue-10)
                    h.SetFillColor(ROOT.kBlue-10)
                    h.SetLineStyle(1)
        
                intResoMC = hResoMC.Integral()
                intResoData = hResoData.Integral()
                hResoMC.Scale(intResoData/intResoMC)        
            
                maxResoData = hResoData.GetMaximum()
                maxResoMC = hResoMC.GetMaximum()
                hResoMC.SetMaximum(1.2*max(maxResoData,maxResoMC))
                hResoMC.SetMinimum(0.)
        
                c1 = ROOT.TCanvas()
                
                hResoMC.Draw('hist')
                hResoData.Draw('e1 sames')
        
                resResoMC, resoMC, resoErrMC, resoChi2MC = fit.doFit('mcfit',hResoMC,'g1',38)
                resResoData, resoData, resoErrData, resoChi2Data = fit.doFit('datafit',hResoData,'g1',1)
                
#                if (resoChi2MC > 2. or resoChi2Data > 2.) and hResoData.Integral() > 1000:
#                    resResoMC, resoMC, resoErrMC, resoChi2MC = fit.doFit('mcfit',hResoMC,'g2',38)
#                    resResoData, resoData, resoErrData, resoChi2Data = fit.doFit('datafit',hResoData,'g2',1)
                
                resResoMC.Draw("same")
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
                
                c1.Print(options.output+'/ipReso_'+x+ktrk+k+'.eps')
                c1.Clear()
                
                # collect fit results
                if k != '':
                    rout['reso']['data'][x][k][ktrk].update([('value',resoData), ('error',resoErrData), ('int',intResoData)])
                    rout['reso']['mc'][x][k][ktrk].update([('value',resoMC), ('error',resoErrMC), ('int',intResoMC)])
                    
    with open("ip.json", "w") as write_file:
        json.dump(rout, write_file)
