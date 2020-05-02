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
    parser.add_option("-d","--data",default="data.root",help="input data file name [default: %default]")
    parser.add_option("-m","--mc",default="mc.root",help="input mc file name [default: %default]")
    parser.add_option("-o","--output",default="pics",help="output directory [default: %default]")
    parser.add_option("-t","--type",default="pv",help="type of the measurement [default: %default]")
    parser.add_option("--qcd",action='store_true',help="Use QCD events [default: %default]")
    parser.add_option("--parampv",default="PVsumTrackPtSq",help="Parameterisation for PV resolution measurement [default: %default]")
#    parser.add_option("-s","--switch",default="",help="use only mc or data [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(1)

    if not os.path.isdir(options.output):
        os.system("mkdir "+options.output)

    ip = options.type

    evt = 'zb'
    if options.qcd: evt = 'qcd'

    parampv = c.PVnTracks[evt]
    if options.parampv == 'PVsumTrackPt': parampv = c.PVsumTrackPt[evt]
    elif options.parampv == 'PVsumTrackPtSq': parampv = c.PVsumTrackPtSq[evt]
    elif options.parampv != 'PVnTracks':
        print 'Unknown PV parameterisation requested:', options.parampv
        sys.exit()        
    
    fHistData = ROOT.TFile.Open(options.data,'read')
    fHistMC = ROOT.TFile.Open(options.mc,'read')
    
#    ParamList = [c.IPpt,c.IPeta,c.IPphi,c.IPnpv,c.IPsz]
#    ParamList = [c.IPpt,c.IPeta,c.IPphi,c.IPnpv,c.IPdr,c.IPsz]
    ParamList = [c.IPpt]

    hIncl = [\
    'evNpv','ipPt','ipEta','ipPhi',\
    'ipDrTrkJet','ipNTrkJet',\
    'ipD0','ipDz','ipSD0','ipSDz',\
    'ippvD0','ippvDz','ippvSD0','ippvSDz',\
    'ipbsD0','ipbsDz','ipbsSD0','ipbsSDz',\
    'ipbszpvD0','ipbszpvSD0',\
    'ipbszpcaD0','ipbszpcaSD0','bsSigmaZ'\
    ]
    for h in hIncl:
        
        c1 = ROOT.TCanvas()
        
        hData = fHistData.Get('h_'+h)
        hMC = fHistMC.Get('h_'+h)
        
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

        if h in ['ipPt','ipDrTrkJet']:
            hData.SetMinimum(10)
            hMC.SetMinimum(10)
            c1.SetLogy(1)
        else: c1.SetLogy(0)
        
        c1.Print(options.output+'/'+h+'.pdf')
        c1.Clear()
      
    #### Fits
    
    pstyle.SetOptFit(1111)
    
    rout = {}
    rout['reso'] = {}
    for t in ['data','mc']:
        rout['reso'][t] = {}
        for x in c.IPmeas:
            rout['reso'][t][x] = {}
            for p in ParamList:
                for k, v in p.iteritems():
                    rout['reso'][t][x][k] = {}
                    for ktrk, vtrk in parampv.iteritems():
                        if ip == 'bs' and ktrk != '': break
                        rout['reso'][t][x][k][ktrk] = {}

    for ktrk, vtrk in parampv.iteritems():

        if ip == 'pv' and ktrk == '': continue
        if ip == 'bs' and ktrk != '': break
        
#        if ktrk != '_nTrks0to15': continue
        
        for param in ParamList:
        
            for k, v in param.iteritems():
    
                if k == '': continue
                
#                if k != '_pt1p3to1p5': continue
                
                for x in c.IPmeas:                    

                    hNameResoData = 'h_ip'+ip+x+ktrk+k
                    hResoData = fHistData.Get(hNameResoData)
#                    func.addbin(hResoData)
        
                    hNameResoMC = 'h_ip'+ip+x+ktrk+k
                    hResoMC = fHistMC.Get(hNameResoMC)
#                    func.addbin(hResoMC)
                    
                    hasStats = True
                    if hResoData.Integral() == 0:
                        hasStats = False
                        print 'No stats in '+hResoData.GetName()+' (data)'                        
                    if hResoMC.Integral() == 0:
                        hasStats = False
                        print 'No stats in '+hResoMC.GetName()+' (mc)'

                    for h in [hResoData]:
                        h.SetMarkerSize(0.7)
                        h.SetMarkerColor(1)
                        h.SetLineColor(1)

                    for h in [hResoMC]:
                        h.SetMarkerSize(0)
                        h.SetMarkerColor(c.mccol)
                        h.SetLineColor(c.mccol)
                        h.SetFillColor(c.mccol)
                        h.SetLineStyle(1)
        
                    intResoMC = hResoMC.Integral()
                    intResoData = hResoData.Integral()
                    
                    if intResoData < 20000 or intResoMC < 20000:
                        hResoMC = hResoMC.Rebin(2)
                        hResoData = hResoData.Rebin(2)
#                        if intResoData < 1000 or intResoMC < 1000:
#                            hResoMC = hResoMC.Rebin(5)
#                            hResoData = hResoData.Rebin(5)
                    
                    if hasStats:
                        hResoMC.Scale(intResoData/intResoMC)        
            
                    maxResoData = hResoData.GetMaximum()
                    maxResoMC = hResoMC.GetMaximum()
                    hResoMC.SetMaximum(1.2*max(maxResoData,maxResoMC))
                    hResoMC.SetMinimum(0.)
        
                    c1 = ROOT.TCanvas()
                    
                    hResoMC.Draw('hist')
                    hResoData.Draw('e1 sames')
                    
                    if intResoData > 100 and intResoMC > 100:

                        ffit = ''
                        if x in ['dz']: 
                            ffit = '3g'
#                            if lowstat: ffit = '1g'
                        
                        resResoMC, resoMC, resoErrMC, resoChi2MC = fit.doFit('mcfit',hResoMC,x,k,c.mcfit,ffit)

                        resResoData, resoData, resoErrData, resoChi2Data = fit.doFit('datafit',hResoData,x,k,1,ffit)

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

                        xLabel = x
                        if x == 'd0': xLabel = 'd_{xy}'
                        elif x == 'dz': xLabel = 'd_{z}'
                        
                        lDataReso = ROOT.TLatex(0.20,0.56,"#sigma^{Data}_{"+xLabel+"} = %.1f #mum" % (resoData))
                        lDataReso.SetNDC(); lDataReso.SetTextFont(43); lDataReso.SetTextSize(20); lDataReso.Draw()
                        
                        lMCReso = ROOT.TLatex(0.20,0.65,"#sigma^{Sim.}_{"+xLabel+"} = %.1f #mum" % (resoMC))
                        lMCReso.SetNDC(); lMCReso.SetTextFont(43); lMCReso.SetTextSize(20); lMCReso.Draw()

                        if ip != 'bs':
                            lSelPV = ROOT.TLatex(0.20,0.80,str(vtrk['bins'][1])+' < N_{trk} < '+str(vtrk['bins'][2]))
                            lSelPV.SetTextSize(0.035)
                            lSelPV.SetNDC()
                            lSelPV.Draw()

                        pLabel = 'p_{T}'
                        pUnits = 'GeV'
                        pPrec = '%.1f'
                        if param == c.IPeta: 
                            pLabel = '#eta'
                            pUnits = ''
                        elif param == c.IPphi: 
                            pLabel = '#phi'
                            pUnits = ''
                        elif param == c.IPnpv: 
                            pLabel = 'N_{PV}'
                            pUnits = ''
                            pPrec = '%d'
                        elif param == c.IPdr: 
                            pLabel = '#DeltaR'
                            pUnits = ''
                        elif param == c.IPsz: 
                            pLabel = '#sigma_{z}'
                            pUnits = '#mum'

                        lSel = ROOT.TLatex(0.20,0.75,(pPrec+' < '+pLabel+' < '+pPrec+' '+pUnits) % (v['bins'][1], v['bins'][2]))
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
                
                        c1.Print(options.output+'/ip'+ip+'Reso_'+x+ktrk+k+'.pdf')
                        c1.Clear()
                        
                    else:
                        
                        resoMC, resoErrMC, resoData, resoErrData = 0, 0, 0, 0
                
                    # collect fit results
                    if k != '':
                        rout['reso']['data'][x][k][ktrk].update([('value',resoData), ('error',resoErrData), ('int',intResoData)])
                        rout['reso']['mc'][x][k][ktrk].update([('value',resoMC), ('error',resoErrMC), ('int',intResoMC)])

    if not os.path.isdir('results'): os.system('mkdir results/')
                        
    foutput = 'results/ip_zb.json'
    if options.qcd: foutput = 'results/ip_qcd.json'
    with open(foutput, "w") as write_file:
        json.dump(rout, write_file, indent=2)
