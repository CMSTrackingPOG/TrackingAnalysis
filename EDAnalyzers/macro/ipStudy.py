import os
import sys
import math
import subprocess
import common as c
import multiprocessing
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

    usage = "usage: %prog [options]\n Analysis script to study IP resolution"
    
    parser = OptionParser(usage)
    parser.add_option("--data", default="data.root", help="input data file name [default: %default]")
    parser.add_option("--mc", default="mc.root", help="input mc file name [default: %default]")
    parser.add_option("--output", default="pics", help="output directory [default: %default]")
    parser.add_option("--type", default="pv", help="type of the measurement [default: %default]")
    parser.add_option("--qcd", action='store_true', help="Use QCD events [default: %default]")
    parser.add_option("--parampv", default="sumTrackPtSq", help="Parameterisation for PV resolution measurement [default: %default]")
    parser.add_option("--image", default="eps", help="Image format [default: %default]")
    parser.add_option("--fit", action='store_true', help="Calculate FWHM from fit function [default: %default]")
    parser.add_option("--method", default="fwhm", help="Method to extract resolution (fwhm or fit) [default: %default]")
    parser.add_option("--selection", default="", help="Selection bin for IP resolution [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

def runFit(evt, ip, vtrk, v, x, ktrkstr, kstr, param, img, hResoData, hResoMC):

    ROOT.gROOT.Reset()
    
    ROOT.gROOT.SetBatch()
    
    pstyle = style.SetPlotStyle(1, 'FIT')
    pstyle.SetOptFit(1111)
    
    figs = []

    if hResoData.GetEntries() < 10:
#        print 'No stats in data: '+hResoData.GetName(), ktrkstr
        return figs, x, ktrkstr, kstr, -1, -1, -1, -1, -1, -1, -1, -1

    if hResoMC.GetEntries() < 10:
#        print 'No stats in mc: '+hResoMC.GetName(), ktrkstr
        return figs, x, ktrkstr, kstr, -1, -1, -1, -1, -1, -1, -1, -1
         
    if param in ['pt', 'phi']:
        if options.selection == '':
            rResoData = fit.rebin(hResoData, 60, 30)
            rResoMC = fit.rebin(hResoMC, 60, 30)
        elif kstr == '_pt7p36to7p6':
            rResoData = fit.rebin(hResoData, 80, 30)
            rResoMC = fit.rebin(hResoMC, 80, 30)
        else:
            rResoData = fit.rebin(hResoData, 50, 30)
            rResoMC = fit.rebin(hResoMC, 50, 30)
        rmax = max(rResoData, rResoMC)
        hResoData = hResoData.Rebin(rmax)
        hResoMC = hResoMC.Rebin(rmax)
    else:
        if options.selection != '_pt3p0to10p0' and options.type == 'pv':
            rResoData = fit.rebin(hResoData, 1000, 30)
            rResoMC = fit.rebin(hResoMC, 1000, 30)
            rmax = max(rResoData, rResoMC)
            hResoData = hResoData.Rebin(rmax)
            hResoMC = hResoMC.Rebin(rmax)
        elif options.type == 'bs' and options.qcd:
            rResoData = fit.rebin(hResoData, 500, 30)
            rResoMC = fit.rebin(hResoMC, 500, 30)
            rmax = max(rResoData, rResoMC)
            hResoData = hResoData.Rebin(rmax)
            hResoMC = hResoMC.Rebin(rmax)
        elif options.type == 'bs' and options.selection == '_pt0p0to1p0':
            rResoData = fit.rebin(hResoData, 1000, 30)
            rResoMC = fit.rebin(hResoMC, 1000, 30)
            rmax = max(rResoData, rResoMC)
            hResoData = hResoData.Rebin(rmax)
            hResoMC = hResoMC.Rebin(rmax)        
        else:
            rResoData = fit.rebin(hResoData, 500, 30)
            rResoMC = fit.rebin(hResoMC, 500, 30)
            rmax = max(rResoData, rResoMC)
            hResoData = hResoData.Rebin(rmax)
            hResoMC = hResoMC.Rebin(rmax)        
            
    for h in [hResoData]:
        h.SetMarkerStyle(c.datamrk)
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
    
    hResoMC.Scale(intResoData/intResoMC)
    
    maxResoData = hResoData.GetMaximum()
    maxResoMC = hResoMC.GetMaximum()
    hResoMC.SetMaximum(1.2*max(maxResoData,maxResoMC))
    hResoMC.SetMinimum(0.)

    selName = 'N_{trk}'
    units = ''
    if options.parampv == 'sumTrackPt':
        selName = '#sump_{T}'
        units = ' GeV'
    elif options.parampv == 'sumTrackPtSq': 
        selName = '#sqrt{#sump^{2}_{T}}'
        units = ' GeV'
            
    # Resolution
        
    c1 = ROOT.TCanvas()

    hResoMC.Draw('hist')
    hResoData.Draw('e1 sames')
    
    fun.adjust(hResoMC, hResoData, nsig=2.0)
    
    if options.fit:

        if options.method == 'fwhm': 
            nsig = 2.0
        else:
            ffit = '3g'
            nsig = 2.0

        if x == 'dz': ffit = '3g'
        
        resoChi2MC = 1e+10
        resResoMC, resoMC, resoErrMC, resoChi2MC = fit.doFit('mcfit', hResoMC, x, kstr, c.mcfit, ffit, nsig=nsig, nTries=3)

        resoChi2Data = 1e+10
        resResoData, resoData, resoErrData, resoChi2Data = fit.doFit('datafit', hResoData, x, kstr, 1, ffit, nsig=nsig, nTries=3)

        sysErrData, sysErrMC = hResoData.GetXaxis().GetBinWidth(2), hResoMC.GetXaxis().GetBinWidth(2)

        resResoMC.Draw("same")
        resResoData.Draw("same")
        
        if options.method == 'fwhm':
            # get resolution estimation from bins using the results of the maximum fit
            resoData, resoErrData, sysErrData = fit.fwhm(hResoData, resResoData, nmin=10000)
            resoMC, resoErrMC, sysErrMC = fit.fwhm(hResoMC, resResoMC, nmin=10000)
        
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
        
    else:
        
        resoData, resoErrData, sysErrData = fit.fwhm(hResoData, nmin=10000)
        resoMC, resoErrMC, sysErrMC = fit.fwhm(hResoMC, nmin=10000)

    if resoData*resoMC == 0.:
        resoData, resoErrData, sysErrData = 0., 0., 0.
        resoMC, resoErrMC, sysErrMC = 0., 0., 0.
        
    xLabel = x
    if x == 'd0': xLabel = 'd_{xy}'
    elif x == 'dz': xLabel = 'd_{z}'
    
    lDataReso = ROOT.TLatex(0.20,0.61,"#sigma^{Data}_{"+xLabel+"} = %.1f #mum" % (resoData))
    lDataReso.SetNDC(); lDataReso.SetTextFont(43); lDataReso.SetTextSize(20); lDataReso.Draw()
                        
    lMCReso = ROOT.TLatex(0.20,0.70,"#sigma^{Sim.}_{"+xLabel+"} = %.1f #mum" % (resoMC))
    lMCReso.SetNDC(); lMCReso.SetTextFont(43); lMCReso.SetTextSize(20); lMCReso.Draw()

    lBin = ROOT.TLatex(0.43,0.20,"Bin = %.1f #mum" % (sysErrData))
    lBin.SetNDC(); lBin.SetTextFont(43); lBin.SetTextSize(14); lBin.Draw()
    
    if ip != 'bs':
        
        lSelPV = ROOT.TLatex(0.20,0.85,str(vtrk['bins'][1])+' < '+selName+' < '+str(vtrk['bins'][2])+units)
        lSelPV.SetTextSize(0.035)
        lSelPV.SetNDC()
        lSelPV.Draw()
        
    else:

        lSelPVx = ROOT.TLatex(0.20,0.87,'Beam width (x) = %.1f #mum' % (vtrk['beamwidthx'][int(ktrkstr)]))
        lSelPVy = ROOT.TLatex(0.20,0.83,'Beam width (y) = %.1f #mum' % (vtrk['beamwidthy'][int(ktrkstr)]))
        lSelPVx.SetTextSize(0.032); lSelPVy.SetTextSize(0.032)
        lSelPVx.SetNDC(); lSelPVy.SetNDC()
        lSelPVx.Draw(); lSelPVy.Draw()
        
    pLabel = 'p_{T}'
    pUnits = 'GeV'
    pPrec = '%.2f'
    if param == 'eta':
        pLabel = '#eta'
        pUnits = ''
    elif param == 'phi':
        pLabel = '#phi'
        pUnits = ''
    elif param == 'npv':
        pLabel = 'N_{PV}'
        pUnits = ''
        pPrec = '%d'
    elif param == 'dr':
        pLabel = '#DeltaR'
        pUnits = ''

    lSel = ROOT.TLatex(0.20,0.78,(pPrec+' < '+pLabel+' < '+pPrec+' '+pUnits) % (v['bins'][1], v['bins'][2]))
    lSel.SetTextSize(0.035)
    lSel.SetNDC()
    lSel.Draw()
        
    c1.Update()

    leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
    leg.SetFillColor(253)
    leg.SetBorderSize(0)
    leg.AddEntry(hResoData,"Data","p")
    leg.AddEntry(hResoMC,"Simulation","f")
    if options.fit:
        leg.AddEntry(resResoData,"Data (fit)","l")
        leg.AddEntry(resResoMC,"Simulation (fit)","l")
    leg.Draw()        

    t1, t2, t3, t4 = style.cmslabel(1, c.year, evt)
    t1.Draw()
    t2.Draw()
    t3.Draw()
    t4.Draw()
            
#    b = fun.isbadfit(resoChi2MC, resoChi2Data)
#    b.Draw()

    foutput = 'ip'+ip+'Reso_'+x+ktrkstr+kstr
    if ip == 'bs': foutput = 'ip'+ip+'Reso_'+x+kstr+'_'+ktrkstr
    figs.append(foutput)
    c1.Print(options.output+'/'+foutput+'.'+img)
    c1.Clear()

    return figs, x, ktrkstr, kstr, resoData, resoErrData, intResoData, resoMC, resoErrMC, intResoMC, sysErrData, sysErrMC

if __name__ == '__main__':
    
    options = main()
    
    pool = multiprocessing.Pool(c.ncores)

    ROOT.gROOT.SetBatch()
    
    pstyle = style.SetPlotStyle(1, 'MAIN')

    figs = []
    
    if not os.path.isdir(options.output):
        os.system("mkdir "+options.output)

    ip = options.type

    evt = 'zb'
    if options.qcd: evt = 'qcd'

    ppath = '/user/kskovpen/analysis/Track/CMSSW_10_5_0_pre2/src/TrackingAnalysis/EDAnalyzers/macro/data/bins/'
    param = {}
    param['bs'] = fun.param(ppath+'zb_bs.json') if evt == 'zb' else fun.param(ppath+'qcd_bs.json')
    param['bsw'] = fun.param(ppath+'zb_bsw.json') if evt == 'zb' else fun.param(ppath+'qcd_bsw.json')
    param['pv'] = fun.param(ppath+'zb_pv.json') if evt == 'zb' else fun.param(ppath+'qcd_pv.json')
    
    ppath = 'data/bins/'    
    fparam = fun.param(ppath+evt+'_pv.json')
    parampv = fparam.get(options.parampv)
    
#    ipParamList = ['pt', 'eta', 'phi', 'npv', 'dr']
    ipParamList = ['eta']

    ParamList = {}
    
    for t in ['bs', 'pv']:
        ParamList[t] = {}        
        for p in ipParamList:
            ParamList[t][p] = param[t].get(p)
            del ParamList[t][p]['allbins']
    
    for t in ['bsw']:
        ParamList[t] = {}
        for p in ['beamwidthx', 'beamwidthy']:
            ParamList[t][p] = param[t].get(p)
            
    img = options.image
    
    fHistData = ROOT.TFile.Open(options.data,'read')
    fHistMC = ROOT.TFile.Open(options.mc,'read')

    print 'Plot cumulative distributions'
    
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
        
        c1.Print(options.output+'/'+h+'.eps')
        c1.Clear()
    
    pstyle.Reset()
    
    #### Fits
    
    print 'Run IP fits'
    
    if 'bs' in ip: c.IPmeas.remove('dz')
    
    rout = {}
    rout['reso'] = {}
    for t in ['data','mc']:
        rout['reso'][t] = {}
        for x in c.IPmeas:
            rout['reso'][t][x] = {}
            for pip in ipParamList:
                bins = ParamList[ip][pip]
                for k, v in bins.iteritems():
                    rout['reso'][t][x][k] = {}
                    if ip == 'pv':
                        for ktrk, vtrk in parampv.iteritems():
                            rout['reso'][t][x][k][ktrk] = {}
                    elif ip == 'bs':
                        for itrk, ktrk in enumerate(ParamList['bsw']['beamwidthx']):
                            rout['reso'][t][x][k][str(itrk)] = {}
                            
    jobs = []

    selName = 'N_{trk}'
    units = ''
    if options.parampv == 'sumTrackPt':
        selName = '#sump_{T}'
        units = ' GeV'
    elif options.parampv == 'sumTrackPtSq': 
        selName = '#sqrt{#sump^{2}_{T}}'
        units = ' GeV'

    if ip == 'pv':
        
        for ktrk, vtrk in parampv.iteritems():

            ktrkstr = str(ktrk)
        
            if ktrkstr in ['allbins', '']: continue

            for x in c.IPmeas:
        
                for pip in ipParamList:
            
                    bins = ParamList[ip][pip]
            
                    for k, v in bins.iteritems():
                    
                        if k in ['allbins', '']: continue
                    
                        kstr = str(k)

                        hNameResoData = 'h_ip'+ip+x+ktrkstr+kstr+options.selection
                        hResoData = fHistData.Get(hNameResoData).Clone('hResoData')
                        
                        hNameResoMC = 'h_ip'+ip+x+ktrkstr+kstr+options.selection
                        hResoMC = fHistMC.Get(hNameResoMC).Clone('hResoMC')
                        
                        jobs.append( pool.apply_async(runFit, (evt, ip, vtrk, v, x, ktrkstr, kstr, pip, img, hResoData, hResoMC)) )

    elif ip == 'bs':
        
        for itrk, ktrk in enumerate(ParamList['bsw']['beamwidthx']):

            ktrkstr = str(itrk)

            for x in c.IPmeas:
        
                for pip in ipParamList:
            
                    bins = ParamList[ip][pip]
            
                    for k, v in bins.iteritems():
                    
                        kstr = str(k)
                        
                        if kstr in ['allbins', '']: continue

                        hNameResoData = 'h_ip'+ip+x+kstr+'_'+ktrkstr+options.selection
                        hResoData = fHistData.Get(hNameResoData).Clone('hResoData')
                       
                        hNameResoMC = 'h_ip'+ip+x+kstr+options.selection
                        hResoMC = fHistMC.Get(hNameResoMC).Clone('hResoMC')
                        
                        jobs.append( pool.apply_async(runFit, (evt, ip, ParamList['bsw'], v, x, ktrkstr, kstr, pip, img, hResoData, hResoMC)) )
                        
    pool.close()

    for job in jobs: 

        result = job.get()
        
        figs += result[0]
        
        x = result[1]
        ktrkstr = result[2]
        kstr = result[3]
        
        resoData = result[4]
        resoErrData = result[5]
        intResoData = result[6]
        resoMC = result[7]
        resoErrMC = result[8]
        intResoMC = result[9]

        sysErrData = result[10]
        sysErrMC = result[11]
        
        rout['reso']['data'][x][kstr][ktrkstr].update([('value',resoData), ('error',resoErrData), ('int',intResoData), ('sys',sysErrData)])
        rout['reso']['mc'][x][kstr][ktrkstr].update([('value',resoMC), ('error',resoErrMC), ('int',intResoMC), ('sys',sysErrMC)])

    if not os.path.isdir('results'): os.system('mkdir results/')
                        
    foutput = 'results/ip_zb'+options.selection+'.json'
    if options.qcd: foutput = 'results/ip_qcd'+options.selection+'.json'
    with open(foutput, "w") as write_file:
        json.dump(rout, write_file, indent=2)

    print 'Generate thumbnails'
    fun.generateThumbs(figs, options.output, img)
    print 'Put plots on a web page'
    fun.createPage('ip'+options.type)
