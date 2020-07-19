import os
import sys
import math
import subprocess
import common as c
from subprocess import call
import style
import functions as fun
import numpy as np
import pickle
import json
import ROOT

import fit

ROOT.PyConfig.IgnoreCommandLineOptions = True
from optparse import OptionParser

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]

    usage = "usage: %prog [options]\n Analysis script to plot the PV resolution results"
    
    parser = OptionParser(usage)
    parser.add_option("--input", default="results/ip.json", help="input file name for IP [default: %default]")
    parser.add_option("--inputpv", default="results/pv.json", help="input file name for PV [default: %default]")
    parser.add_option("--output", default="pics", help="output directory [default: %default]")
    parser.add_option("--mode", default="pv", help="measurement mode [default: %default]")
    parser.add_option("--process", default="zb", help="type of process [default: %default]")
    parser.add_option("--parampv", default="sumTrackPtSq", help="Parameterisation for PV resolution measurement [default: %default]")
    parser.add_option("--sys", action='store_true', help="Add bin width as systematics [default: %default]")
    parser.add_option("--selection", default="", help="Selection criteria in IP measurement [default: %default]")
    parser.add_option("--pub", action='store_true', help="Publish plots [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

def plot(c1, hData, hMC, mode, m, x, isDeconv = False):

    hData.SetMarkerStyle(20)
    hData.SetMarkerSize(1.0)
    hData.SetMarkerColor(1)
    hData.SetLineColor(1)
    
    hMC.SetMarkerStyle(22)
    hMC.SetMarkerSize(1.0)
    hMC.SetMarkerColor(c.mcfit)
    hMC.SetLineColor(c.mcfit)
    
    if mode == 'pv':
        if m == 'reso':
            if options.process == 'zb': h0 = c1.DrawFrame(0., 0., 20., 70.)
            else: h0 = c1.DrawFrame(0., 0., 400., 30.)
        else:
            if options.process == 'zb': h0 = c1.DrawFrame(0., 0., 20., 1.7)
            else: h0 = c1.DrawFrame(0., 0., 400., 1.7)
    
    if mode == 'pv': 
        hMC.Draw('same')
        h0.GetXaxis().SetTitle(hMC.GetXaxis().GetTitle())
        h0.GetYaxis().SetTitle(hMC.GetYaxis().GetTitle())
    else: hMC.Draw('')
    hMC.Draw('pe1 same')
    hData.Draw('same')
    hData.Draw('pe1 same')
    
    if mode != 'pv':
        hMC.GetYaxis().SetRangeUser(0.,450.)
        if 'eta' in mode and x == 'dz': hMC.GetYaxis().SetRangeUser(0.,3000.)
            
    leg = ROOT.TLegend(0.62,0.55,0.86,0.73)
    leg.SetFillColor(253)
    leg.SetBorderSize(0)
    leg.AddEntry(hData,"Data","p")
    leg.AddEntry(hMC,"Simulation","p")
    
    if mode == 'pv' and m == 'pull':
        leg = ROOT.TLegend(0.67,0.25,0.910,0.43)
        leg.SetFillColor(253)
        leg.SetBorderSize(0)
        leg.AddEntry(hData,"Data","p")
        leg.AddEntry(hMC,"Simulation","p")
    
    leg.Draw()
    
    t1, t2, t3, t4 = style.cmslabel(2, c.year, evt, False)
    t1.Draw()
    t2.Draw()
    t3.Draw()
    t4.Draw()

    if isDeconv == False:
        for outdir in [options.output, 'pub']:
            if not options.pub and outdir == 'pub': continue
            c1.Print(outdir+'/'+mode+'_'+m+'_'+x+'_'+options.process+'.pdf')
        if options.pub: c1.SaveAs('pub/'+mode+'_'+m+'_'+x+'_'+options.process+'.root')
    else:
        for outdir in [options.output, 'pub']:
            if not options.pub and outdir == 'pub': continue
            c1.Print(outdir+'/'+mode+'_'+m+'_'+x+'_deconv_'+options.process+'.pdf')
        if options.pub: c1.SaveAs('pub/'+mode+'_'+m+'_'+x+'_deconv_'+options.process+'.root')
    c1.Clear()    

if __name__ == '__main__':
    
    options = main()
    
    mode = options.mode

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(2)
#    pstyle.SetErrorX(0.5)

#    if os.path.isdir(options.output):
#        os.system("rm -rf "+options.output)

#    os.system("mkdir "+options.output)

    if os.path.isfile(options.input):
        with open(options.input, "r") as read_file:
            data = json.load(read_file)
    elif mode != 'pv':
        print 'Input file', options.input, 'not found'

    if os.path.isfile(options.inputpv):
        with open(options.inputpv, "r") as read_file:
            datapv = json.load(read_file)
    else:
        print 'Input file', options.inputpv, 'not found'
        
    c1 = ROOT.TCanvas()
    
    h = {}
            
    meas = c.PVmeas
    if mode != 'pv': meas = c.IPmeas
    if 'bs' in mode: meas.remove('dz')
    
    typ = ['reso','pull']
    if mode != 'pv': typ = ['reso']

    evt = options.process
    
    ppath = 'data/bins/'
    fparampv = fun.param(ppath+evt+'_pv.json')
#    fparampv = fun.param(ppath+evt+'_bs.json')
    parampv = fparampv.get(options.parampv)
    parambinspv = np.array(list(parampv['allbins'].replace('[','').replace(']','').split(',')), dtype='float64')
    
    ppath = '/user/kskovpen/analysis/Track/CMSSW_10_5_0_pre2/src/TrackingAnalysis/EDAnalyzers/macro/data/bins/'
    param = {}
    param['bs'] = fun.param(ppath+'zb_bs.json') if evt == 'zb' else fun.param(ppath+'qcd_bs.json')
    param['bsw'] = fun.param(ppath+'zb_bsw.json') if evt == 'zb' else fun.param(ppath+'qcd_bsw.json')
    param['pv'] = fun.param(ppath+'zb_pv.json') if evt == 'zb' else fun.param(ppath+'qcd_pv.json')

    paramip = None
    parambsw = {}
    if mode != 'pv':
        if 'ippv' in mode: paramip = param['pv'].get(mode.replace('ippv',''))
        elif 'ipbs' in mode: 
            paramip = param['bs'].get(mode.replace('ipbs',''))
            parambsw['x'] = param['bsw'].get('beamwidthx')
            parambsw['y'] = param['bsw'].get('beamwidthy')
        parambinsip = np.array(list(paramip['allbins'].replace('[','').replace(']','').split(',')), dtype='float64')
    
    for m in typ:
        for x in meas:
            
            for t in ['data','mc']:
                
                hname = 'h_'+m+'_'+t+'_'+x
                
                if mode == 'pv':
                    h[hname] = ROOT.TH1F(hname, hname, len(parambinspv)-1, parambinspv)
                    if options.parampv == 'nTracks': h[hname].GetXaxis().SetTitle('Number of tracks')
                    elif options.parampv == 'sumTrackPt': h[hname].GetXaxis().SetTitle('#sum p_{T} [GeV]')
                    elif options.parampv == 'sumTrackPtSq': h[hname].GetXaxis().SetTitle('#sqrt{#sum p_{T}^{2}} [GeV]')
                    ytit = 'PV resolution in '+x+' [#mum]'
                    if m == 'pull': ytit = 'Primary vertex pull in '+x
                    h[hname].GetYaxis().SetTitle(ytit)
                    h[hname].Sumw2()
                else:                        
                    for hn in [hname,hname+'_deconv']:
                        h[hn] = ROOT.TH1F(hn, hn, len(parambinsip)-1, parambinsip)
                        if 'pt' in mode: h[hn].GetXaxis().SetTitle('Track p_{T} [GeV]')
                        elif 'eta' in mode: h[hn].GetXaxis().SetTitle('Track #eta')
                        elif 'phi' in mode: h[hn].GetXaxis().SetTitle('Track #phi')
                        elif 'npv' in mode: h[hn].GetXaxis().SetTitle('Number of primary vertices')
                        elif 'dr' in mode: h[hn].GetXaxis().SetTitle('#DeltaR(track,jet axis)')
                            
                        lab = x.replace('0','_{xy}').replace('z','_{z}')
                        ytit = 'Track IP resolution ('+lab+') [#mum]'
                        h[hn].GetYaxis().SetTitle(ytit)
                        h[hn].Sumw2()

            hnameData = 'h_'+m+'_data_'+x
            hnameMC = 'h_'+m+'_mc_'+x
            
            param = parampv
            if 'ip' in mode: param = paramip
            
            for kparam, vparam in param.iteritems():
                
                if kparam in ['allbins', '']: continue

                if mode == 'pv':
                    
                    resData = datapv[m]['data'][x][kparam]
                    resMC = datapv[m]['mc'][x][kparam]
                
                    if bool(resData) and bool(resMC):
                    
                        vData = resData['value']
                        eData = resData['error']
                        sData = resData['sys']
                        if options.sys: esData = math.sqrt(eData*eData+sData*sData)
                        else: esData = math.sqrt(eData*eData)

                        vMC = resMC['value']
                        eMC = resMC['error']
                        sMC = resMC['sys']
                        if options.sys: esMC = math.sqrt(eMC*eMC+sMC*sMC)
                        else: esMC = math.sqrt(eMC*eMC)
                    
                        bidx = param[kparam]['bins'][0]

                        h[hnameData].SetBinContent(bidx, vData)
                        h[hnameData].SetBinError(bidx, esData)

                        h[hnameMC].SetBinContent(bidx, vMC)
                        h[hnameMC].SetBinError(bidx, esMC)
                            
                else:
                        
                    if 'ippv' in mode:
                        
                        vData = 0; vMC = 0;
                        eData = 0; eMC = 0;
                        vDataDeconv = 0; vMCDeconv = 0;
                        eDataDeconv = 0; eMCDeconv = 0;                        
                        wDataSum = 0; wMCSum = 0;
                        
                        for ktrk, vtrk in parampv.iteritems():
                        
                            if ktrk in ['allbins', '']: continue

                            resData = data[m]['data'][x][kparam][ktrk]
                            resMC = data[m]['mc'][x][kparam][ktrk]
                            pvx = 'x'
                            if x == 'dz': pvx = 'z'
                            resPVData = datapv[m]['data'][pvx][ktrk]
                            resPVMC = datapv[m]['mc'][pvx][ktrk]

                            vDataTrk = resData['value']
                            eDataTrk = resData['error']
                            sDataTrk = resData['sys']
                            iDataTrk = resData['int']
                            vPVDataTrk = resPVData['value']
                            ePVDataTrk = resPVData['error']
                            sPVDataTrk = resPVData['sys']
                            
                            vMCTrk = resMC['value']
                            eMCTrk = resMC['error']
                            sMCTrk = resMC['sys']
                            iMCTrk = resMC['int']
                            vPVMCTrk = resPVMC['value']
                            ePVMCTrk = resPVMC['error']
                            sPVMCTrk = resPVMC['sys']
                                
                            wDataTrk = float(iDataTrk)
                            wMCTrk = float(iMCTrk)                                                        
                            
                            if vDataTrk*vMCTrk == 0.: continue
                            
                            wDataSum += wDataTrk
                            wMCSum += wMCTrk

                            # conv
                            vData += vDataTrk*wDataTrk
                            vMC += vMCTrk*wMCTrk

                            if options.sys:
                                eData += math.pow(eDataTrk*wDataTrk,2) + math.pow(sDataTrk*wDataTrk,2)
                                eMC += math.pow(eMCTrk*wMCTrk,2) + math.pow(sMCTrk*wMCTrk,2)
                            else:
                                eData += math.pow(eDataTrk*wDataTrk,2)
                                eMC += math.pow(eMCTrk*wMCTrk,2)

                            # deconv
                            if (vDataTrk < vPVDataTrk or vMCTrk < vPVMCTrk):
                                if vDataTrk != -1 and vMCTrk != -1:
                                    print x, ktrk, kparam, vDataTrk, vPVDataTrk, vMCTrk, vPVMCTrk
                                continue
                            
                            vSigmaData = math.sqrt(vDataTrk*vDataTrk-vPVDataTrk*vPVDataTrk)
                            vDataDeconv += vSigmaData*wDataTrk
                            vSigmaMC = math.sqrt(vMCTrk*vMCTrk-vPVMCTrk*vPVMCTrk)
                            vMCDeconv += vSigmaMC*wMCTrk
                            
                            if options.sys:
                                eSigmaData = math.sqrt(math.pow(vDataTrk*eDataTrk,2)+math.pow(vDataTrk*sDataTrk,2)+math.pow(vPVDataTrk*ePVDataTrk,2)+math.pow(vPVDataTrk*sPVDataTrk,2))/vSigmaData
                                eSigmaMC = math.sqrt(math.pow(vMCTrk*eMCTrk,2)+math.pow(vMCTrk*sMCTrk,2)+math.pow(vPVMCTrk*ePVMCTrk,2)+math.pow(vPVMCTrk*sPVMCTrk,2))/vSigmaMC
                            else:
                                eSigmaData = math.sqrt(math.pow(vDataTrk*eDataTrk,2)+math.pow(vPVDataTrk*ePVDataTrk,2))/vSigmaData
                                eSigmaMC = math.sqrt(math.pow(vMCTrk*eMCTrk,2)+math.pow(vPVMCTrk*ePVMCTrk,2))/vSigmaMC
                                
                            eDataDeconv += math.pow(eSigmaData*wDataTrk,2)                            
                            eMCDeconv += math.pow(eSigmaMC*wMCTrk,2)

                        if wDataSum == 0 or wMCSum == 0: continue

                        vData /= wDataSum
                        vMC /= wMCSum
                        
                        eData = math.sqrt(eData)/wDataSum
                        eMC = math.sqrt(eMC)/wMCSum

                        vDataDeconv /= wDataSum
                        vMCDeconv /= wMCSum
                        
                        eDataDeconv = math.sqrt(eDataDeconv)/wDataSum
                        eMCDeconv = math.sqrt(eMCDeconv)/wMCSum
                    
                        bidx = param[kparam]['bins'][0]
                    
                        h[hnameData].SetBinContent(bidx, vData)
                        h[hnameData].SetBinError(bidx, eData)
                        
                        h[hnameMC].SetBinContent(bidx, vMC)
                        h[hnameMC].SetBinError(bidx, eMC)

                        h[hnameData+'_deconv'].SetBinContent(bidx, vDataDeconv)
                        h[hnameData+'_deconv'].SetBinError(bidx, eDataDeconv)
                        
                        h[hnameMC+'_deconv'].SetBinContent(bidx, vMCDeconv)
                        h[hnameMC+'_deconv'].SetBinError(bidx, eMCDeconv)

                    elif 'ipbs' in mode:
                        
                        vData = 0; vMC = 0;
                        eData = 0; eMC = 0;
                        vDataDeconv = 0; vMCDeconv = 0;
                        eDataDeconv = 0; eMCDeconv = 0;                        
                        wDataSum = 0; wMCSum = 0;
                        
                        bwx = datapv['widthx']['mc']['x'].itervalues()
                        bwy = datapv['widthy']['mc']['y'].itervalues()
                        bwx.next()
                        bwy.next()
#                        beamWidthXMC = datapv['widthx']['mc']['x'].itervalues().next()['value']
#                        beamWidthYMC = datapv['widthy']['mc']['y'].itervalues().next()['value']
                        beamWidthXMC = bwx.next()['value']
                        beamWidthYMC = bwy.next()['value']
                        beamWidthXYMC = (beamWidthXMC+beamWidthYMC)/2.
                        beamWidthErrorXYMC = abs(beamWidthXYMC-beamWidthXMC)
                        
                        for itrk, ktrk in enumerate(parambsw['x']):

                            resData = data[m]['data'][x][kparam][str(itrk)]
                            resMC = data[m]['mc'][x][kparam].itervalues().next()
                            resXBSData = ktrk
                            resYBSData = parambsw['y'][itrk]
                            resXYBSData = (resXBSData+resYBSData)/2.
                            resErrorXYBSData = abs(resXYBSData-resXBSData)

                            vDataTrk = resData['value']
                            eDataTrk = resData['error']
                            sDataTrk = resData['sys']
                            iDataTrk = resData['int']
                            
                            vMCTrk = resMC['value']
                            eMCTrk = resMC['error']
                            sMCTrk = resMC['sys']
                            iMCTrk = resMC['int']

                            wDataTrk = float(iDataTrk)
                            wMCTrk = float(iMCTrk)
                            
                            if vDataTrk*vMCTrk == 0.: continue
                            
                            wDataSum += wDataTrk
                            wMCSum += wMCTrk

                            # conv
                            vData += vDataTrk*wDataTrk
                            vMC += vMCTrk*wMCTrk

                            if options.sys:
                                eData += math.pow(eDataTrk*wDataTrk,2) + math.pow(sDataTrk*wDataTrk,2)
                                eMC += math.pow(eMCTrk*wMCTrk,2) + math.pow(sMCTrk*wMCTrk,2)
                            else:
                                eData += math.pow(eDataTrk*wDataTrk,2)
                                eMC += math.pow(eMCTrk*wMCTrk,2)

                            # deconv
                            if (vDataTrk < resXYBSData or vMCTrk < beamWidthXYMC):
                                if vDataTrk != -1 and vMCTrk != -1:
                                    print ktrk, kparam, vDataTrk, resXYBSData, vMCTrk, beamWidthXYMC
                                continue
                            
                            vSigmaData = math.sqrt(vDataTrk*vDataTrk-resXYBSData*resXYBSData)
                            vDataDeconv += vSigmaData*wDataTrk
                            vSigmaMC = math.sqrt(vMCTrk*vMCTrk-beamWidthXYMC*beamWidthXYMC)
                            vMCDeconv += vSigmaMC*wMCTrk

                            if options.sys:
                                eSigmaData = math.sqrt(math.pow(vDataTrk*eDataTrk,2)+math.pow(vDataTrk*sDataTrk,2)+math.pow(resXYBSData*resErrorXYBSData,2))/vSigmaData
                                eSigmaMC = math.sqrt(math.pow(vMCTrk*eMCTrk,2)+math.pow(vMCTrk*sMCTrk,2)+math.pow(beamWidthXYMC*beamWidthErrorXYMC,2))/vSigmaMC
                            else:
                                eSigmaData = math.sqrt(math.pow(vDataTrk*eDataTrk,2)+math.pow(resXYBSData*resErrorXYBSData,2))/vSigmaData
                                eSigmaMC = math.sqrt(math.pow(vMCTrk*eMCTrk,2)+math.pow(beamWidthXYMC*beamWidthErrorXYMC,2))/vSigmaMC
                                
                            eDataDeconv += math.pow(eSigmaData*wDataTrk,2)
                            eMCDeconv += math.pow(eSigmaMC*wMCTrk,2)

                        if wDataSum == 0 or wMCSum == 0: continue

                        vData /= wDataSum
                        vMC /= wMCSum
                        
                        eData = math.sqrt(eData)/wDataSum
                        eMC = math.sqrt(eMC)/wMCSum

                        vDataDeconv /= wDataSum
                        vMCDeconv /= wMCSum
                        
                        eDataDeconv = math.sqrt(eDataDeconv)/wDataSum
                        eMCDeconv = math.sqrt(eMCDeconv)/wMCSum
                    
                        bidx = param[kparam]['bins'][0]
                    
                        h[hnameData].SetBinContent(bidx, vData)
                        h[hnameData].SetBinError(bidx, eData)

                        h[hnameMC].SetBinContent(bidx, vMC)
                        h[hnameMC].SetBinError(bidx, eMC)

                        h[hnameData+'_deconv'].SetBinContent(bidx, vDataDeconv)
                        h[hnameData+'_deconv'].SetBinError(bidx, eDataDeconv)
                        
                        h[hnameMC+'_deconv'].SetBinContent(bidx, vMCDeconv)
                        h[hnameMC+'_deconv'].SetBinError(bidx, eMCDeconv)
                        
            hMC = h[hnameMC]; hData = h[hnameData]
            pickle.dump(hMC,open('results/'+mode+'_'+m+'_'+x+options.selection+'_'+evt+'_mc.pkl','wb'))
            pickle.dump(hData,open('results/'+mode+'_'+m+'_'+x+options.selection+'_'+evt+'_data.pkl','wb'))
            plot(c1, hData, hMC, mode, m, x)
            
            if mode != 'pv':
                
                hMC = h[hnameMC+'_deconv']; hData = h[hnameData+'_deconv']
                pickle.dump(hMC,open('results/'+mode+'_'+m+'_'+x+options.selection+'_'+evt+'_mc_deconv.pkl','wb'))
                pickle.dump(hData,open('results/'+mode+'_'+m+'_'+x+options.selection+'_'+evt+'_data_deconv.pkl','wb'))
                plot(c1, hData, hMC, mode, m, x, True)
