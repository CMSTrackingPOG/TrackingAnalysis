import os
import sys
import math
import subprocess
import common as c
from subprocess import call
import style
import functions as func
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
    parser.add_option("-i","--input",default="results/pv.json",help="input file name [default: %default]")
    parser.add_option("-o","--output",default="pics",help="output directory [default: %default]")
    parser.add_option("-m","--mode",default="pv",help="measurement mode [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

def plot(c1, hData, hMC, mode, m, x, isDeconv = False):
    
    hMC.Draw('e1')
    hData.Draw('e1 same')

    hData.SetMarkerStyle(20)
    hData.SetMarkerSize(0.7)
    hData.SetMarkerColor(1)
    hData.SetLineColor(1)
    
    hMC.SetMarkerStyle(21)
    hMC.SetMarkerSize(0.7)
    hMC.SetMarkerColor(ROOT.kBlue-10)
    hMC.SetLineColor(ROOT.kBlue-10)
    
    if mode == 'pv':
        if m == 'reso' and x != 'z': hMC.GetYaxis().SetRangeUser(0.,100.)
        elif m == 'reso' and x == 'z': hMC.GetYaxis().SetRangeUser(0.,100.)
        else: hMC.GetYaxis().SetRangeUser(0.,1.4)
    else:
        hMC.GetYaxis().SetRangeUser(0.,350.)
            
    leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
    leg.SetFillColor(253)
    leg.SetBorderSize(0)
    leg.AddEntry(hData,"Data","p")
    leg.AddEntry(hMC,"Simulation","p")
    leg.Draw()        
    
    t1, t2 = style.cmslabel(1,777)
    t1.Draw()

    if isDeconv == False:
        c1.Print(options.output+'/'+mode+'_'+m+'_'+x+'.pdf')
    else:
        c1.Print(options.output+'/'+mode+'_'+m+'_'+x+'_deconv.pdf')
    c1.Clear()    

if __name__ == '__main__':
    
    options = main()
    
    mode = options.mode

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(1)
    pstyle.SetErrorX(0.5)

#    if os.path.isdir(options.output):
#        os.system("rm -rf "+options.output)                
#    os.system("mkdir "+options.output)

    with open(options.input, "r") as read_file:
        data = json.load(read_file)        

    with open('results/pv.json', "r") as read_file:
        datapv = json.load(read_file)
        
    c1 = ROOT.TCanvas()
    
    h = {}
            
    meas = c.PVmeas
    if mode != 'pv': meas = c.IPmeas
    
    typ = ['reso','pull']
    if mode != 'pv': typ = ['reso']

    addParam = c.PVnTracks if mode != 'pv' else {'':[]}
    
    for m in typ:
        for x in meas:
            
            for t in ['data','mc']:
                
                hname = 'h_'+m+'_'+t+'_'+x
                
                if mode == 'pv':
                    h[hname] = ROOT.TH1F(hname,hname,len(c.PVnTracksBins)-1,c.PVnTracksBins)
                    h[hname].GetXaxis().SetTitle('Number of tracks')
                    ytit = 'PV resolution in '+x+' [#mum]'
                    if m == 'pull': ytit = 'Primary vertex pull in '+x
                    h[hname].GetYaxis().SetTitle(ytit)
                    h[hname].Sumw2()
                else:                        
                    for hn in [hname,hname+'_deconv']:
                        if 'pt' in mode:
                            h[hn] = ROOT.TH1F(hn,hn,len(c.IPptBins)-1,c.IPptBins)
                            h[hn].GetXaxis().SetTitle('Track p_{T} [GeV]')
                        elif 'eta' in mode:
                            h[hn].GetXaxis().SetTitle('Track #eta')
                            h[hn] = ROOT.TH1F(hn,hn,len(c.IPetaBins)-1,c.IPetaBins)
                            
                        lab = x.replace('0','_{xy}').replace('z','_{z}')
                        ytit = 'Track IP resolution ('+lab+') [#mum]'
                        h[hn].GetYaxis().SetTitle(ytit)
                        h[hn].Sumw2()

            hnameData = 'h_'+m+'_data_'+x
            hnameMC = 'h_'+m+'_mc_'+x
            
            param = c.PVnTracks
            if 'pt' in mode: param = c.IPpt
            elif 'eta' in mode: param = c.IPeta
            
            for kparam, vparam in param.iteritems():

                if mode == 'pv':
                        
                    resData = data[m]['data'][x][kparam]
                    resMC = data[m]['mc'][x][kparam]
                
                    if bool(resData) and bool(resMC):
                    
                        vData = resData['value']
                        eData = resData['error']

                        vMC = resMC['value']
                        eMC = resMC['error']
                    
                        bidx = param[kparam]['bins'][0]

                        h[hnameData].SetBinContent(bidx, vData)
                        h[hnameData].SetBinError(bidx, eData)

                        h[hnameMC].SetBinContent(bidx, vMC)
                        h[hnameMC].SetBinError(bidx, eMC)
                            
                else:
                        
                    vData = 0; vMC = 0;
                    eData = 0; eMC = 0;                        
                    vDataDeconv = 0; vMCDeconv = 0;
                    eDataDeconv = 0; eMCDeconv = 0;                        
                    wDataSum = 0; wMCSum = 0;
                    
                    for ktrk, vtrk in c.PVnTracks.iteritems():
                   
                        if 'ipbs' in mode and (ktrk != '' or kparam == ''): break

                        resData = data[m]['data'][x][kparam][ktrk]
                        resMC = data[m]['mc'][x][kparam][ktrk]
                        if 'ippv' in mode:
                            pvx = 'x'
                            if x == 'dz': pvx = 'z'                        
                            resPVData = datapv[m]['data'][pvx][ktrk]
                            resPVMC = datapv[m]['mc'][pvx][ktrk]
                        elif 'ipbs' in mode:
                            pvx = 'x'
                            if x == 'dz': pvx = 'z'
                            resPVData = datapv['widthx']['data'][pvx][ktrk]
                            resPVMC = datapv['widthx']['mc'][pvx][ktrk]                            

                        if (bool(resData) and bool(resMC) and ktrk != '') or 'ipbs' in mode:

#                            print kparam, ktrk, resData, resMC
                            vDataTrk = resData['value']
                            eDataTrk = resData['error']
                            iDataTrk = resData['int']
                            vPVDataTrk = resPVData['value']
                            ePVDataTrk = resPVData['error']
                            
                            vMCTrk = resMC['value']
                            eMCTrk = resMC['error']
                            iMCTrk = resMC['int']
                            vPVMCTrk = resPVMC['value']
                            ePVMCTrk = resPVMC['error']
                                
                            wDataTrk = iDataTrk
                            wMCTrk = iMCTrk
                            
                            wDataSum += wDataTrk
                            wMCSum += wMCTrk

                            # conv
                            vData += vDataTrk*wDataTrk
                            vMC += vMCTrk*wMCTrk

                            eData += math.pow(eDataTrk*wDataTrk,2)
                            eMC += math.pow(eMCTrk*wMCTrk,2)

                            # deconv
                            if vDataTrk < vPVDataTrk or vMCTrk < vPVMCTrk:
                                print ktrk, kparam, vDataTrk, vPVDataTrk, vMCTrk, vPVMCTrk
                                continue
                            vSigmaData = math.sqrt(vDataTrk*vDataTrk-vPVDataTrk*vPVDataTrk)
                            vDataDeconv += vSigmaData*wDataTrk
                            vSigmaMC = math.sqrt(vMCTrk*vMCTrk-vPVMCTrk*vPVMCTrk)
                            vMCDeconv += vSigmaMC*wMCTrk
                            
                            eSigmaData = math.sqrt(math.pow(vDataTrk*eDataTrk,2)+math.pow(vPVDataTrk*ePVDataTrk,2))/vSigmaData
                            eDataDeconv += math.pow(eSigmaData*wDataTrk,2)
                            eSigmaMC = math.sqrt(math.pow(vMCTrk*eMCTrk,2)+math.pow(vPVMCTrk*ePVMCTrk,2))/vSigmaMC
                            eMCDeconv += math.pow(eSigmaMC*wMCTrk,2)
                            
                        else:
                            continue

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
            pickle.dump(hMC,open('results/'+mode+'_'+m+'_'+x+'_mc.pkl','wb'))
            pickle.dump(hData,open('results/'+mode+'_'+m+'_'+x+'_data.pkl','wb'))
            plot(c1, hData, hMC, mode, m, x)
            
            hMC = h[hnameMC+'_deconv']; hData = h[hnameData+'_deconv']
            pickle.dump(hMC,open('results/'+mode+'_'+m+'_'+x+'_mc_deconv.pkl','wb'))
            pickle.dump(hData,open('results/'+mode+'_'+m+'_'+x+'_data_deconv.pkl','wb'))
            plot(c1, hData, hMC, mode, m, x, True)
