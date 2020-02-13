import os
import sys
import math
import subprocess
import common as c
from subprocess import call
import array
import json
import pickle
import ROOT

import style
import fit

ROOT.PyConfig.IgnoreCommandLineOptions = True
from optparse import OptionParser

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]

    usage = "usage: %prog [options]\n Plotting script to overlay histograms"
    
    parser = OptionParser(usage)
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(1)
    pstyle.SetErrorX(0.5)
        
    c1 = ROOT.TCanvas()
    
    h = {}

    leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
    leg.SetFillColor(253)
    leg.SetBorderSize(0)
    
    lst = [20,20,21,21]
    lab = ['IPBS (Data)','IPBS (MC)','IPPV (Data)','IPPV (MC)']
    hname = ['ipbspt_reso_d0_data_deconv','ipbspt_reso_d0_mc_deconv',\
    'ippvpt_reso_d0_data_deconv','ippvpt_reso_d0_mc_deconv']
    
    for i, hn in enumerate(hname):

        h[hn] = pickle.load(open('results/'+hn+'.pkl','rb'))

        h[hn].SetMarkerStyle(lst[i])
        h[hn].SetMarkerSize(0.7)
        
        if 'data' in hn:
            h[hn].SetMarkerColor(1)
            h[hn].SetLineColor(1)
        else:    
            h[hn].SetMarkerColor(ROOT.kBlue-10)
            h[hn].SetLineColor(ROOT.kBlue-10)

        if i == 0: h[hn].Draw('hist e1p')
        else: h[hn].Draw('hist e1p same')
            
        leg.AddEntry(h[hn],lab[i],'p')

    leg.Draw()
    
    h[hname[0]].GetYaxis().SetRangeUser(0,350)
    
    t1, t2 = style.cmslabel(1,777)
    t1.Draw()
    c1.Print('pics/comp.pdf')
    c1.Clear()    

    hfit = h['ipbspt_reso_d0_data_deconv']

    hfit.Draw('e1p')
    
    x = array.array('d')
    y = array.array('d')
    for b in range(hfit.GetXaxis().GetNbins()):
        x.append(hfit.GetXaxis().GetBinCenter(b+1))
        y.append(hfit.GetBinContent(b+1))
        
    gr = ROOT.TGraph(int(len(x)),x,y)
        
    res, p0, p1, chi2 = fit.doFitIP('ipfit',gr,ROOT.kRed)
    print p0, p1
    res.Draw('same')
    
    t1, t2 = style.cmslabel(1,777)
    t1.Draw()
    c1.Print('pics/ipfit.pdf')
    c1.Clear()    
    
