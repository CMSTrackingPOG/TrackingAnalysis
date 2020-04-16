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
    parser.add_option("-p","--param",default="pt",help="parametrization [default: %default]")
    parser.add_option("-q","--qcd",action='store_true',help="compare between ZeroBias and QCD [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(1)
    pstyle.SetErrorX(0.5)
        
    c1 = ROOT.TCanvas()
    c1.SetLogy(1)
    
    h = {}

    leg = ROOT.TLegend(0.82,0.92,0.990,0.75)
    leg.SetFillColor(253)
    leg.SetBorderSize(0)
    
    lst = [20,20,22,22]
    lab = ['IPBS (Data)','IPBS (Sim.)','IPPV (Data)','IPPV (Sim.)']
    if options.qcd:
        lab = ['ZeroBias (Data)','ZeroBias (Sim.)','QCD (Data)','QCD (Sim.)']
        
    hname = ['ipbs'+options.param+'_reso_d0_zb_data_deconv','ipbs'+options.param+'_reso_d0_zb_mc_deconv',\
    'ippv'+options.param+'_reso_d0_zb_data_deconv','ippv'+options.param+'_reso_d0_zb_mc_deconv']
    if options.qcd:
        hname = ['ippv'+options.param+'_reso_d0_zb_data_deconv','ippv'+options.param+'_reso_d0_zb_mc_deconv',\
        'ippv'+options.param+'_reso_d0_qcd_data_deconv','ippv'+options.param+'_reso_d0_qcd_mc_deconv']
    
    for i, hn in enumerate(hname):

        h[hn] = pickle.load(open('results/'+hn+'.pkl','rb'))

        h[hn].SetMarkerStyle(lst[i])
        h[hn].SetMarkerSize(1.0)
        
        if 'data' in hn:
            h[hn].SetMarkerColor(1)
            h[hn].SetLineColor(1)
        else:    
            h[hn].SetMarkerColor(c.mccol)
            h[hn].SetLineColor(c.mccol)

        if i == 0: 
            h[hn].Draw('')
            h[hn].Draw('p same')
        else: 
            h[hn].Draw('same')
            h[hn].Draw('p same')
            
        leg.AddEntry(h[hn],lab[i],'p')

    leg.Draw()

    h[hname[0]].GetYaxis().SetRangeUser(0,350)
    if c1.GetLogy(): h[hname[0]].GetYaxis().SetRangeUser(10,350)
    
    t1, t2, t3 = style.cmslabel(1,c.year)
    t1.Draw()
    t2.Draw()
    t3.Draw()
    c1.Print('pics/comp.pdf')
    c1.Clear()

    hfit = h['ippv'+options.param+'_reso_d0_zb_data_deconv']

    hfit.Draw('e1p')
    
    x = array.array('d')
    y = array.array('d')
    for b in range(hfit.GetXaxis().GetNbins()):
        x.append(hfit.GetXaxis().GetBinCenter(b+1))
        y.append(hfit.GetBinContent(b+1))
        
    gr = ROOT.TGraph(int(len(x)),x,y)
        
    c1.SetLogy(0)
    hfit.SetMinimum(0.)
    res, p0, p1, chi2 = fit.doFitIP('ipfit',gr,ROOT.kRed)
    print p0, p1
    res.Draw('same')
    
    t1, t2, t3 = style.cmslabel(1,c.year)
    t1.Draw()
    t2.Draw()
    t3.Draw()
    c1.Print('pics/ipfit.pdf')
    c1.Clear()    
    
