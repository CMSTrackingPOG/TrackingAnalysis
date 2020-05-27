#!/usr/bin/env python

import os
import sys
import math
import array
import common as c
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
    parser.add_option("--param", default='pt', help="parametrization [default: %default]")
    parser.add_option("--measurement", default='ip', help="ip or pv [default: %default]")
    parser.add_option("--type", default='bs', help="type of IP measurement (pv or bs) [default: %default]")
    parser.add_option("--pkl", default="reso_d0_zb_data,reso_d0_zb_mc", help="list of measurements [default: %default]")
    parser.add_option("--fit", action='store_true', help="Perform a pT-fit [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(2)

    files = options.pkl.split(',')
    
    c1 = ROOT.TCanvas()
    c1.SetLogx(1)

    h = {}

    leg = ROOT.TLegend(0.62,0.72,0.890,0.55)
    leg.SetFillColor(253)
    leg.SetBorderSize(0)
    
    for f in files:
        
        print '  ->', f

        pref = 'pv_'
        if options.measurement == 'ip': pref = options.measurement+options.type+options.param+'_'
            
        h[f] = pickle.load(open('results/'+pref+f+'.pkl','rb'))

        h[f].SetMarkerStyle(22)
        h[f].SetMarkerSize(0.7)
        
        if 'qcd' in f:
            h[f].SetMarkerStyle(20)
            h[f].SetMarkerColor(2)
            h[f].SetLineColor(2)
        
        if 'mc' in f: 
            if 'qcd' in f: h[f].SetMarkerStyle(24)
            else: h[f].SetMarkerStyle(26)
            h[f].SetMarkerColor(9)
            h[f].SetLineColor(9)

        lab = ''
        
        if 'zb' in f: lab += 'ZeroBias ('
        else: lab += ' QCD ('

        if 'data' in f: lab += 'Data)'
        else: lab += 'MC)'
        
        leg.AddEntry(h[f], lab, 'p')
    
    xmax = 0.; ymax = 0.
    xmin = 1E+10; ymin = 1E+10
    
    for k, v in h.iteritems():
        
        xmaxx = v.GetXaxis().GetXmax()
        if xmaxx > xmax: xmax = xmaxx
        xminn = v.GetXaxis().GetXmin()
        if xminn < xmin: xmin = xminn
    
        ymaxx = v.GetMaximum()
        if ymaxx > ymax: ymax = ymaxx
        yminn = v.GetMinimum()
        if yminn < ymin: ymin = yminn
       
    xmin = min(xmin, 0.1)
    xmax = xmax*2.0
    ymin = min(ymin, 0.)
    ymax = ymax*1.2
    
    h0 = c1.DrawFrame(xmin, ymin, xmax, ymax)

    if 'pt' in options.param: h0.GetXaxis().SetTitle('Track p_{T} [GeV]')
    elif 'eta' in options.param: h0.GetXaxis().SetTitle('Track #eta')
    elif 'phi' in options.param: h0.GetXaxis().SetTitle('Track #phi')
    elif 'npv' in options.param: h0.GetXaxis().SetTitle('Number of primary vertices')
    elif 'dr' in options.param: h0.GetXaxis().SetTitle('#DeltaR(jet axis,track)')
    
    if options.measurement == 'pv': h0.GetXaxis().SetTitle('#sqrt{#sum p_{T}^{2}} [GeV]')
    
    for i, f in enumerate(files):

        h[f].Draw('pe1 same')
        
        if 'd0' in f: h0.GetYaxis().SetTitle('Track IP resolution (d_{0}) [#mum]')
        elif 'dz' in f: h0.GetYaxis().SetTitle('Track IP resolution (d_{z}) [#mum]')
        elif 'reso_x' in f: h0.GetYaxis().SetTitle('PV resolution in x [#mum]')
        elif 'reso_y' in f: h0.GetYaxis().SetTitle('PV resolution in y [#mum]')
        elif 'reso_z' in f: h0.GetYaxis().SetTitle('PV resolution in z [#mum]')
        elif 'pull_x' in f: h0.GetYaxis().SetTitle('PV pull in x')
        elif 'pull_y' in f: h0.GetYaxis().SetTitle('PV pull in y')
        elif 'pull_z' in f: h0.GetYaxis().SetTitle('PV pull in z')
        
    leg.Draw()

    t1, t2, t3, t4 = style.cmslabel(1, c.year, '', False)
    t1.Draw()
    t2.Draw()
    t3.Draw()
    t4.Draw()
    
    # fit

    if options.param == 'pt' and options.fit:

        ptmin_qcd = 0.4
        ptmin_zb = 0.2
        
        gr, fres = {}, {}
        
        for i, f in enumerate(files):

            x = array.array('d')
            y = array.array('d')
            for b in range(h[f].GetXaxis().GetNbins()):
                x.append(h[f].GetXaxis().GetBinCenter(b+1))
                y.append(h[f].GetBinContent(b+1))
        
            gr[f] = ROOT.TGraph(int(len(x)),x,y)
            
            if 'qcd' in f:
                res = fit.doFitIP('ipfit', gr[f], h[f].GetMarkerColor(), ptmin_qcd, xmax=h0.GetXaxis().GetXmax())
            else:
                res = fit.doFitIP('ipfit', gr[f], h[f].GetMarkerColor(), ptmin_zb, xmax=h[f].GetXaxis().GetXmax())
                
            res.Draw('same')
            
            a = res.GetParameter(0)
            aError = res.GetParError(0)
            b = res.GetParameter(1)
            bError = res.GetParError(1)
            
            fres[f] = ROOT.TLatex(0.62,0.52-i*0.10,'#splitline{a = %.2f #pm %.2f #mum}{b = %.2f #pm %.2f #mum}' %(a, aError, b, bError))
            fres[f].SetNDC()
            fres[f].SetTextAlign(13)
            fres[f].SetTextFont(42)
            fres[f].SetTextSize(0.035)
            fres[f].SetTextColor(h[f].GetMarkerColor())
            fres[f].Draw()

        tex = ROOT.TLatex(0.22,0.35,'#sigma_{d_{0}}(p_{T})  =  #sqrt{a^{2} + #frac{b^{2}}{p_{T}^{2}}}')
        tex.SetNDC()
        tex.SetTextAlign(13)
        tex.SetTextFont(42)
        tex.SetTextSize(0.04)
        tex.Draw()
            
    c1.Print('pics/comp.pdf')
    c1.Clear()
            
