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
    parser.add_option("--type", default='pv,bs,pv,bs', help="type of IP measurement (pv or bs) [default: %default]")
    parser.add_option("--pkl", default="reso_d0_zb_data_deconv,reso_d0_zb_data_deconv,reso_d0_zb_mc_deconv,reso_d0_zb_mc_deconv", help="list of measurements [default: %default]")
    parser.add_option("--fit", action='store_true', help="Perform a pT-fit [default: %default]")
    parser.add_option("--selection", action='store_true', help="Multi parametrisation draw settings [default: %default]")
    parser.add_option("--log", action='store_true', help="Use log scale on y axis [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()

    ROOT.gROOT.SetBatch()

    pstyle = style.SetPlotStyle(2)

    files = options.pkl.split(',')
    types = options.type.split(',')
    
    c1 = ROOT.TCanvas()
    
#    if 'pt' in options.param: c1.SetLogx(1)

    h = {}

    leg = ROOT.TLegend(0.52,0.72,0.85,0.55)
    if options.param in ['phi', 'npv']:
#        leg = ROOT.TLegend(0.52,0.42,0.85,0.25)
        leg = ROOT.TLegend(0.22,0.82,0.55,0.65)
#    if options.selection:
#        leg = ROOT.TLegend(0.42,0.72,0.75,0.55)
    if options.selection:
        if 'pt' in options.param:
            leg = ROOT.TLegend(0.42,0.72,0.65,0.55)
        else:
            leg = ROOT.TLegend(0.32,0.72,0.75,0.55)
    leg.SetFillColor(253)
    leg.SetBorderSize(0)
    
    for i, f in enumerate(files):
        
        fname = f+'_'+types[i]
        
        print '  ->', fname

        pref = 'pv_'
        if options.measurement == 'ip': pref = options.measurement+types[i]+options.param+'_'
            
        h[fname] = pickle.load(open('results/'+pref+f+'.pkl','rb'))

        h[fname].SetMarkerStyle(26)
        if 'dz' in f and not options.selection:
            h[fname].SetMarkerStyle(22)
            
        h[fname].SetMarkerSize(0.7)
        
        if 'qcd' in fname:
            h[fname].SetMarkerStyle(24)
            if 'dz' in f and not options.selection:
                h[fname].SetMarkerStyle(20)
            h[fname].SetMarkerColor(2)
            h[fname].SetLineColor(2)
        
        if 'mc' in fname:
            h[fname].SetMarkerColor(9)
            h[fname].SetLineColor(9)

        if 'bs' in types[i]:
            if 'qcd' in f: h[fname].SetMarkerStyle(20)
            else: h[fname].SetMarkerStyle(22)
            
        if options.selection:
            mstydata = [20, 21, 22, 23]
            mstymc = [24, 25, 26, 27]
            mcol = [2, 4, 6, 8] 
            if i % 2 == 0:
                h[fname].SetMarkerStyle(mstydata[int(i/2)])
            else:
                h[fname].SetMarkerStyle(mstymc[int(i/2)])
            h[fname].SetMarkerColor(mcol[int(i/2)])
            h[fname].SetLineColor(mcol[int(i/2)])

        lab = ''
        
#        if 'zb' in fname: lab += 'ZeroBias '+types[i].upper()+' ('
#        else: lab += ' QCD '+types[i].upper()+' ('
        if options.measurement == 'pv': lab += types[i].upper()+' ('
        else: lab += options.measurement.upper()+types[i].upper()+' ('

        if 'data' in fname: lab += 'Data)'
        else: lab += 'Simulation)'
        
        if options.selection:
            
            cvar = 'pt'
            if options.param == 'pt': cvar = 'eta'

            cuts = fname.split('_')[2].replace(cvar,'').replace('p','.').replace('abs','').split('to')
            
            if len(fname.split('_')) > 6:
                cmin = cuts[0]
                cmax = cuts[1]
            else:
                cmin = '0.0'
                cmax = '2.5'
            
            addPar = False
            if len(fname.split('_')) > 7:
                cutseta = fname.split('_')[3].replace('eta','').replace('p','.').replace('abs','').split('to')
                addPar = True
            
            if options.param == 'pt':
                lab = cmin + ' < |#eta| < ' + cmax
            else:
                if cmin == '0.0': cmin = '0.1'
                lab = cmin + ' < p_{T} < ' + cmax + ' GeV'
                if addPar: lab += ', ' + cutseta[0] + ' < |#eta| < ' + cutseta[1]
                
#            if 'mc' in f: lab += ' (Simulation)'
#            else: lab += ' (Data)'
        
        leg.AddEntry(h[fname], lab, 'p')
        
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
       
    if options.measurement == 'ip' and options.param == 'pt': xmin = min(xmin, 0.1)
    if 'pt' in options.param: 
        xmin = xmin/2.
        xmax = xmax*1.05
#        xmax = xmax*2.0
        
    ymin = min(ymin, 0.)
    if options.param not in ['phi', 'npv']:
        ymax = ymax*1.2
    elif options.param == 'phi': ymax = 250.
    else: 
        ymax = 180.
        if options.param == 'npv' and 'dz' in f:
            ymax *= 1.2
    
    if options.param == 'eta':
        if options.selection and 'd0' in f: ymax *= 1.2
        else:
            if options.selection: ymin = 10.
            else: ymin = 20.
            if 'dz' in f:
                ymax *= 10.
        
    # temporary
    if options.param == 'pt':
        if 'zb' in f: ymax = 180.
        if 'qcd' in f: ymax = 180.
    if options.param == 'eta':
        if 'zb' in f: ymax = 150.
        if 'qcd' in f: ymax = 150.

#    if 'npv' in options.param: ymax *= 1.5
    
    h0 = c1.DrawFrame(xmin, ymin, xmax, ymax)    

    if 'pt' in options.param: h0.GetXaxis().SetTitle('Track p_{T} [GeV]')
    elif 'eta' in options.param: h0.GetXaxis().SetTitle('Track #eta')
    elif 'phi' in options.param: h0.GetXaxis().SetTitle('Track #phi')
    elif 'npv' in options.param: h0.GetXaxis().SetTitle('Number of primary vertices')
    elif 'dr' in options.param: h0.GetXaxis().SetTitle('#DeltaR(jet axis,track)')
    
    if options.measurement == 'pv': h0.GetXaxis().SetTitle('#sqrt{#sum p_{T}^{2}} [GeV]')
    
    for i, f in enumerate(files):

        fname = f+'_'+types[i]
        
        h[fname].Draw('pe1 same')
        
        if 'd0' in f: h0.GetYaxis().SetTitle('Track IP resolution (d_{xy}) [#mum]')
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

    if 'eta' in options.param and 'dz' in files[0]: c1.SetLogy(1)
    if 'pt' in options.param:
        if options.log:
            h0.SetMinimum(5.)
            h0.SetMaximum(1000.)
            c1.SetLogy(1)
    
#    if options.selection:

    if options.param in ['phi', 'npv']:
        if 'zb' in files[0]:
            lab = ROOT.TLatex(0.70,0.70,'ZeroBias')
        else:
            lab = ROOT.TLatex(0.70,0.70,'QCD')
    else:
        if 'zb' in files[0]:
            lab = ROOT.TLatex(0.30,0.80,'ZeroBias')
        else:
            lab = ROOT.TLatex(0.30,0.80,'QCD')
    lab.SetTextSize(0.045)
    lab.SetNDC()
    lab.Draw()

    # fit

    if options.param == 'pt' and options.fit:

        ptmin_qcd = 0.4
        ptmin_zb = 0.2
        
        gr, fres = {}, {}
        
        for i, f in enumerate(files):

            fname = f+'_'+types[i]
            
            x = array.array('d')
            y = array.array('d')
            for b in range(h[fname].GetXaxis().GetNbins()):
                x.append(h[fname].GetXaxis().GetBinCenter(b+1))
                y.append(h[fname].GetBinContent(b+1))
        
            gr[fname] = ROOT.TGraph(int(len(x)),x,y)
            
            if 'qcd' in fname:
                res = fit.doFitIP('ipfit', gr[fname], h[fname].GetMarkerColor(), ptmin_qcd, xmax=h0.GetXaxis().GetXmax())
            else:
                res = fit.doFitIP('ipfit', gr[fname], h[fname].GetMarkerColor(), ptmin_zb, xmax=h[fname].GetXaxis().GetXmax())
                
            res.Draw('same')
            
            a = res.GetParameter(0)
            aError = res.GetParError(0)
            b = res.GetParameter(1)
            bError = res.GetParError(1)
            
            fres[fname] = ROOT.TLatex(0.62,0.52-i*0.10,'#splitline{a = %.2f #pm %.2f #mum}{b = %.2f #pm %.2f #mum}' %(a, aError, b, bError))
            fres[fname].SetNDC()
            fres[fname].SetTextAlign(13)
            fres[fname].SetTextFont(42)
            fres[fname].SetTextSize(0.035)
            fres[fname].SetTextColor(h[fname].GetMarkerColor())
            fres[fname].Draw()

        tex = ROOT.TLatex(0.22,0.35,'#sigma_{d_{0}}(p_{T})  =  #sqrt{a^{2} + #frac{b^{2}}{p_{T}^{2}}}')
        tex.SetNDC()
        tex.SetTextAlign(13)
        tex.SetTextFont(42)
        tex.SetTextSize(0.04)
        tex.Draw()
            
    c1.Print('pics/comp.pdf')
    c1.Print('pub/comp.pdf')
    c1.SaveAs('pub/comp.root')
    c1.Clear()
            
