import ROOT
import math
import sys

def doFit(name, hist, func='g1', color=38):
    
    rms = hist.GetRMS()

    xminHist = hist.GetXaxis().GetBinLowEdge(2)
    xmaxHist = hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins()-1)
    
    xmin = -rms*3
    xmax = rms*3
    
    if xmin <= xminHist: xmin = xminHist
    if xmax >= xmaxHist: xmax = xmaxHist
    
#    a = hist.GetMaximum()*math.sqrt(2*3.14)*rms
    a = hist.GetMaximum()
    
    if func == 'g1':
        ffunc = ROOT.TF1(name,"[0]/sqrt(2*pi)/[1]*exp(-x*x/2/[1]/[1])",xmin,xmax)
    elif func == 'g2':
        ffunc = ROOT.TF1(name,"[0]/sqrt(2*pi)/[1]*exp(-x*x/2/[1]/[1])+[2]/sqrt(2*pi)/[3]*exp(-x*x/2/[3]/[3])",xmin,xmax)
#        ffunc.SetParLimits(0,a*(1-0.5),a*(1+0.5))
#        ffunc.SetParLimits(2,a*0.1*(1-0.5),a*0.1*(1+0.5))
    elif func == 'g3':
        ffunc = ROOT.TF1(name,"[0]/sqrt(2*pi)/[1]*exp(-x*x/2/[1]/[1])+[2]/sqrt(2*pi)/[3]*exp(-x*x/2/[3]/[3])+[4]/sqrt(2*pi)/[5]*exp(-x*x/2/[5]/[5])",xmin,xmax)
    else:
        print 'Function is unknown'
        sys.exit()
        
    ffunc.SetParameter(0,a)
    ffunc.SetParameter(1,rms)
    if func == 'g2' or func == 'g3':
        ffunc.SetParameter(2,a*0.05)
        ffunc.SetParameter(3,rms*3.)
    if func == 'g3':
        ffunc.SetParameter(4,a*0.05)
        ffunc.SetParameter(5,rms*5.)        
        
    hist.Fit(name,"QR")
    
    res = hist.GetFunction(name)
    
    chi2 = res.GetChisquare()/res.GetNDF()
    
    res.SetLineColor(color)

#    reso = res.GetParameter(1)
#    resoErr = res.GetParError(1)
    reso = 10E+10
    resoErr = 10E+10
    for p in range(res.GetNpar()):
        if res.GetParName(p) in ['p1','p3','p5']:
            par = res.GetParameter(p)
            err = res.GetParError(p)
            if par < reso and par >= 0:
#            if par < reso and par >= 0:
                reso = par
                resoErr = err
    
    return res, reso, resoErr, chi2
