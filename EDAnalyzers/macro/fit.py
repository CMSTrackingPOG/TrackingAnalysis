import ROOT
import math
import sys

def setFunction(gfit, fname, xmin, xmax):
    
    if gfit == '1g':
        ffunc = ROOT.TF1(fname,"[0]/sqrt(2*pi)/[1]*exp(-x*x/2/[1]/[1])",xmin,xmax)
    elif gfit == '2g':
        ffunc = ROOT.TF1(fname,"[0]*exp(-x*x/2/[1]/[1])+[2]*exp(-x*x/2/[3]/[3])",xmin,xmax)
    elif gfit == '3g':
        ffunc = ROOT.TF1(fname,"[0]*exp(-x*x/2/[1]/[1])+[2]*exp(-x*x/2/[3]/[3])+[4]*exp(-x*x/2/[5]/[5])",xmin,xmax)
    else:
        print 'Function is unknown'
        sys.exit()
        
    return ffunc

def doFit(name, hist, var, param, color=38, func=''):

    rand = ROOT.TRandom3()
    
    rms = hist.GetRMS()

    xminHist = hist.GetXaxis().GetBinLowEdge(2)
    xmaxHist = hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins()-1)
    
    xmin = -rms*3
    xmax = rms*3
    
    if xmin <= xminHist: xmin = xminHist
    if xmax >= xmaxHist: xmax = xmaxHist
    
    a = hist.GetMaximum()

    nTries = 100
    chi2Min = 10E+10
    tryMin = -1
    fMin = ''

    p0, p1, p2, p3, p4, p5 = {}, {}, {}, {}, {}, {}
    p0Err, p1Err, p2Err, p3Err, p4Err, p5Err = {}, {}, {}, {}, {}, {}
        
    for f in ['3g','2g','1g']:

        if func != '' and func != f: continue
        
        ffunc = setFunction(f, name, xmin, xmax)
            
        ffunc.SetParameter(0,a)
        ffunc.SetParameter(1,0.5*rms)
        if f == '2g':
            ffunc.SetParameter(2,a*0.2)
            ffunc.SetParameter(3,rms*5.)
        elif f == '3g':
            ffunc.SetParameter(0,a)
            ffunc.SetParameter(1,0.5*rms)
            ffunc.SetParameter(2,0.2*a)
            ffunc.SetParameter(3,5.*rms)
            ffunc.SetParameter(4,0.1*a)
            ffunc.SetParameter(5,30.*rms)
            
        hist.Fit(name,"QR")
        res = hist.GetFunction(name)
        if res.GetNDF() < 1: continue
        chi2 = res.GetChisquare()/res.GetNDF()
        
        if chi2 > 1:

            p0[f], p1[f], p2[f], p3[f], p4[f], p5[f] = ({} for i in range(6))
            p0Err[f], p1Err[f], p2Err[f], p3Err[f], p4Err[f], p5Err[f] = ({} for i in range(6))
            
            for it in range(nTries):
                
                ra1 = rand.Uniform(0.,1.)
                rs1 = rand.Uniform(0.,1.)
                ffunc.SetParameter(0,a*ra1)
                ffunc.SetParameter(1,rms*rs1)
                if f in ['2g']:
                    ra1 = rand.Uniform(0.,1.)
                    rs1 = rand.Uniform(0.,1.)
                    ffunc.SetParameter(0,a*ra1*0.7)
                    ffunc.SetParameter(1,rms*rs1)
                    ra2 = rand.Uniform(0.,1.)
                    rs2 = rand.Uniform(0.,1.)
                    ffunc.SetParameter(2,a*ra2*0.3)
                    ffunc.SetParameter(3,rms*rs2*5.)
                if f in ['3g']:
                    ra1 = rand.Uniform(0.,1.)
                    rs1 = rand.Uniform(0.,1.)
                    ffunc.SetParameter(0,a*ra1*0.5)
                    ffunc.SetParameter(1,rms*rs1)
                    ra2 = rand.Uniform(0.,1.)
                    rs2 = rand.Uniform(0.,1.)
                    ffunc.SetParameter(2,a*ra2*0.3)
                    ffunc.SetParameter(3,rms*rs2*5.)
                    ra3 = rand.Uniform(0.,1.)
                    rs3 = rand.Uniform(0.,1.)
                    ffunc.SetParameter(4,a*ra3*0.2)
                    ffunc.SetParameter(5,rms*rs3*30.)                    
                
                hist.Fit(name,"QR")
                
                if param == '': break
                
                res = hist.GetFunction(name)
                chi2 = res.GetChisquare()/res.GetNDF()
                
                if chi2 < chi2Min: 
                                
                    p0[f][it] = ffunc.GetParameter(0)
                    p0Err[f][it] = ffunc.GetParError(0)
                    p1[f][it] = ffunc.GetParameter(1)
                    p1Err[f][it] = ffunc.GetParError(1)
                    
                    if f in ['2g','3g']:
                        if p0Err[f][it] > p0[f][it] or p1Err[f][it] > p1[f][it]: continue
                    
                    if f in ['2g','3g']:
                        p2[f][it] = ffunc.GetParameter(2)
                        p2Err[f][it] = ffunc.GetParError(2)
                        p3[f][it] = ffunc.GetParameter(3)
                        p3Err[f][it] = ffunc.GetParError(3)

                        if p2Err[f][it] > p2[f][it] or p3Err[f][it] > p3[f][it]: continue
                        
                    if f in ['3g']:
                        p4[f][it] = ffunc.GetParameter(4)
                        p4Err[f][it] = ffunc.GetParError(4)
                        p5[f][it] = ffunc.GetParameter(5)
                        p5Err[f][it] = ffunc.GetParError(5)
                        
                        if p4Err[f][it] > p4[f][it] or p5Err[f][it] > p5[f][it]: continue

                    chi2Min = chi2
                    tryMin = it
                    fMin = f

                    if chi2 < 1: break
    
    if tryMin >= 0:
        
        ffunc = setFunction(fMin, name, xmin, xmax)
        
        ffunc.SetParameter(0,p0[fMin][tryMin])
        ffunc.SetParameter(1,p1[fMin][tryMin])
        if fMin in ['2g','3g']:
            ffunc.SetParameter(2,p2[fMin][tryMin])
            ffunc.SetParameter(3,p3[fMin][tryMin])
        if fMin in ['3g']:
            ffunc.SetParameter(4,p4[fMin][tryMin])
            ffunc.SetParameter(5,p5[fMin][tryMin])
        hist.Fit(name,"QR")
            
    res = hist.GetFunction(name)
        
    chi2 = res.GetChisquare()/res.GetNDF()
    
    res.SetLineColor(color)
    
    ymax = res.GetMaximum()
    xmax = res.GetMaximumX()
    xlow = res.GetXmin()
    xhigh = res.GetXmax()
    wl = res.GetX(ymax/2.,xlow,xmax)
    wr = res.GetX(ymax/2.,xmax,xhigh)
    FWHM = wr - wl
    reso = FWHM/2.36
    resoErr = 0.

    return res, reso, resoErr, chi2

def doFitIP(name, gr, color=38):

    ffunc = ROOT.TF1(name,"sqrt([0]*[0]+[1]*[1]/(x*x))",0.5,10.)
        
    ffunc.SetParameter(0,40)
    ffunc.SetParameter(1,80)

    gr.Fit(name,"QR")

    res = gr.GetFunction(name)
    chi2 = res.GetChisquare()/res.GetNDF()

    p0 = ffunc.GetParameter(0)
    p1 = ffunc.GetParameter(1)
    
    res.SetLineColor(color)
    
    return res, p0, p1, chi2
