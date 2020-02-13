import ROOT
import math
import sys

def doFit(name, hist, var, param, color=38):

    rand = ROOT.TRandom3()
    
    rms = hist.GetRMS()

    xminHist = hist.GetXaxis().GetBinLowEdge(2)
    xmaxHist = hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins()-1)
    
    xmin = -rms*3
    xmax = rms*3
    
    if xmin <= xminHist: xmin = xminHist
    if xmax >= xmaxHist: xmax = xmaxHist
    
    a = hist.GetMaximum()

    if var in ['x','y','z']:
        ffunc = ROOT.TF1(name,"[0]/sqrt(2*pi)/[1]*exp(-x*x/2/[1]/[1])",xmin,xmax)
    elif var == 'd0':
        ffunc = ROOT.TF1(name,"[0]*exp(-x*x/2/[1]/[1])+[2]*exp(-x*x/2/[3]/[3])",xmin,xmax)
#        ffunc.SetParLimits(0,0.5*a,1.1*a)
#        ffunc.SetParLimits(1,0.,1.5*rms)
#        ffunc.SetParLimits(2,0.,0.3*a)
#        ffunc.SetParLimits(3,rms,4*rms)
    elif var == 'dz':
        ffunc = ROOT.TF1(name,"[0]*exp(-x*x/2/[1]/[1])+[2]*exp(-x*x/2/[3]/[3])+[4]*exp(-x*x/2/[5]/[5])",xmin,xmax)
#        ffunc.SetParLimits(0,0.4*a,1.0*a)
#        ffunc.SetParLimits(1,0.,1.0*rms)
#        ffunc.SetParLimits(2,0.,0.5*a)
#        ffunc.SetParLimits(3,0.5*rms,2*rms)
#        ffunc.SetParLimits(4,0.,0.3*a)
#        ffunc.SetParLimits(5,2*rms,1E+10)
    else:
        print 'Function is unknown'
        sys.exit()
        
    ffunc.SetParameter(0,a)
    ffunc.SetParameter(1,rms)
    if var == 'd0':
        ffunc.SetParameter(2,a*0.01)
        ffunc.SetParameter(3,rms*2.)
    if var == 'dz':
        ffunc.SetParameter(0,a)
        ffunc.SetParameter(1,0.5*rms)
        ffunc.SetParameter(2,0.3*a)
        ffunc.SetParameter(3,1.0*rms)
        ffunc.SetParameter(4,0.1*a)
        ffunc.SetParameter(5,10*rms)
        
    hist.Fit(name,"QR")
    res = hist.GetFunction(name)
    chi2 = res.GetChisquare()/res.GetNDF()
        
    if (var == 'd0' or var == 'dz') and chi2 > 1:
        nTries = 100
        tryMin = -1
        chi2Min = 10E+10
        p0, p1, p2, p3, p4, p5 = ({} for i in range(6))
        for f in range(nTries):
            hist.Fit(name,"QR")
            if param == '': break
            res = hist.GetFunction(name)
            chi2 = res.GetChisquare()/res.GetNDF()

            if chi2 < chi2Min: 
                chi2Min = chi2
                tryMin = f
                
                p0[f] = ffunc.GetParameter(0)
                p1[f] = ffunc.GetParameter(1)
                p2[f] = ffunc.GetParameter(2)
                p3[f] = ffunc.GetParameter(3)
                if var == 'dz':
                    p4[f] = ffunc.GetParameter(4)
                    p5[f] = ffunc.GetParameter(5)
                        
            if var == 'd0':
                ra1 = rand.Uniform(0,1)
                rs1 = rand.Uniform(0,1)
                ra2 = rand.Uniform(0,1)
                rs2 = rand.Uniform(0,1)
                ffunc.SetParameter(0,a*ra1)
                ffunc.SetParameter(1,rms*rs1)
                ffunc.SetParameter(2,a*ra2*0.1)
                ffunc.SetParameter(3,rms*rs2*4.)
            elif var == 'dz':
                ra1 = rand.Uniform(0,1)
                rs1 = rand.Uniform(0,1)
                ra2 = rand.Uniform(0,1)
                rs2 = rand.Uniform(0,1)
                ra3 = rand.Uniform(0,1)
                rs3 = rand.Uniform(0,1)
                ffunc.SetParameter(0,a*ra1)
                ffunc.SetParameter(1,rms*rs1*0.7)
                ffunc.SetParameter(2,a*ra2*0.4)
                ffunc.SetParameter(3,rms*rs2*2.)
                ffunc.SetParameter(4,a*ra3*0.2)
                ffunc.SetParameter(5,rms*rs3*20.)
                
            if chi2 < 1: break
    
        if tryMin >= 0:
            ffunc.SetParameter(0,p0[tryMin])
            ffunc.SetParameter(1,p1[tryMin])
            ffunc.SetParameter(2,p2[tryMin])
            ffunc.SetParameter(3,p3[tryMin])
            if var == 'dz':
                ffunc.SetParameter(4,p4[tryMin])
                ffunc.SetParameter(5,p5[tryMin])
            hist.Fit(name,"QR")
    
    res = hist.GetFunction(name)
    
    chi2 = res.GetChisquare()/res.GetNDF()
    
    res.SetLineColor(color)

    reso = res.GetParameter(1)
    resoErr = res.GetParError(1)
    if var == 'd0' or var == 'dz':
        ymax = res.GetMaximum()
        xmax = res.GetMaximumX()
        xlow = res.GetXmin()
        xhigh = res.GetXmax()
        wl = res.GetX(ymax/2.,xlow,xmax)
        wr = res.GetX(ymax/2.,xmax,xhigh)
        FWHM = wr - wl
        reso = FWHM/2.36
        resoErr = 0.

#    reso = 10E+10
#    resoErr = 10E+10
#    for p in range(res.GetNpar()):
#        if res.GetParName(p) in ['p1','p3','p5']:
#            par = res.GetParameter(p)
#            err = res.GetParError(p)
#            if par < reso and par >= 0:
#                reso = par
#                resoErr = err
    
    return res, reso, resoErr, chi2

def doFitIP(name, gr, color=38):

    ffunc = ROOT.TF1(name,"sqrt([0]*[0]+[1]*[1]/(x*x))",0.,350.)
        
    ffunc.SetParameter(0,24)
    ffunc.SetParameter(1,90)

    gr.Fit(name,"QR")

    res = gr.GetFunction(name)
    chi2 = res.GetChisquare()/res.GetNDF()

    p0 = ffunc.GetParameter(0)
    p1 = ffunc.GetParameter(1)
    
    res.SetLineColor(color)
    
    return res, p0, p1, chi2
