import ROOT
import math
import sys

class fitFunction():

    def __init__(self, gfit, fname, xmin, xmax, par = None, chi2 = None, ndof = None, parFit = None, errFit = None):

        self.gfit = gfit
        
        if gfit == '1g':
            self.f = ROOT.TF1(fname, "[0]/sqrt(2*pi)/[1]*exp(-x*x/2/[1]/[1])", xmin, xmax)
        elif gfit == '2g':
            self.f = ROOT.TF1(fname, "[0]*exp(-x*x/2/[1]/[1])+[2]*exp(-x*x/2/[3]/[3])", xmin, xmax)
        elif gfit == '3g':
            self.f = ROOT.TF1(fname, "[0]*exp(-x*x/2/[1]/[1])+[2]*exp(-x*x/2/[3]/[3])+[4]*exp(-x*x/2/[5]/[5])", xmin, xmax)
        else:
            print 'Function is unknown'
            sys.exit()
            
        self.par = par
        self.parFit = parFit
        self.errFit = errFit
        self.chi2 = chi2
        self.ndof = ndof
        
def doFit(name, hist, var, param, color=38, func=''):
    
    rand = ROOT.TRandom3()
    
    rms = hist.GetRMS()

    xminHist = hist.GetXaxis().GetBinLowEdge(2)
    xmaxHist = hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins()-1)
    
    xmin = -rms*4
    xmax = rms*4
    
    if xmin <= xminHist: xmin = xminHist
    if xmax >= xmaxHist: xmax = xmaxHist
    
    a = float(hist.GetMaximum())

    nTries = 3
    chi2Min = 10E+10
    fMin = []
    
    nev = float(hist.GetEntries())

    for f in ['3g','2g','1g']:

        if func != '' and func != f: continue
        
        for it in range(nTries):

            fname = name+'_'+f+'_'+str(it)
            
            ffunc = fitFunction(f, fname, xmin, xmax)
                
            ra1 = rand.Uniform(0.,1.)
            rs1 = rand.Uniform(0.,1.)
            ffunc.f.SetParameter(0,a*ra1)            
            ffunc.f.SetParameter(1,rms*rs1)
            
            p0i = ffunc.f.GetParameter(0)
            p1i = ffunc.f.GetParameter(1)

            if f in ['2g']:
                
                ra1 = rand.Uniform(0.,1.)
                rs1 = rand.Uniform(0.,1.)
                ffunc.f.SetParameter(0,a*ra1*0.7)
                ffunc.f.SetParameter(1,rms*rs1)
                ra2 = rand.Uniform(0.,1.)
                rs2 = rand.Uniform(0.,1.)
                ffunc.f.SetParameter(2,a*ra2*0.3)
                ffunc.f.SetParameter(3,rms*rs2*5.)
                
                p2i = ffunc.f.GetParameter(2)
                p3i = ffunc.f.GetParameter(3)
            
            if f in ['3g']:
                
                ra1 = rand.Uniform(0.,1.)
                rs1 = rand.Uniform(0.,1.)
                ffunc.f.SetParameter(0,a*ra1*0.5)
                ffunc.f.SetParameter(1,rms*rs1)
                ra2 = rand.Uniform(0.,1.)
                rs2 = rand.Uniform(0.,1.)
                ffunc.f.SetParameter(2,a*ra2*0.3)
                ffunc.f.SetParameter(3,rms*rs2*5.)
                ra3 = rand.Uniform(0.,1.)
                rs3 = rand.Uniform(0.,1.)
                ffunc.f.SetParameter(4,a*ra3*0.2)
                ffunc.f.SetParameter(5,rms*rs3*30.)                    

                p2i = ffunc.f.GetParameter(2)
                p3i = ffunc.f.GetParameter(3)
                p4i = ffunc.f.GetParameter(4)
                p5i = ffunc.f.GetParameter(5)

            hist.Fit(fname, 'QR0+')

            if param == '': break

            res = hist.GetFunction(fname)
            
            chi2 = res.GetChisquare()
            ndof = res.GetNDF()
                
            if chi2/ndof < chi2Min:
                    
                par, parFit, errFit = [], [], []

                unc = 1./math.sqrt(nev)
                frac = 0.01

                p0 = res.GetParameter(0); p0err = res.GetParError(0)
                p1 = res.GetParameter(1); p1err = res.GetParError(1)
                p0rUnc = abs(p0err/p0) if p0 != 0 else 1.; p1rUnc = abs(p1err/p1) if p1 != 0 else 1.
                if any(x < 0 for x in [p0, p0err, p1, p1err]) or (p0rUnc > unc*p0) or (p1rUnc > unc*p1): continue
                if p0/a < frac: continue
                par.append(p0i); parFit.append(p0); errFit.append(p0err)
                par.append(p1i); parFit.append(p1); errFit.append(p1err)
                
                if f in ['2g', '3g']:
                    p2 = res.GetParameter(2); p2err = res.GetParError(2)
                    p3 = res.GetParameter(3); p3err = res.GetParError(3)
                    p2rUnc = abs(p2err/p2) if p2 != 0 else 1.; p3rUnc = abs(p3err/p3) if p3 != 0 else 1.
                    if any(x < 0 for x in [p2, p2err, p3, p3err]) or (p2rUnc > unc*p2) or (p3rUnc > unc*p3): continue
                    if p2/a < frac: continue
                    par.append(p2i); parFit.append(p2); errFit.append(p2err)
                    par.append(p3i); parFit.append(p3); errFit.append(p3err)

                if f in ['3g']:
                    p4 = res.GetParameter(4); p4err = res.GetParError(4)
                    p5 = res.GetParameter(5); p5err = res.GetParError(5)
                    p4rUnc = abs(p4err/p4) if p4 != 0 else 1.; p5rUnc = abs(p5err/p5) if p5 != 0 else 1.
                    if any(x < 0 for x in [p4, p4err, p5, p5err]) or (p4rUnc > unc*p4) or (p5rUnc > unc*p5): continue
                    if p4/a < frac: continue
                    par.append(p4i); parFit.append(p4); errFit.append(p4err)
                    par.append(p5i); parFit.append(p5); errFit.append(p5err)

                chi2Min = chi2/ndof

                fres = fitFunction(f, fname+'_sel', xmin, xmax, par=par, chi2=chi2, ndof=ndof, parFit=parFit, errFit=errFit)
                fMin.append(fres)

                if chi2Min < 1.: break

    if len(fMin) == 0:
        print 'Failed to find a suitable fit function'
        sys.exit()
 
    fMin.sort(key=lambda x: x.chi2/x.ndof, reverse=False)

    f = fitFunction(fMin[0].gfit, 'finalFit', xmin, xmax)
    for ip in range(len(fMin[0].parFit)):
        f.f.SetParameter(ip, fMin[0].parFit[ip])

    hist.Fit('finalFit', 'QR')
    res = hist.GetFunction('finalFit')
    
    chi2 = res.GetChisquare()/res.GetNDF()

    res.SetLineColor(color)
    
    xmax = res.GetMaximumX()
    xlow = res.GetXmin()
    xhigh = res.GetXmax()
    ymax = res.GetMaximum(xlow, xhigh)
    wl = res.GetX(ymax/2., xlow, xmax)
    wr = res.GetX(ymax/2., xmax, xhigh)

    FWHM = wr - wl
    reso = FWHM/2.36
    resoErr = reso/math.sqrt(2.*nev)

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
