import ROOT
import math
import sys

class fitFunction():

    def __init__(self, gfit, fname, xmin, xmax, par = None, chi2 = None, ndof = None, parFit = None, errFit = None):

        self.gfit = gfit
        self.par = par
        self.parFit = parFit
        self.errFit = errFit
        self.chi2 = chi2
        self.ndof = ndof
        
def doFit(name, hist, var, param, color=38, func='', nsig=4, nTries=3):

    formula = {}
    formula['1g'] = "[0]/sqrt(2*pi)/[1]*exp(-x*x/2/[1]/[1])"
    formula['2g'] = "[0]*exp(-x*x/2/[1]/[1])+[2]*exp(-x*x/2/[3]/[3])"
    formula['3g'] = "[0]*exp(-x*x/2/[1]/[1])+[2]*exp(-x*x/2/[3]/[3])+[4]*exp(-x*x/2/[5]/[5])"
    
    rand = ROOT.TRandom3()
    
    rms = hist.GetRMS()

    xminHist = hist.GetXaxis().GetBinLowEdge(2)
    xmaxHist = hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins()-1)
    
    xmin = -rms*nsig
    xmax = rms*nsig
    
    if xmin <= xminHist: xmin = xminHist
    if xmax >= xmaxHist: xmax = xmaxHist
    
    a = float(hist.GetMaximum())

    bmin = hist.GetXaxis().FindBin(xmin)
    bmax = hist.GetXaxis().FindBin(xmax)
    nbfit = bmax-bmin

    nbfitmin = 10
    while nbfit < nbfitmin:
        nsig += 0.1
        xmin = -rms*nsig
        xmax = rms*nsig

        if xmin <= xminHist or xmax >= xmaxHist:
            print 'Please extend the initial x axis range in the displayed histogram'
            sys.exit()
        
        bmin = hist.GetXaxis().FindBin(xmin)
        bmax = hist.GetXaxis().FindBin(xmax)
        nbfit = bmax-bmin

    binw = hist.GetXaxis().GetBinWidth(1)
        
    chi2Min = 10E+10
    fMin = []
    
    nev = float(hist.GetEntries())

    for f in ['3g','2g','1g']:

        if func != '' and func != f: continue

        for it in range(nTries):
            
            fname = name+'_'+f+'_'+str(it)

            f1 = ROOT.TF1(fname, formula[f], xmin, xmax)
            
            if f in ['1g', '2g', '3g']:
                ra1 = rand.Uniform(0.,1.)
                rs1 = rand.Uniform(0.,1.)
                f1.SetParameter(0, a*ra1)
                f1.SetParameter(1, binw+rms*rs1)
                f1.SetParLimits(0, 0., 3*a)
                f1.SetParLimits(1, binw, 1000*rms)
                p0i = f1.GetParameter(0)
                p1i = f1.GetParameter(1)
                
            if f in ['2g', '3g']:
                ra1 = rand.Uniform(0.,1.)
                rs1 = rand.Uniform(0.,1.)
                f1.SetParameter(2, a*ra1)
                f1.SetParameter(3, binw+rms*rs1)
                f1.SetParLimits(2, 0., 3*a)
                f1.SetParLimits(3, binw, 1000*rms)
                p2i = f1.GetParameter(2)
                p3i = f1.GetParameter(3)
                
            if f in ['3g']:
                ra1 = rand.Uniform(0.,1.)
                rs1 = rand.Uniform(0.,1.)
                f1.SetParameter(4, a*ra1)
                f1.SetParameter(5, binw+rms*rs1)
                f1.SetParLimits(4, 0., 3*a)
                f1.SetParLimits(5, binw, 1000*rms)
                p4i = f1.GetParameter(4)
                p5i = f1.GetParameter(5)
                
            ress = hist.Fit(f1, 'QRS0+')

            if param == '': break

            res = hist.GetFunction(fname)
            
            chi2 = res.GetChisquare()
            ndof = res.GetNDF()

            if chi2/ndof < chi2Min:
                
                par, parFit, errFit = [], [], []

                p0 = res.GetParameter(0); p0err = res.GetParError(0)
                p1 = res.GetParameter(1); p1err = res.GetParError(1)

                par.append(p0i); parFit.append(p0); errFit.append(p0err)
                par.append(p1i); parFit.append(p1); errFit.append(p1err)
                
                if f in ['2g', '3g']:
                    p2 = res.GetParameter(2); p2err = res.GetParError(2)
                    p3 = res.GetParameter(3); p3err = res.GetParError(3)                    

                    par.append(p2i); parFit.append(p2); errFit.append(p2err)
                    par.append(p3i); parFit.append(p3); errFit.append(p3err)

                if f in ['3g']:
                    p4 = res.GetParameter(4); p4err = res.GetParError(4)
                    p5 = res.GetParameter(5); p5err = res.GetParError(5)

                    par.append(p4i); parFit.append(p4); errFit.append(p4err)
                    par.append(p5i); parFit.append(p5); errFit.append(p5err)

                chi2Min = chi2/ndof

                fres = fitFunction(f, fname+'_sel', xmin, xmax, par=par, chi2=chi2, ndof=ndof, parFit=parFit, errFit=errFit)
                fMin.append(fres)

    if len(fMin) == 0:
        
        print 'Failed to find a suitable fit function'
        
        return None, -1, -1, 10e+10
 
    fMin.sort(key=lambda x: x.chi2/x.ndof, reverse=False)

    f1 = ROOT.TF1('finalFit', formula[fMin[0].gfit], xmin, xmax)
    for ip in range(len(fMin[0].par)):
        pfit = fMin[0].parFit[ip]
        f1.SetParameter(ip, pfit)
#        print ip, pfit
        if (ip % 2 == 0): f1.SetParLimits(ip, 0., 3*a)
        else: f1.SetParLimits(ip, 0., 1000*rms)

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

def doFitIP(name, gr, color=38, xmin=0.0, xmax=10.0, p0s = 40, p1s = 80):

    ffunc = ROOT.TF1(name, "sqrt([0]*[0]+[1]*[1]/(x*x))", xmin, xmax)
        
    ffunc.SetParameter(0, p0s)
    ffunc.SetParameter(1, p1s)

    gr.Fit(name, "QR")

    res = gr.GetFunction(name)

    res.SetLineColor(color)
    
    return res

def fwhm(h, ffit = None, nmin = 10000):
    
    nev = float(h.GetEntries())

    bscanmin = h.GetXaxis().FindBin(h.GetXaxis().GetXmin())
    bscanmax = h.GetXaxis().FindBin(h.GetXaxis().GetXmax())
    
    if ffit is None:
        
        max = h.GetMaximum()/2.
    
        b1, b2 = -1, -1
    
        for b in range(bscanmin, bscanmax+1):
            if h.GetBinContent(b) > max:
                b1 = b
                for bb in range(b1+1, bscanmax+1):
                    if h.GetBinContent(bb) < max:
                        b2 = bb
                        break
            if b2 >= 0: break
    
            #    b1 = h.FindFirstBinAbove(h.GetMaximum()/2.)
            #    b2 = h.FindLastBinAbove(h.GetMaximum()/2.)

        if (b1 < 0) or (b2 < 0) or (b1 == b2):
            print 'FWHM: Use finer binning and more stats, bins =', h.GetXaxis().GetNbins(), ' max =', max
            b1 = h.GetMaximumBin()-1
            b2 = h.GetMaximumBin()+1
                
        reso = (h.GetBinCenter(b2) - h.GetBinCenter(b1))/2.    
        resoErr = reso/math.sqrt(2.*nev)
        sys = ax.GetBinWidth(2)

        return reso, resoErr, sys
    
    else:
        
        ax = h.GetXaxis()
        
        xmin = ax.GetXmin()
        xmax = ax.GetXmax()
        xmaxfpos = ffit.GetMaximumX(xmin, xmax)
        xmaxf = ffit.GetMaximum(xmin, xmax)/2.
        bmax = ax.FindBin(xmaxfpos)

        lmax = -1
        for ib in range(bscanmin, bmax):
            if h.GetBinContent(ib) > xmaxf:
                lmax = ib
                break

        rmax = -1
        for ib in range(bmax+1, bscanmax+1):
            if h.GetBinContent(ib) < xmaxf:
                rmax = ib-1
                break
            
        if lmax >= 0 and rmax >= 0:

            lmaxdist = h.GetBinCenter(lmax)
            rmaxdist = h.GetBinCenter(rmax)
            
            reso = (rmaxdist-lmaxdist)/2.            
            resoErr = reso/math.sqrt(2.*nev)
            sys = ax.GetBinWidth(bmax)
            
        else: return 0., 0., 0.
            
        if h.Integral() < nmin: return 0., 0., 0.
        
        return reso, resoErr, sys

def rebin(h, nmin, nbinmin):

    hcl = h.Clone('reb')
    
    nbins = hcl.GetXaxis().GetNbins()
    hcl.GetXaxis().SetRange(1, nbins)
    max = hcl.GetMaximum()
    
    fac = 1
    
    while (max < nmin) and (nbins > nbinmin):
        
        facc = -1
        for f in [2, 3]:
            r = nbins % f
            if r == 0:
                facc = f
                break

        if facc < 0: return fac
        else: fac *= facc
        
        hcl = hcl.Rebin(facc)
        nbins = hcl.GetXaxis().GetNbins()
        hcl.GetXaxis().SetRange(1, nbins)
        max = hcl.GetMaximum()
        
    return fac
    
    
