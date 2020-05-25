import ROOT
import numpy
import math
import os, sys
import json
import glob
import time
from os.path import expanduser
from PIL import Image
from array import array
import multiprocessing
import common

def addbin(h):
    
   x_nbins = h.GetXaxis().GetNbins()
   h.SetBinContent(1,h.GetBinContent(0)+h.GetBinContent(1))
   h.SetBinError(1,ROOT.TMath.Sqrt(ROOT.TMath.Power(h.GetBinError(0),2)+ROOT.TMath.Power(h.GetBinError(1),2)))
   h.SetBinContent(x_nbins,h.GetBinContent(x_nbins)+h.GetBinContent(x_nbins+1))
   h.SetBinError(x_nbins,ROOT.TMath.Sqrt(ROOT.TMath.Power(h.GetBinError(x_nbins),2)+ROOT.TMath.Power(h.GetBinError(x_nbins+1),2)))

   h.SetBinContent(0,0.);
   h.SetBinError(0,0.);
   h.SetBinContent(x_nbins+1,0.);
   h.SetBinError(x_nbins+1,0.);

def addbin2D(h):

    x_nbins = h.GetXaxis().GetNbins()
    y_nbins = h.GetYaxis().GetNbins()
            
    for iy in range(1,y_nbins+1):
        
        h.SetBinContent(1,iy,h.GetBinContent(0,iy)+h.GetBinContent(1,iy))
        h.SetBinError(1,iy,math.sqrt(math.pow(h.GetBinError(0,iy),2)+math.pow(h.GetBinError(1,iy),2)))
        h.SetBinContent(x_nbins,iy,h.GetBinContent(x_nbins,iy)+h.GetBinContent(x_nbins+1,iy))
        h.SetBinError(x_nbins,iy,math.sqrt(math.pow(h.GetBinError(x_nbins,iy),2)+math.pow(h.GetBinError(x_nbins+1,iy),2)))

        h.SetBinContent(0,iy,0.)
        h.SetBinError(0,iy,0.)
        h.SetBinContent(x_nbins+1,iy,0.)
        h.SetBinError(x_nbins+1,iy,0.)
         
    for ix in range(1,x_nbins+1):
                          
        h.SetBinContent(ix,1,h.GetBinContent(ix,0)+h.GetBinContent(ix,1))
        h.SetBinError(ix,1,math.sqrt(math.pow(h.GetBinError(ix,0),2)+math.pow(h.GetBinError(ix,1),2)))
        h.SetBinContent(ix,y_nbins,h.GetBinContent(ix,y_nbins)+h.GetBinContent(ix,y_nbins+1))
        h.SetBinError(ix,y_nbins,math.sqrt(math.pow(h.GetBinError(ix,y_nbins),2)+math.pow(h.GetBinError(ix,y_nbins+1),2)))

        h.SetBinContent(ix,0,0.)
        h.SetBinError(ix,0,0.)
        h.SetBinContent(ix,y_nbins+1,0.)
        h.SetBinError(ix,y_nbins+1,0.)
   
def combSysLinear(h_nom,h_sys_up,h_sys_down,h_comb,h_sys_down_comb,h_sys_up_comb,isys):

    sys_up = numpy.empty(1000,dtype=float)
    sys_down = numpy.empty(1000,dtype=float)
    stat = numpy.empty(1000,dtype=float)
    
    nBins = h_nom[0].GetXaxis().GetNbins()
 
    for sys in range(1,nBins+1):
        
        sys_up[sys-1] = 0.
        sys_down[sys-1] = 0.
        stat[sys-1] = 0.

    for ib in range(1,nBins+1):

        delta_down = 0.
	delta_up = 0.

        for i in range(len(h_nom)):
            
	     nom = h_nom[i].GetBinContent(ib)
	     up = h_sys_up[i][isys].GetBinContent(ib)
	     down = h_sys_down[i][isys].GetBinContent(ib)
	     err = h_nom[i].GetBinError(ib)

	     delta_down += (nom-down)
	     delta_up += (up-nom)
	     
	     stat[ib-1] += math.pow(err,2)
	
	sys_up[ib-1] = sys_up[ib-1] + delta_up
	sys_down[ib-1] = sys_down[ib-1] + delta_down

	stat[ib-1] = math.sqrt(stat[ib-1])

    for ib in range(1,nBins+1):
        
        nom = 0.
        
        for i in range(len(h_nom)):
            
            nom += h_nom[i].GetBinContent(ib)

	h_comb.SetBinContent(ib,nom)
	h_comb.SetBinError(ib,stat[ib-1])
	h_sys_down_comb[isys].SetBinContent(ib,nom-sys_down[ib-1])
	h_sys_up_comb[isys].SetBinContent(ib,nom+sys_up[ib-1])

def combSys(h_nom,h_sys_down,h_sys_up,h_sys_down_comb,h_sys_up_comb,nSys):

    sys_up = numpy.empty(1000,dtype=float)
    sys_down = numpy.empty(1000,dtype=float)
   
    nBins = h_nom.GetXaxis().GetNbins()

    for sys in range(1,nBins+1):

        sys_up[sys-1] = 0.
        sys_down[sys-1] = 0.
   
    for j in range(nSys):

        for ib in range(1,nBins+1):

            nom = h_nom.GetBinContent(ib)
            down = h_sys_down[j].GetBinContent(ib)
            up = h_sys_up[j].GetBinContent(ib)
            delta_down = nom-down
            delta_up = up-nom
            sys_up[ib-1] = math.sqrt(math.pow(sys_up[ib-1],2) + math.pow(delta_up,2))
            sys_down[ib-1] = math.sqrt(math.pow(sys_down[ib-1],2) + math.pow(delta_down,2))
   
    for ib in range(1,nBins+1):

        nom = h_nom.GetBinContent(ib)
        err = h_nom.GetBinError(ib)
	
	sys_up[ib-1] = math.sqrt(math.pow(sys_up[ib-1],2) + math.pow(err,2))
	sys_down[ib-1] = math.sqrt(math.pow(sys_down[ib-1],2) + math.pow(err,2))

    for ib in range(1,nBins+1):

        nom = h_nom.GetBinContent(ib)
        h_sys_down_comb.SetBinContent(ib,nom-sys_down[ib-1])
        h_sys_up_comb.SetBinContent(ib,nom+sys_up[ib-1])
        
def totSys(h_nom,h_down,h_up):

   nbins = h_nom.GetXaxis().GetNbins()
   
   for ib in range(1,nbins+1):

       b_nom = h_nom.GetBinContent(ib)
       b_cur_down = h_down.GetBinContent(ib)
       b_cur_up = h_up.GetBinContent(ib)
       
       del_down = b_cur_down - b_nom
       del_up = b_cur_up - b_nom
       del_down_res = del_down
       del_up_res = del_up
       
       del_min = del_down if del_down < del_up else del_up
       del_max = del_down if del_down > del_up else del_up
       
       del_down_res = 0. if del_min > 0. else del_min
       del_up_res = 0. if del_max < 0. else del_max
       
       sys_down = b_nom + del_down_res
       sys_up = b_nom + del_up_res
       
       h_down.SetBinContent(ib,sys_down)
       h_up.SetBinContent(ib,sys_up)

def makeErrorBand(tot,plus,minus):

    nbins = tot.GetNbinsX()

    x = array('f')
    xerr = array('f')
    y = array('f')
    yp = array('f')
    ym = array('f')
    
    for bin in range(1,nbins+1):
        
        index = bin-1
        xerr.append(tot.GetBinWidth(bin)/2.0)
        x.append(tot.GetBinLowEdge(bin) + xerr[index])
        
        y.append(tot.GetBinContent(bin))
        yp.append(plus.GetBinContent(bin)-y[index])
        ym.append(y[index]-minus.GetBinContent(bin))
        if y[index] - ym[index] < 0: ym[index] = y[index]
        
    error = ROOT.TGraphAsymmErrors(nbins, x, y, xerr, xerr, ym, yp)
        
    return error

class pileup():
    
    def __init__(self, dpath, evt):
        
        fName = dpath+'puCentral.root'
        
        if evt == 'zb': fName = dpath+'ZeroBias_puCentral.root'
        else: fName = dpath+'JetHT_puCentral.root'
        
        fData = ROOT.TFile.Open(fName,'READ')
        hData = fData.Get('pileup')
        hData.Scale(1./hData.Integral())
        
        hMC = ROOT.TH1D('pileupMC', 'pileupMC', 100, 0, 100)
        
        sys.path.append(dpath)
        from mix_2017_25ns_UltraLegacy_PoissonOOTPU_cfi import probValue
        for i, value in enumerate(probValue): hMC.SetBinContent(i+1, value)
        
        hMC.Scale(1./hMC.Integral())
        
        ROOT.gROOT.cd()
        self.pileupRW = hData.Clone('pileupRW')
        self.pileupRW.Divide(hMC)
        
    def getWeight(self, nTrue):
        
        return self.pileupRW.GetBinContent(self.pileupRW.FindBin(nTrue))

class reweight():
    
    def __init__(self, dpath, rwvar):
        
        frw = ROOT.TFile.Open(dpath+rwvar+'.root','READ')
        
        hData = frw.Get('h_data')
        hData.Scale(1./hData.Integral())

        hMC = frw.Get('h_mc')
        hMC.Scale(1./hMC.Integral())
        
        ROOT.gROOT.cd()
        self.rw = hData.Clone('varRW')
        self.rw.Divide(hMC)
        
        frw.Close()
        
    def getWeight(self, var):
        
        return self.rw.GetBinContent(self.rw.FindBin(var))

class param():
    
    def __init__(self, fname):

        with open(fname, "r") as read_file:
            self.data = json.load(read_file)

    def get(self, pname):
        
        if pname not in self.data:
            print 'Parametrisation '+pname+' not found'
            sys.exit()
        
        p = self.data[pname]
        
        return p    

def copy(finput, foutput):

    for f in finput:
        os.system('cp '+f+' '+foutput)

def convert(finput, foutput):
    
    size = (400,400)
    
    im = Image.open(finput)
    im.thumbnail(size)
    im.save(foutput)
    
    return foutput
    
def generateThumbs(files, fdir, img):
    
    pool = multiprocessing.Pool(common.ncores)
    
    jobs = []
    output = []
    
    for image in files:
        
        finput = fdir+'/'+image+'.'+img
        foutput = fdir+'/'+image+'_thumb.png'
        
        jobs.append( pool.apply_async(convert, (finput, foutput)) )
        
    output += [job.get() for job in jobs]
    
    pool.close()
        
    return output

def createPage(typ):
    
    pool = multiprocessing.Pool(common.ncores)
    
    home = expanduser("~")
    webpath = home+'/public_html/tracking/'
    
    figs = {}
    
    if typ == 'pv':
        cat = ['pvPull_x', 'pvPull_y', 'pvPull_z', \
        'pvReso_x', 'pvReso_y', 'pvReso_z', \
        'bs', 'Others']
    else:
        cat = [typ+'Reso_d0', typ+'Reso_dz', \
        'Others']
    
    os.system('rm -rf '+webpath)
    os.system('mkdir '+webpath)
    os.system('mkdir '+webpath+'pics/')
    
    page = '<html>\n'    
    page += '<body>\n'    
    
    files = {}
    filesAll = glob.glob('pics/*thumb*')
    
    jobs = []
    
    for c in cat:
        
        files[c] = glob.glob('pics/'+c+'*thumb*')
        
        if len(files[c]) == 0: continue
        
        filesLeft = list(set(filesAll)^set(files[c]))
        filesAll = filesLeft
    
#        page += '<section><h2>'+c+'</h2>\n'
        page += '<a href="'+c+'.html">'+c+'</a>\n'
        
        subpage = '<html>\n'
        subpage += '<body>\n'
        
        for f in files[c]:
        
            finput = [f, f.replace('_thumb','').replace('.png','.eps')]
            foutput = webpath+'pics/'
            
            jobs.append( pool.apply_async(copy, (finput, foutput)) )
            
            subpage += '<a href="'+f.replace('_thumb','').replace('.png','.eps')+'">'+'<img src="'+f+'" height="100"/>'+'</a>\n'

        subpage += '</body>\n'
        subpage += '</html>'

        sfile = open(webpath+c+'.html', 'w')
        sfile.write(subpage)
        sfile.close()
        
#        page += '</section>\n'

    if len(filesAll) > 0:
 
        page += '<a href="others.html">Others</a>\n'
#        page += '<section><h2>Others</h2>\n'

        subpage = '<html>\n'
        subpage += '<body>\n'

        for f in filesAll:
            
            finput = [f, f.replace('_thumb','').replace('.png','.eps')]
            foutput = webpath+'pics/'
            
            jobs.append( pool.apply_async(copy, (finput, foutput)) )
            
            subpage += '<a href="'+f.replace('_thumb','').replace('.png','.eps')+'">'+'<img src="'+f+'" height="100"/>'+'</a>\n'
            
#        page += '</section>\n'

        subpage += '</body>\n'
        subpage += '</html>'

        sfile = open(webpath+'others.html', 'w')
        sfile.write(subpage)
        sfile.close()
        
    page += '</body>\n'
    page += '</html>'

    for job in jobs: job.get()
    
    pool.close()
    
    wfile = open(webpath+'index.html', 'w')
    wfile.write(page)
    wfile.close()

def adjust(h1, h2, nsig=5):
    
    rms = max(h1.GetRMS(), h2.GetRMS())
    
    for h in [h1, h2]: h.GetXaxis().SetRangeUser(-nsig*rms,nsig*rms)

def isbadfit(r1, r2):
    
    tex = ROOT.TLatex(0.3,0.4,'')
    
    rmax = 10.
    
    if (r1 > rmax) or (r2 > rmax):
        tex = ROOT.TLatex(0.3,0.4,'BAD')
        tex.SetNDC()
        tex.SetTextSize(0.2)
        
    if r1 > rmax:
        tex.SetTextColor(common.mcfit)
        
    return tex
    
