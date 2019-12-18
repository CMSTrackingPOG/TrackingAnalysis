import ROOT
import numpy
import math
from array import array

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
       
