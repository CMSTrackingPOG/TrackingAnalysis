#!/usr/bin/env python

import ROOT
import re
import math

#eta = [-3.0, -2.5, -1.3, -0.5, 0., 0.5, 1.3, 2.5, 2.8]
#pt = [0.0, 1.0, 3.0, 10.0]

eta = [-2.5, -1.4, 0., 1.4, 2.4]
pt = [1.0, 3.0, 10.0]

#figIP = ['ipbsetasel_reso_d0_qcd_deconv', \
#'ippvetasel_reso_dz_qcd_deconv', \
#'ipeta_reso_d0_qcd_deconv', \
#'ippveta_reso_dz_qcd_deconv', \
#'ippvpt_reso_dz_qcd_deconv', \
#'ippt_reso_d0_qcd_deconv']

figIP = ['ipbsetasel_reso_d0_qcd_deconv', \
'ipbsetasel_reso_d0_zb_deconv', \
'ipeta_reso_d0_qcd_deconv', \
'ipeta_reso_d0_zb_deconv']
#'ippt_reso_d0_qcd_deconv',\
#'ippt_reso_d0_zb_deconv']

for fig in figIP:

    f = ROOT.TFile.Open('pub_dps/'+fig+'.root')

    c1 = f.Get('c1')
    
    print '----------', fig, '----------'

    next = ROOT.TIter(c1.GetListOfPrimitives())
    
    for i in next:
        
        titl = i.GetTitle()
        name = i.GetName()
        
        if 'eta' in fig:
            if 'h_reso_data' in titl:
                nbins = i.GetXaxis().GetNbins()
                etares = []
                for b in range(nbins+1):
                    c = i.GetXaxis().GetBinCenter(b+1)
                    v = i.GetBinContent(b+1)
                    for e in range(len(etares), len(eta)):
                        if c > eta[e]:
                            etares.append(v)
                            break
                print fig, titl
                if 'd0' in titl:
                    print '|eta|<1.4:', str(etares[2])+'-'+str((etares[1]+etares[3])/2.)
                    print '|eta|<2.5:', str(etares[2])+'-'+str((etares[0]+etares[4])/2.)
#                    print '|eta|<2.5:', str((etares[1]+etares[7])/2.)+'-'+str((etares[0]+etares[8])/2.)
                else:
                    print '|eta|<1.4:', str((etares[3]+etares[5])/2.)+'-'+str((etares[2]+etares[6])/2.)
                    print '|eta|<2.5:', str((etares[1]+etares[7])/2.)+'-'+str((etares[0]+etares[8])/2.)

        if 'pt' in fig:
            if 'h_reso_data' in titl:
                nbins = i.GetXaxis().GetNbins()
                ptres = []
                for b in range(nbins+1):
                    c = i.GetXaxis().GetBinCenter(b+1)
                    v = i.GetBinContent(b+1)
                    for p in range(len(ptres), len(pt)):
                        if c > pt[p]:
                            ptres.append(v)
                            break
                print fig, titl
                print '0.1<pt<1.0:', str(ptres[1])+'-'+str(ptres[0])
                print '1.0<pt<3.0:', str(ptres[2])+'-'+str(ptres[1])
                print '3.0<pt<10.0:', str(ptres[3])+'-'+str(ptres[2])
                    
