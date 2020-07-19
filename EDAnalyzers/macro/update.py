#!/usr/bin/env python

import ROOT
import re
import math

figBS = ['pvXY_data', 'pvXZ_data', 'pvYZ_data']

#figIP = ['ipbsetasel_reso_d0_qcd_deconv', 'ipbsetasel_reso_d0_zb_deconv', \
#'ipeta_reso_d0_qcd_deconv', 'ipeta_reso_d0_zb_deconv', \
#'ipnpv_reso_d0_qcd_deconv', 'ipnpv_reso_d0_zb_deconv', \
#'ipphi_reso_d0_qcd_deconv', 'ipphi_reso_d0_zb_deconv', \
#'ippt_reso_d0_qcd_deconv', 'ippt_reso_d0_zb_deconv', \
#'ippveta_reso_dz_qcd_deconv', 'ippveta_reso_dz_zb_deconv', \
#'ippvetasel_reso_d0_qcd_deconv', 'ippvetasel_reso_d0_zb_deconv', \
#'ippvetasel_reso_dz_qcd_deconv', 'ippvetasel_reso_dz_zb_deconv', \
#'ippvnpv_reso_dz_qcd_deconv', 'ippvnpv_reso_dz_zb_deconv', \
#'ippvphi_reso_dz_qcd_deconv', 'ippvphi_reso_dz_zb_deconv', \
#'ippvpt_reso_dz_qcd_deconv', 'ippvpt_reso_dz_zb_deconv']

#figIP = ['ipbsetasel_reso_d0_qcd_deconv', 'ipbsetasel_reso_d0_zb_deconv', \
#'ippt_reso_d0_qcd_deconv', 'ippt_reso_d0_zb_deconv']
figIP = ['ippt_reso_d0_zb_deconv', 'ippt_reso_d0_zb_deconv_log', 'ipbsetasel_reso_d0_zb_deconv']

figPV = []
#figPV = ['pv_pull_x_zb', 'pv_pull_y_zb', 'pv_pull_z_zb', \
#'pv_pull_x_qcd', 'pv_pull_y_qcd', 'pv_pull_z_qcd', \
#'pv_reso_x_zb', 'pv_reso_y_zb', 'pv_reso_z_zb', \
#'pv_reso_x_qcd', 'pv_reso_y_qcd', 'pv_reso_z_qcd' ]

for fig in figBS+figIP+figPV:

    f = ROOT.TFile.Open('pub_dps2/'+fig+'.root')

    if fig != 'pvXZ_data': c1 = f.Get('c1')
    else: c1 = f.Get('c1_n2')

    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetErrorX(0.0001)

    prelim = None

    next = ROOT.TIter(c1.GetListOfPrimitives())
    for i in next:
        titl = i.GetTitle()
        name = i.GetName()
#        if 'ippt' in fig: 
#            if 'Legend' in titl:
#                i.Clear()
        if titl == 'Work in progress': prelim = i
        if name == 'hframe' and 'pull' in fig:
            i.GetYaxis().SetTitle(i.GetYaxis().GetTitle().replace('Primary vertex', 'PV'))
        if 'x = ' in titl or 'y = ' in titl or 'z = ' in titl:
            i.Clear()
#            fl = re.findall("\d+\.\d+", titl)[0]
#            rou = int(round(float(fl)))
#            i.SetText(i.GetX(), i.GetY()+0.02, titl.replace(fl, str(rou)))
        if 'RMS' in titl:
            i.Clear()
#            i.SetText(i.GetX(), i.GetY()+0.02, titl)
        if name in 'TFrame' and ('pvXY' in fig or 'pvXZ' in fig or 'pvYZ' in fig):
            ROOT.gPad.SetLogz()
        if 'h2' in name and ('pvXY' in fig or 'pvXZ' in fig or 'pvYZ' in fig):
            i.SetMinimum(3.)
            i.SetMaximum(300000.)
            i.GetZaxis().SetTitle('Number of events')
        if 'ZeroBias' in titl:
            if 'ip' in fig:
                if 'phi' not in fig and 'npv' not in fig:
                    i.SetText(i.GetX()-0.08, i.GetY()+0.04, 'Unbiased collision events (Data)')
                else:
                    i.SetText(i.GetX()-0.12, i.GetY(), 'Unbiased collision events (Data)')                    
            else:
                i.SetText(i.GetX()+0.2, i.GetY(), 'Unbiased collision events (Data)')
                if 'XY' in fig or 'XZ' in fig or 'YZ' in fig: 
                    i.SetText(i.GetX()+0.22, i.GetY()-0.58, 'Unbiased collision events (Data)')
                    i.SetTextSize(i.GetTextSize()*0.9)
        if 'QCD' in titl:
            if 'ip' in fig:
                if 'phi' not in fig and 'npv' not in fig:
                    i.SetText(i.GetX()-0.08, i.GetY()+0.04, 'High-q^{2} multi jet events (Data)')
                else:
                    i.SetText(i.GetX()-0.12, i.GetY(), 'High-q^{2} multi jet events (Data)')
            else:
                i.SetText(i.GetX()+0.2, i.GetY(), 'High-q^{2} multi jet events (Data)')
                if 'XY' in fig or 'XZ' in fig or 'YZ' in fig:
                    i.SetText(i.GetX()+0.22, i.GetY()-0.58, 'High-q^{2} multi jet events (Data)')
                    i.SetTextSize(i.GetTextSize()*0.9)

    prelim.SetText(prelim.GetX()+0.08, prelim.GetY(), 'Preliminary')
    
    c1.Update()
    c1.Modified()
    
    c1.Print('pub_dps2/'+fig+'_updated.pdf')    
