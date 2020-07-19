#!/bin/env zsh

# d0 vs pt, zb, with selection
./comp.py \
--param=pt \
--measurement=ip \
--type=bs,bs \
--pkl=reso_d0_abseta0p0to1p4_zb_data_deconv,reso_d0_zb_data_deconv \
--selection

mv pub/comp.pdf pub_dps2/ippt_reso_d0_zb_deconv.pdf
mv pub/comp.root pub_dps2/ippt_reso_d0_zb_deconv.root

# d0 vs pt, zb, with selection, log scale
./comp.py \
--param=pt \
--measurement=ip \
--type=bs,bs \
--pkl=reso_d0_abseta0p0to1p4_zb_data_deconv,reso_d0_zb_data_deconv \
--selection \
--log

mv pub/comp.pdf pub_dps2/ippt_reso_d0_zb_deconv_log.pdf
mv pub/comp.root pub_dps2/ippt_reso_d0_zb_deconv_log.root

# d0 vs eta, zb, bs, with selection
./comp.py \
--param=eta \
--measurement=ip \
--type=bs,bs \
--pkl=reso_d0_pt1p0to10p0_abseta0p0to1p4_zb_data_deconv,reso_d0_pt1p0to10p0_zb_data_deconv \
--selection

mv pub/comp.pdf pub_dps2/ipbsetasel_reso_d0_zb_deconv.pdf
mv pub/comp.root pub_dps2/ipbsetasel_reso_d0_zb_deconv.root




# d0 vs pt, zb
#./comp.py \
#--param=pt \
#--measurement=ip \
#--type=bs \
#--pkl=reso_d0_zb_data_deconv

#mv pub/comp.pdf pub_dps/ippt_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub_dps/ippt_reso_d0_zb_deconv.root

# d0 vs pt, qcd
#./comp.py \
#--param=pt \
#--measurement=ip \
#--type=bs \
#--pkl=reso_d0_qcd_data_deconv

#mv pub/comp.pdf pub_dps/ippt_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub_dps/ippt_reso_d0_qcd_deconv.root

# d0 vs eta, zb, bs, with selection
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=bs,bs,bs,bs \
#--pkl=reso_d0_pt1p0to3p0_abseta0p0to1p4_zb_data_deconv,reso_d0_pt1p0to3p0_zb_data_deconv \
#--selection

#mv pub/comp.pdf pub_dps/ipbsetasel_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub_dps/ipbsetasel_reso_d0_zb_deconv.root

# d0 vs eta, qcd, bs, with selection
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=bs,bs,bs,bs \
#--pkl=reso_d0_pt1p0to10p0_abseta0p0to1p4_qcd_data_deconv,reso_d0_pt1p0to10p0_qcd_data_deconv,reso_d0_pt3p0to10p0_abseta0p0to1p4_qcd_data_deconv,reso_d0_pt3p0to10p0_qcd_data_deconv \
#--selection

#mv pub/comp.pdf pub_dps/ipbsetasel_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub_dps/ipbsetasel_reso_d0_qcd_deconv.root

# d0 vs eta, zb
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=bs \
#--pkl=reso_d0_pt1p0to3p0_zb_data_deconv

#mv pub/comp.pdf pub_dps/ipeta_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub_dps/ipeta_reso_d0_zb_deconv.root

# d0 vs eta, qcd
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=bs \
#--pkl=reso_d0_pt1p0to10p0_qcd_data_deconv

#mv pub/comp.pdf pub_dps/ipeta_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub_dps/ipeta_reso_d0_qcd_deconv.root








# d0 vs pt, zb
#./comp.py \
#--param=pt \
#--measurement=ip \
#--type=pv,bs,pv,bs \
#--pkl=reso_d0_zb_data_deconv,reso_d0_zb_data_deconv,reso_d0_zb_mc_deconv,reso_d0_zb_mc_deconv

#mv pub/comp.pdf pub/ippt_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub/ippt_reso_d0_zb_deconv.root

# d0 vs pt, qcd
#./comp.py \
#--param=pt \
#--measurement=ip \
#--type=pv,bs,pv,bs \
#--pkl=reso_d0_qcd_data_deconv,reso_d0_qcd_data_deconv,reso_d0_qcd_mc_deconv,reso_d0_qcd_mc_deconv

#mv pub/comp.pdf pub/ippt_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub/ippt_reso_d0_qcd_deconv.root

# d0 vs eta, zb
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=pv,bs,pv,bs \
#--pkl=reso_d0_zb_data_deconv,reso_d0_zb_data_deconv,reso_d0_zb_mc_deconv,reso_d0_zb_mc_deconv

#mv pub/comp.pdf pub/ipeta_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub/ipeta_reso_d0_zb_deconv.root

# d0 vs eta, qcd
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=pv,bs,pv,bs \
#--pkl=reso_d0_qcd_data_deconv,reso_d0_qcd_data_deconv,reso_d0_qcd_mc_deconv,reso_d0_qcd_mc_deconv

#mv pub/comp.pdf pub/ipeta_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub/ipeta_reso_d0_qcd_deconv.root

# d0 vs phi, zb
#./comp.py \
#--param=phi \
#--measurement=ip \
#--type=pv,bs,pv,bs \
#--pkl=reso_d0_zb_data_deconv,reso_d0_zb_data_deconv,reso_d0_zb_mc_deconv,reso_d0_zb_mc_deconv

#mv pub/comp.pdf pub/ipphi_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub/ipphi_reso_d0_zb_deconv.root

# d0 vs phi, qcd
#./comp.py \
#--param=phi \
#--measurement=ip \
#--type=pv,bs,pv,bs \
#--pkl=reso_d0_qcd_data_deconv,reso_d0_qcd_data_deconv,reso_d0_qcd_mc_deconv,reso_d0_qcd_mc_deconv

#mv pub/comp.pdf pub/ipphi_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub/ipphi_reso_d0_qcd_deconv.root

# d0 vs npv, zb
#./comp.py \
#--param=npv \
#--measurement=ip \
#--type=pv,bs,pv,bs \
#--pkl=reso_d0_zb_data_deconv,reso_d0_zb_data_deconv,reso_d0_zb_mc_deconv,reso_d0_zb_mc_deconv

#mv pub/comp.pdf pub/ipnpv_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub/ipnpv_reso_d0_zb_deconv.root

# d0 vs npv, qcd
#./comp.py \
#--param=npv \
#--measurement=ip \
#--type=pv,bs,pv,bs \
#--pkl=reso_d0_qcd_data_deconv,reso_d0_qcd_data_deconv,reso_d0_qcd_mc_deconv,reso_d0_qcd_mc_deconv

#mv pub/comp.pdf pub/ipnpv_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub/ipnpv_reso_d0_qcd_deconv.root

# d0 vs eta, zb, bs, with selection
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=bs,bs,bs,bs \
#--pkl=reso_d0_pt0p0to1p0_zb_data_deconv,reso_d0_pt0p0to1p0_zb_mc_deconv,reso_d0_pt1p0to3p0_zb_data_deconv,reso_d0_pt1p0to3p0_zb_mc_deconv \
#--selection

#mv pub/comp.pdf pub/ipbsetasel_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub/ipbsetasel_reso_d0_zb_deconv.root

# d0 vs eta, zb, pv, with selection
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=pv,pv,pv,pv \
#--pkl=reso_d0_pt0p0to1p0_zb_data_deconv,reso_d0_pt0p0to1p0_zb_mc_deconv,reso_d0_pt1p0to3p0_zb_data_deconv,reso_d0_pt1p0to3p0_zb_mc_deconv \
#--selection

#mv pub/comp.pdf pub/ippvetasel_reso_d0_zb_deconv.pdf
#mv pub/comp.root pub/ippvetasel_reso_d0_zb_deconv.root

# d0 vs eta, qcd, bs, with selection
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=bs,bs,bs,bs,bs,bs \
#--pkl=reso_d0_pt0p0to1p0_qcd_data_deconv,reso_d0_pt0p0to1p0_qcd_mc_deconv,reso_d0_pt1p0to3p0_qcd_data_deconv,reso_d0_pt1p0to3p0_qcd_mc_deconv,reso_d0_pt3p0to10p0_qcd_data_deconv,reso_d0_pt3p0to10p0_qcd_mc_deconv \
#--selection

#mv pub/comp.pdf pub/ipbsetasel_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub/ipbsetasel_reso_d0_qcd_deconv.root

# d0 vs eta, qcd, pv, with selection
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=pv,pv,pv,pv,pv,pv \
#--pkl=reso_d0_pt0p0to1p0_qcd_data_deconv,reso_d0_pt0p0to1p0_qcd_mc_deconv,reso_d0_pt1p0to3p0_qcd_data_deconv,reso_d0_pt1p0to3p0_qcd_mc_deconv,reso_d0_pt3p0to10p0_qcd_data_deconv,reso_d0_pt3p0to10p0_qcd_mc_deconv \
#--selection

#mv pub/comp.pdf pub/ippvetasel_reso_d0_qcd_deconv.pdf
#mv pub/comp.root pub/ippvetasel_reso_d0_qcd_deconv.root

# dz vs eta, zb
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=pv,pv \
#--pkl=reso_dz_zb_data_deconv,reso_dz_zb_mc_deconv

#mv pub/comp.pdf pub/ippveta_reso_dz_zb_deconv.pdf
#mv pub/comp.root pub/ippveta_reso_dz_zb_deconv.root

# dz vs eta, qcd
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=pv,pv \
#--pkl=reso_dz_qcd_data_deconv,reso_dz_qcd_mc_deconv

#mv pub/comp.pdf pub/ippveta_reso_dz_qcd_deconv.pdf
#mv pub/comp.root pub/ippveta_reso_dz_qcd_deconv.root

# dz vs pt, zb
#./comp.py \
#--param=pt \
#--measurement=ip \
#--type=pv,pv \
#--pkl=reso_dz_zb_data_deconv,reso_dz_zb_mc_deconv

#mv pub/comp.pdf pub/ippvpt_reso_dz_zb_deconv.pdf
#mv pub/comp.root pub/ippvpt_reso_dz_zb_deconv.root

# dz vs pt, qcd
#./comp.py \
#--param=pt \
#--measurement=ip \
#--type=pv,pv \
#--pkl=reso_dz_qcd_data_deconv,reso_dz_qcd_mc_deconv

#mv pub/comp.pdf pub/ippvpt_reso_dz_qcd_deconv.pdf
#mv pub/comp.root pub/ippvpt_reso_dz_qcd_deconv.root

# dz vs phi, zb
#./comp.py \
#--param=phi \
#--measurement=ip \
#--type=pv,pv \
#--pkl=reso_dz_zb_data_deconv,reso_dz_zb_mc_deconv

#mv pub/comp.pdf pub/ippvphi_reso_dz_zb_deconv.pdf
#mv pub/comp.root pub/ippvphi_reso_dz_zb_deconv.root

# dz vs phi, qcd
#./comp.py \
#--param=phi \
#--measurement=ip \
#--type=pv,pv \
#--pkl=reso_dz_qcd_data_deconv,reso_dz_qcd_mc_deconv

#mv pub/comp.pdf pub/ippvphi_reso_dz_qcd_deconv.pdf
#mv pub/comp.root pub/ippvphi_reso_dz_qcd_deconv.root

# dz vs npv, zb
#./comp.py \
#--param=npv \
#--measurement=ip \
#--type=pv,pv \
#--pkl=reso_dz_zb_data_deconv,reso_dz_zb_mc_deconv

#mv pub/comp.pdf pub/ippvnpv_reso_dz_zb_deconv.pdf
#mv pub/comp.root pub/ippvnpv_reso_dz_zb_deconv.root

# dz vs npv, qcd
#./comp.py \
#--param=npv \
#--measurement=ip \
#--type=pv,pv \
#--pkl=reso_dz_qcd_data_deconv,reso_dz_qcd_mc_deconv

#mv pub/comp.pdf pub/ippvnpv_reso_dz_qcd_deconv.pdf
#mv pub/comp.root pub/ippvnpv_reso_dz_qcd_deconv.root

# dz vs eta, zb, with selection
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=pv,pv,pv,pv \
#--pkl=reso_dz_pt0p0to1p0_zb_data_deconv,reso_dz_pt0p0to1p0_zb_mc_deconv,reso_dz_pt1p0to3p0_zb_data_deconv,reso_dz_pt1p0to3p0_zb_mc_deconv \
#--selection

#mv pub/comp.pdf pub/ippvetasel_reso_dz_zb_deconv.pdf
#mv pub/comp.root pub/ippvetasel_reso_dz_zb_deconv.root

# dz vs eta, qcd, with selection
#./comp.py \
#--param=eta \
#--measurement=ip \
#--type=pv,pv,pv,pv,pv,pv \
#--pkl=reso_dz_pt0p0to1p0_qcd_data_deconv,reso_dz_pt0p0to1p0_qcd_mc_deconv,reso_dz_pt1p0to3p0_qcd_data_deconv,reso_dz_pt1p0to3p0_qcd_mc_deconv,reso_dz_pt3p0to10p0_qcd_data_deconv,reso_dz_pt3p0to10p0_qcd_mc_deconv \
#--selection

#mv pub/comp.pdf pub/ippvetasel_reso_dz_qcd_deconv.pdf
#mv pub/comp.root pub/ippvetasel_reso_dz_qcd_deconv.root
