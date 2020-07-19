#!/usr/bin/env python

import os, sys

process = 'zb'
#process = 'qcd'

#selection = ''
#selection = '_pt1p0to10p0'
selection = '_pt1p0to10p0_abseta0p0to1p4'
#selection = '_abseta0p0to1p4'

#mode = 'pv'
#mode = 'ippveta'
mode = 'ipbseta'

inputpv = 'results/pv_'+process+'.json'
input = 'results/ip_'+process+selection+'.json'

#sys=' '
sys=' --sys'

pub=' '
#pub=' --pub'

os.system('python draw.py --mode='+mode+' --input='+input+' --inputpv='+inputpv+' --process='+process+sys+' --selection='+selection+' '+pub)
