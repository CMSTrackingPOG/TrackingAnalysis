#!/usr/bin/env python

import os, sys

#process = 'zb'
process = 'qcd'

#selection = ''
selection = '_pt3p0to10p0'

#mode = 'pv'
#mode = 'ippvdr'
mode = 'ipbseta'

inputpv = 'results/pv_'+process+'.json'
input = 'results/ip_'+process+selection+'.json'

sys=' '
#sys=' --sys'

os.system('python draw.py --mode='+mode+' --input='+input+' --inputpv='+inputpv+' --process='+process+sys+' --selection='+selection)
