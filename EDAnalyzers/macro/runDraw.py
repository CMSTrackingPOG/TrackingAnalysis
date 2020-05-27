#!/usr/bin/env python

import os, sys

#process = 'zb'
process = 'qcd'

#mode = 'pv'
mode = 'ipbspt'
#mode = 'ipbsdr'
inputpv = 'results/pv_'+process+'.json'
input = 'results/ip_'+process+'.json'

sys=' '
#sys=' --sys'

os.system('python draw.py --mode='+mode+' --input='+input+' --inputpv='+inputpv+' --process='+process+sys)
