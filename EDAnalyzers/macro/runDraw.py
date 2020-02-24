#!/usr/bin/env python

import os, sys

mode = 'pv'
#mode = 'ipbspt'
#mode = 'ippvpt'
input = 'results/pv.json'
#input = 'results/ip.json'

os.system('python draw.py --mode='+mode+' --input='+input)
