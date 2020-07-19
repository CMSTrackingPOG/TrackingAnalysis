#!/usr/bin/env python

import os

inputjson = ['data.json', 'mc.json']
#inputjson = ['data_zb.json', 'mc_zb.json']

param = 'sumTrackPtSq'

outdir = 'jobs_test/'
if os.path.isdir(outdir):
    os.system("rm -rf "+outdir)
os.system("mkdir "+outdir)

for j in inputjson:
    outfile = outdir+j.replace('.json','.root')
#    os.system('python plot.py --input='+j+' --output='+outfile+' --param='+param+' --time')
    os.system('python plot.py --input='+j+' --output='+outfile+' --param='+param)
