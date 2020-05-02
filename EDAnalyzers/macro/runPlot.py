#!/usr/bin/env python

import os

inputjson=['data.json','mc.json']

param='PVnTracks,PVsumTrackPtSq'

outdir = 'jobs/'
if os.path.isdir(outdir):
    os.system("rm -rf "+outdir)
os.system("mkdir "+outdir)

for j in inputjson:
    outfile = 'jobs/'+j.replace('.json','.root')
    os.system('python plot.py --input='+j+' --output='+outfile+' --param='+param)
