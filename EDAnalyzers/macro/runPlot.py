#!/usr/bin/env python

import os

inputjson=['data.json','mc.json']

#param='PVnTracks'
param='PVsumTrackPt'

for j in inputjson:
    outfile = 'jobs/'+j.replace('.json','.root')
    os.system('python plot.py --input='+j+' --output='+outfile+' --param='+param)
