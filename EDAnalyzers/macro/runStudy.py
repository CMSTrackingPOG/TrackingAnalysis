#!/usr/bin/env python

import os, sys

jobDir = 'jobs/'

inputMC = ['MinBiasTune4C8TeVpythia8_Summer12DR53XPUS10START53V7Av1AODSIM']
inputData = ['MinimumBias_Run2012A22Jan2013v1AOD','MinimumBias_Run2012B22Jan2013v1AOD',\
'MinimumBias_Run2012C22Jan2013v1AOD','MinimumBias_Run2012D22Jan2013v1AOD']

for f in range(len(inputMC)): inputMC[f] = jobDir+inputMC[f]+'.root'
for f in range(len(inputData)): inputData[f] = jobDir+inputData[f]+'.root'

fMC = jobDir+'mc.root'
fData = jobDir+'data.root'

filesMC = " ".join(inputMC)
#os.system('hadd -f '+fMC+' '+filesMC)

filesData = " ".join(inputData)
#os.system('hadd -f '+fData+' '+filesData)

print "Run PV study"
os.system('python pvStudy.py --data='+fData+' --mc='+fMC)

#print "Run IP study"
#os.system('python ipStudy.py --data='+fData+' --mc='+fMC)
