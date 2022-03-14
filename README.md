# TrackingAnalysis

Measurement of the resolution of the primary vertex and track impact
parameters.

Install:
```
cmsrel CMSSW_10_6_28
cd CMSSW_10_6_28/src
cmsenv

git cms-init

git clone git@github.com:kskovpen/TrackingAnalysis.git

scram b -j 8
```

Submit jobs with CRAB:
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
cd crab/
./submit.zsh
```
