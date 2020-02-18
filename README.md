# TrackingAnalysis

Measurement of the resolution of the primary vertex and track impact
parameter (IP) reconstruction.

Inspired by the IP resolution studies done by B. Mangano et al. (see http://cmscvs.web.cern.ch/cmscvs/cgi/viewvc.cgi/cvsroot/UserCode/Mangano/IpResoStudies/)

Install (Run 2) on SL6:
```
cmsrel CMSSW_10_2_18
cd CMSSW_10_2_18/src
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
