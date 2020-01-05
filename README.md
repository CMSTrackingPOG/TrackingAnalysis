# TrackingAnalysis

Measurement of the resolution of the primary vertex and track impact
parameter (IP) reconstruction.

Inspired by the IP resolution studies done by B. Mangano et al. (see http://cmscvs.web.cern.ch/cmscvs/cgi/viewvc.cgi/cvsroot/UserCode/Mangano/IpResoStudies/)

Install (Run I) on SL6:
```
cmsrel CMSSW_5_3_39
cd CMSSW_5_3_39/src
cmsenv

git cms-init

git cms-addpkg RecoVertex/PrimaryVertexProducer
git checkout 4e3a97503283ff2d2218b2da8e23a8e345e11e86

git clone git@github.com:kskovpen/TrackingAnalysis.git

scram b -j 8
```

Submit jobs with CRAB:
```
source /cvmfs/cms.cern.ch/crab3/crab_slc6.sh
source /cvmfs/cms.cern.ch/slc6_amd64_gcc493/external/python/2.7.11/etc/profile.d/init.sh
cd crab/
./submit.zsh
```
