# TrackingAnalysis

Measurement of the resolution of the primary vertex and track impact
parameters.

Install:
```
cmsrel CMSSW_11_0_3
cd CMSSW_11_0_3/src
cmsenv

git cms-init

#git cms-addpkg SimDataFormats/Associations
#git cms-addpkg SimTracker/TrackAssociatorProducers
#git cms-addpkg SimTracker/TrackHistory
#git cms-addpkg SimTracker/VertexAssociation

#git cherry-pick 2089bd1669d57a87530d1b208666178b6797fc9a

git clone git@github.com:kskovpen/TrackingAnalysis.git

scram b -j 8
```

Submit jobs with CRAB:
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
cd crab/
./submit.zsh
```
