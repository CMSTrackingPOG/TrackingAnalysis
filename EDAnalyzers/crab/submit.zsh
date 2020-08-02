#!/bin/env zsh

# source /cvmfs/cms.cern.ch/crab3/crab.sh

slist="digi.txt"
pver="1" # production tentative
pset="crabConfigTemplate.py"
psetData="crabConfigTemplateData.py"
doTruth="1"
ver="Track-v20200802"
prodv="/store/user/kskovpen/Track/Ntuple/${ver}/"

rm -f crabConfig.py*

samp=()
sec=()
is=1
cat ${slist} | while read line
do
  if [[ ${line[1]} == '#' ]]; then
    continue
  fi
  ds=(${(ps: :)${line}})
  samp[${is}]=${ds[1]}
  sec[${is}]=${ds[2]}
  is=$[$is+1]
done

is=1
for i in ${samp}
do
  spl=($(echo $i | tr "/" "\n"))
  pubdn=$(echo "${spl[2]}_${spl[3]}" | sed 's%-%_%g' | \
  sed 's%106X.*%%g' | sed 's%09Aug2019.*%%g')
  isdata=$(echo "${spl[2]}" | grep -b -o Run2 | awk '{print $1}')
  run=$(echo "${spl[2]}" | cut -c$[${isdata%:*}+1]-$[${isdata%:*}+8])
  nam=$(echo "${spl[1]}" | sed 's%-%_%g')
  
  if [[ ${isdata} == "" ]]; then
    pubdn=$(echo $pubdn | sed 's%kskovpen.*%MC%g')
  else
    pset=${psetData}
    pubdn=$(echo $pubdn | sed "s%kskovpen.*%${run}%g")
  fi

  cat ${pset} | sed "s%INPUTDATASET%${i}%g" \
  | sed "s%SECONDARYDATASET%${sec[${is}]}%g" \
  | sed "s%OUTLFN%${prodv}%g" \
  | sed "s%doTruth=0%doTruth=${doTruth}%g" \
  | sed "s%REQUESTNAME%${nam}_${pubdn}_${pver}%g" \
  | sed "s%PUBLISHDATANAME%${pubdn}%g" \
  > crabConfig.py
  
  if [[ ${doTruth} == "0" ]]; then
    mv crabConfig.py crabConfigCopy.py
    cat crabConfigCopy.py | sed "s%config.Data.secondaryInputDataset%#config.Data.secondaryInputDataset%g" > crabConfig.py
    rm crabConfigCopy.py
  fi
  
  echo "${nam} ${pubdn}"
#  crab submit -c crabConfig.py --dryrun
  crab submit -c crabConfig.py

  is=$[$is+1]
  
done

rm -f crabConfig.py*
