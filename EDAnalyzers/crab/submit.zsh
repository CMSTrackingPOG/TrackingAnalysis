#!/bin/env zsh

# source /cvmfs/cms.cern.ch/crab3/crab.sh

slist="list.txt"
pver="1" # production tentative
pset="crabConfigTemplate.py"
psetData="crabConfigTemplateData.py"
ver="Track-v20191229"
prodv="/store/user/kskovpen/Track/Ntuple/${ver}/"

rm -f crabConfig.py*

samp=()
is=1
cat ${slist} | while read line
do
  if [[ ${line[1]} == '#' ]]; then
    continue
  fi
  samp[${is}]=${line}
  is=$[$is+1]
done

for i in ${samp}
do
  spl=($(echo $i | tr "/" "\n"))
  pubdn=$(echo "${spl[2]}_${spl[3]}" | sed 's%-%_%g')
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
  | sed "s%OUTLFN%${prodv}%g" \
  | sed "s%REQUESTNAME%${nam}_${pubdn}_${pver}%g" \
  | sed "s%PUBLISHDATANAME%${pubdn}%g" \
  > crabConfig.py

  echo "${nam} ${pubdn}"
  crab submit
  
done

rm -f crabConfig.py*
