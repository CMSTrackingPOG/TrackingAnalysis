#!/bin/bash

pileup=Collisions17_13TeV_PileUp_pileup_latest.txt
json=../../../crab/JSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt

pileupCalc.py --verbose -i $json --inputLumiJSON $pileup --calcMode true --minBiasXsec 66000 --maxPileupBin 100 --numPileupBins 100 puDown.root
pileupCalc.py --verbose -i $json --inputLumiJSON $pileup --calcMode true --minBiasXsec 69200 --maxPileupBin 100 --numPileupBins 100 puCentral.root
pileupCalc.py --verbose -i $json --inputLumiJSON $pileup --calcMode true --minBiasXsec 72400 --maxPileupBin 100 --numPileupBins 100 puUp.root
