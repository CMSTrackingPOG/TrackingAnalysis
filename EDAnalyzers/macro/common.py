import numpy as np

datapath = '../test/'

PVmeas = ['x','y','z']
IPmeas = ['d0','dz']

PVnTracksBins = np.array([20.,30.,40.,50.,60.,70.,100.], dtype='float64')
PVnTracks = {'':[0,20,1000],'_nTrks20to30':[1,20,30],\
'_nTrks30to40':[2,30,40],'_nTrks40to50':[3,40,50],'_nTrks50to60':[4,50,60],\
'_nTrks60to70':[5,60,70],'_nTrks70to100':[6,70,100]}

PVsumTrackPtBins = np.array([20.,30.,40.,50.,60.,70.,100.], dtype='float64')
PVsumTrackPt = {'':[0,20,1000],'_sumTrackPt20to30':[1,20,30],\
'_sumTrackPt30to40':[2,30,40],'_sumTrackPt40to50':[3,40,50],'_sumTrackPt50to60':[4,50,60],\
'_sumTrackPt60to70':[5,60,70],'_sumTrackPt70to100':[6,70,100]}

IPptBins = np.array([0.3,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.5,4.,5.], dtype='float64')
IPpt = {'':[0,0.3,1000],'_pt0p3to0p5':[1,0.3,0.5],'_pt0p5to0p75':[2,0.5,0.75],'_pt0p75to1':[3,0.75,1],\
'_pt1to1p25':[4,1,1.25],'_pt1p25to1p5':[5,1.25,1.5],'_pt1p5to1p75':[6,1.5,1.75],'_pt1p75to2':[7,1.75,2],\
'_pt2to2p25':[8,2,2.25],'_pt2p25to2p5':[9,2.25,2.5],'_pt2p5to2p75':[10,2.5,2.75],\
'_pt2p75to3':[11,2.75,3],'_pt3to3p5':[12,3,3.5],'_pt3p5to4':[13,3.5,4],'_pt4to5':[14,4,5]}

IPetaBins = np.array([0.,0.5,1.,1.5,2.,3.], dtype='float64')
IPeta = {'':[0,0,1000],'_eta0to0p5':[1,0,0.5],'_eta0p5to1':[2,0.5,1],'_eta1to1p5':[3,1,1.5],\
'_eta1p5to2':[4,1.5,2],'_eta2to3':[5,2,3]}
