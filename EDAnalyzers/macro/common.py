import numpy as np

datapath = '../test/'

PVmeas = ['x','y','z']
IPmeas = ['d0','dz']

PVData = {'x':[0.8,1.1],'y':[-0.9,-0.4],'z':[-25.,25.]}
PVMC = {'x':[-0.05,0.25],'y':[0.3,0.55],'z':[-25.,25.]}

BSData = {'bwx':[100,0.,30.],'bwy':[100,0.,30.],'bwz':[100,0.,10.]}
BSMC = {'bwx':[100,0.,30.],'bwy':[100,0.,30.],'bwz':[100,0.,10.]}

PVnTracksBins = np.array([0.,20.,30.,40.,50.,60.,70.,80.,100.,150.], dtype='float64')

PVnTracks = {\
''              :{'bins':[0,20,1000],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-200.,200.],'resoy':[50,-200.,200.],'resoz':[50,-400.,400.]},\
'_nTrks0to20'   :{'bins':[1,0,20],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-200.,200.],'resoy':[50,-200.,200.],'resoz':[50,-400.,400.]},\
'_nTrks20to30'  :{'bins':[2,20,30],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-150.,150.],'resoy':[50,-150.,150.],'resoz':[50,-300.,300.]},\
'_nTrks30to40'  :{'bins':[3,30,40],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-150.,150.],'resoy':[50,-150.,150.],'resoz':[50,-200.,200.]},\
'_nTrks40to50'  :{'bins':[4,40,50],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-150.,150.],'resoy':[50,-150.,150.],'resoz':[50,-150.,150.]},\
'_nTrks50to60'  :{'bins':[5,50,60],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-120.,120.],'resoy':[50,-120.,120.],'resoz':[50,-150.,150.]},\
'_nTrks60to70'  :{'bins':[6,60,70],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-100.,100.],'resoy':[50,-100.,100.],'resoz':[50,-120.,120.]},\
'_nTrks70to80'  :{'bins':[7,70,80],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-100.,100.],'resoy':[50,-100.,100.],'resoz':[50,-120.,120.]},\
'_nTrks80to100' :{'bins':[8,80,100],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-100.,100.],'resoy':[50,-100.,100.],'resoz':[50,-120.,120.]},\
'_nTrks100to150':{'bins':[9,100,150],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-80.,80.],'resoy':[50,-80.,80.],'resoz':[50,-100.,100.]}\
}

PVsumTrackPtBins = np.array([0.,20.,30.,40.,50.,60.,70.,80.,100.,150.], dtype='float64')

PVsumTrackPt = {\
''                   :{'bins':[0,20,1000],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-200.,200.],'resoy':[50,-200.,200.],'resoz':[50,-300.,300.]},\
'_sumTrackPt0to20'   :{'bins':[1,0,20],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-200.,200.],'resoy':[50,-200.,200.],'resoz':[50,-300.,300.]},\
'_sumTrackPt20to30'  :{'bins':[2,20,30],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-150.,150.],'resoy':[50,-150.,150.],'resoz':[50,-200.,200.]},\
'_sumTrackPt30to40'  :{'bins':[3,30,40],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-150.,150.],'resoy':[50,-150.,150.],'resoz':[50,-200.,200.]},\
'_sumTrackPt40to50'  :{'bins':[4,40,50],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-150.,150.],'resoy':[50,-150.,150.],'resoz':[50,-150.,150.]},\
'_sumTrackPt50to60'  :{'bins':[5,50,60],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-120.,120.],'resoy':[50,-120.,120.],'resoz':[50,-150.,150.]},\
'_sumTrackPt60to70'  :{'bins':[6,60,70],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-100.,100.],'resoy':[50,-100.,100.],'resoz':[50,-120.,120.]},\
'_sumTrackPt70to80'  :{'bins':[7,70,80],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-100.,100.],'resoy':[50,-100.,100.],'resoz':[50,-120.,120.]},\
'_sumTrackPt80to100' :{'bins':[8,80,100],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-100.,100.],'resoy':[50,-100.,100.],'resoz':[50,-120.,120.]},\
'_sumTrackPt100to150':{'bins':[9,100,150],'pullx':[50,-4.0,4.0],'pully':[50,-4.0,4.0],'pullz':[50,-4.0,4.0],'resox':[50,-80.,80.],'resoy':[50,-80.,80.],'resoz':[50,-100.,100.]}\
}

IPptBins = np.array([0.3,0.5,0.75,1.,1.25,1.5,1.75,2.,2.25,2.5,2.75,3.,3.5,4.,5.], dtype='float64')

IPpt = {\
''            :{'bins':[0,0.3,1000],'d0':[100,-1800,1800],'dz':[100,-3000,3000]},\
'_pt0p3to0p5' :{'bins':[1,0.3,0.5],'d0':[100,-1800,1800],'dz':[100,-3000,3000]},\
'_pt0p5to0p75':{'bins':[2,0.5,0.75],'d0':[100,-1500,1500],'dz':[100,-2000,2000]},\
'_pt0p75to1'  :{'bins':[3,0.75,1],'d0':[100,-1000,1000],'dz':[100,-1500,1500]},\
'_pt1to1p25'  :{'bins':[4,1,1.25],'d0':[100,-1000,1000],'dz':[100,-1500,1500]},\
'_pt1p25to1p5':{'bins':[5,1.25,1.5],'d0':[100,-800,800],'dz':[100,-1000,1000]},\
'_pt1p5to1p75':{'bins':[6,1.5,1.75],'d0':[100,-700,700],'dz':[100,-1000,1000]},\
'_pt1p75to2'  :{'bins':[7,1.75,2],'d0':[100,-600,600],'dz':[100,-1000,1000]},\
'_pt2to2p25'  :{'bins':[8,2,2.25],'d0':[100,-600,600],'dz':[100,-800,800]},\
'_pt2p25to2p5':{'bins':[9,2.25,2.5],'d0':[100,-500,500],'dz':[100,-800,800]},\
'_pt2p5to2p75':{'bins':[10,2.5,2.75],'d0':[100,-500,500],'dz':[100,-800,800]},\
'_pt2p75to3'  :{'bins':[11,2.75,3],'d0':[100,-400,400],'dz':[100,-800,800]},\
'_pt3to3p5'   :{'bins':[12,3,3.5],'d0':[100,-400,400],'dz':[100,-800,800]},\
'_pt3p5to4'   :{'bins':[13,3.5,4],'d0':[100,-300,300],'dz':[100,-800,800]},\
'_pt4to5'     :{'bins':[14,4,5],'d0':[100,-250,250],'dz':[100,-800,800]}\
}

IPetaBins = np.array([0.,0.5,1.,1.5,2.,3.], dtype='float64')

IPeta = {\
''          :{'bins':[0,0,1000],'d0':[100,-1000,1000],'dz':[100,-2000,2000]},\
'_eta0to0p5':{'bins':[1,0,0.5],'d0':[100,-1000,1000],'dz':[100,-1000,1000]},\
'_eta0p5to1':{'bins':[2,0.5,1],'d0':[100,-1000,1000],'dz':[100,-1000,1000]},\
'_eta1to1p5':{'bins':[3,1,1.5],'d0':[100,-1500,1500],'dz':[100,-2000,2000]},\
'_eta1p5to2':{'bins':[4,1.5,2],'d0':[100,-1700,1700],'dz':[100,-4000,4000]},\
'_eta2to3'  :{'bins':[5,2,3],'d0':[100,-2000,2000],'dz':[100,-8000,8000]}\
}
