import numpy as np

datapath = '../test/'
year = '2017'

mccol = 45
mcfit = 2
datamrk = 20
mcmrk = 21

ncores = 8

PVmeas = ['x','y','z']
IPmeas = ['d0','dz']

hPVData = {'x':[150,0.6,1.1],'y':[150,-0.65,-0.05],'z':[150,-25.,25.]}
hPVMC = {'x':[150,-0.5,0.0],'y':[150,0.5,0.9],'z':[150,-25.,25.]}

hBSData = {'bwx':[100,0.,30.],'bwy':[100,0.,30.],'bwz':[100,25.,45.]}
hBSMC = {'bwx':[100,0.,30.],'bwy':[100,0.,30.],'bwz':[100,25.,45.]}

hPVnTrks = {'zb':[50,0.,200.],'qcd':[50,0.,250.]}
hPVsumTrackPt = {'zb':[70,0.,200.],'qcd':[70,0.,1700.]}
hPVsumTrackPt2 = {'zb':[70,0.,40.],'qcd':[70,0.,800.]}

hJetPtMax = {'zb':[70,0.,500.],'qcd':[70,0.,2000.]}
hJetHT = {'zb':[70,0.,1000.],'qcd':[70,0.,4000.]}

sel = {\
'pt':{'_pt0p0to1p0':[0.0,1.0],'_pt1p0to3p0':[1.0,3.0],'_pt3p0to10p0':[3.0,10.0]},\
'eta':{'_abseta0p0to1p0':[0.0,1.0],'_abseta1p0to1p5':[1.0,1.5],'_abseta1p5to2p0':[1.5,2.0],'_abseta2p0to3p0':[2.0,3.0]}\
}

