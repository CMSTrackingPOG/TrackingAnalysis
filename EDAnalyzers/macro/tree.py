import os
import sys
import math
from array import array
import ROOT

class trackTree():

    def __init__(self):
                        
        self.sumTrackPtSq, self.beamWidthX, self.beamWidthY, self.beamSigmaZ, \
        self.beamWidthXError, self.beamWidthYError, self.beamSigmaZError, \
        self.beamEmittanceX, self.beamEmittanceY, self.beamBetaStar \
        = (array( 'f', [ -777 ] ) for _ in range(10))
        
        self.run, self.lumi \
        = (array( 'i', [ -777 ] ) for _ in range(2))
        
        self.pt, self.eta, self.phi, self.dr = (ROOT.vector('float')() for _ in range(4))
        
        self.npv = ROOT.vector('int')()
        
        bsize = 1024*1024 # 1 MB
        nbrs = 17 # Number of branches

        self.t = ROOT.TTree( 'trackTree', 'Track tree' )

        self.t.Branch( 'run', self.run, 'run/I', bsize )
        self.t.Branch( 'lumi', self.lumi, 'lumi/I', bsize )
        
        self.t.Branch( 'sumTrackPtSq', self.sumTrackPtSq, 'sumTrackPtSq/F', bsize )
        self.t.Branch( 'beamWidthX', self.beamWidthX, 'beamWidthX/F', bsize )
        self.t.Branch( 'beamWidthY', self.beamWidthY, 'beamWidthY/F', bsize )
        self.t.Branch( 'beamSigmaZ', self.beamSigmaZ, 'beamSigmaZ/F', bsize )
        self.t.Branch( 'beamWidthXError', self.beamWidthXError, 'beamWidthXError/F', bsize )
        self.t.Branch( 'beamWidthYError', self.beamWidthYError, 'beamWidthYError/F', bsize )
        self.t.Branch( 'beamSigmaZError', self.beamSigmaZError, 'beamSigmaZError/F', bsize )
        self.t.Branch( 'beamEmittanceX', self.beamEmittanceX, 'beamEmittanceX/F', bsize )
        self.t.Branch( 'beamEmittanceY', self.beamEmittanceY, 'beamEmittanceY/F', bsize )
        self.t.Branch( 'beamBetaStar', self.beamBetaStar, 'beamBetaStar/F', bsize )

        self.t.Branch( 'pt', 'std::vector<float>', self.pt, bsize )
        self.t.Branch( 'eta', 'std::vector<float>', self.eta, bsize )
        self.t.Branch( 'phi', 'std::vector<float>', self.phi, bsize )
        self.t.Branch( 'dr', 'std::vector<float>', self.dr, bsize )
        self.t.Branch( 'npv', 'std::vector<int>', self.npv, bsize )
        
        self.t.SetAutoFlush(-bsize*nbrs)

    def clear(self):

        self.pt.clear()
        self.eta.clear()
        self.phi.clear()
        self.dr.clear()
        self.npv.clear()
        
    def fill(self):

        self.t.Fill()
        
