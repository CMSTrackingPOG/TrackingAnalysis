#!/bin/env python

from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import numpy as np
import os, sys
import json
import uproot
import ROOT

from optparse import OptionParser

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]
        
    usage = "usage: %prog [options]\n Bin optimisation tool"
    
    parser = OptionParser(usage)
    
    parser.add_option("--method",default='v',help="Method for optimisation (variable or constant bin width) [default: %default]")
    parser.add_option("--input",default='ZeroBias.root',help="Input file name [default: %default]")
    parser.add_option("--output",default='data/bins/output.json',help="Output file name [default: %default]")
    parser.add_option("--tree",default='trackTree',help="Input tree name [default: %default]")
    parser.add_option("--nmin",default=10000,help="Minimum number of events per bin in the first-level parameterisation [default: %default]")
    parser.add_option("--nbins",default=100,help="Maximum number of bins [default: %default]")
    parser.add_option("--threads",default=8,help="Number of threads [default: %default]")
    parser.add_option("--crop",default=0.1,help="Crop factor [default: %default]")
    parser.add_option("--meas",default='bs',help="Measurement type (bs, pv, or 1d) [default: %default]")
    parser.add_option("--pv",default='sumTrackPtSq',help="PV parameterisation [default: %default]")
    parser.add_option("--param",default='pt,eta,phi,npv,dr',help="List of track parameterisations [default: %default]")

    (options, args) = parser.parse_args(sys.argv[1:])

    return options

def done():

    sys.stdout.write('\x1b[\33[32mdone\x1b[0m\n'); sys.stdout.flush()

def flush(s, flush=True):

    sys.stdout.write(s)
    if flush: sys.stdout.flush()

def optimise(data, nbins, nmin):
    
    nbinsc = nbins
    
    opt = False; res = []
            
    while not opt:
             
        bl = []
        bh = []
        
        if options.method == 'c': splt = pd.cut(data, nbinsc, duplicates='drop')
        else: splt = pd.qcut(data, nbinsc, duplicates='drop')
        bins = splt.value_counts(sort=False)
        nb = len(zip(bins))
            
        goodstats = True
        for i, (v, nev) in enumerate(bins.iteritems()):
            
            if nev < nmin:
                goodstats = False
                nbinsc -= 1
                break                
            
            b = str(v)
            for s in ['(',')','[',']',' ']: b = b.replace(s,'')
            b = b.split(',')
            bl.append(float(b[0]))
            bh.append(float(b[1]))
                 
        if goodstats:
            opt = True
            res = [bl, bh]
            
    if not opt:
        print ''
        print 'Optimisation failed'
        sys.exit()
            
    return opt, res
    
def save(res):
    
    d = {}
    
    for p, r in res.iteritems():
        
        bl = r[0]
        bh = r[1]
        
        d[p] = {}
        d[p]['allbins'] = np.array([], dtype='float64')

        nb = len(bl)
        
        for i in range(nb):
            
            bname = '_'+p+str(bl[i]).replace('.','p').replace('-','m')+'to'+str(bh[i]).replace('.','p').replace('-','m')
            d[p][bname] = {}
            d[p][bname]['bins'] = [i+1, bl[i], bh[i]]
            
            if p in [options.pv]:
                d[p][bname]['pullx'] = [200, -5.0, 5.0]
                d[p][bname]['pully'] = [200, -5.0, 5.0]
                d[p][bname]['pullz'] = [200, -5.0, 5.0]
                d[p][bname]['resox'] = [200, -300.0, 300.0]
                d[p][bname]['resoy'] = [200, -300.0, 300.0]
                d[p][bname]['resoz'] = [200, -500.0, 500.0]
            else:
                d[p][bname]['d0'] = [200, -1800.0, 1800.0]
                d[p][bname]['dz'] = [200, -3000.0, 3000.0]
    
            d[p]['allbins'] = np.append(d[p]['allbins'], bl[i])
            
            if i == nb-1:
        
                d[p][''] = {}
                d[p]['']['bins'] = [0, bl[0], bh[nb-1]]
                
                if p in [options.pv]:
                    d[p]['']['pullx'] = [200, -5.0, 5.0]
                    d[p]['']['pully'] = [200, -5.0, 5.0]
                    d[p]['']['pullz'] = [200, -5.0, 5.0]
                    d[p]['']['resox'] = [200, -300.0, 300.0]
                    d[p]['']['resoy'] = [200, -300.0, 300.0]
                    d[p]['']['resoz'] = [200, -500.0, 500.0]
                else:
                    d[p]['']['d0'] = [200, -1800.0, 1800.0]
                    d[p]['']['dz'] = [200, -3000.0, 3000.0]
                    
                d[p]['allbins'] = np.append(d[p]['allbins'], bh[i])
    
        d[p]['allbins'] = pd.Series(d[p]['allbins']).to_json(orient='values')

    with open(options.output, 'w') as write_file:
        json.dump(d, write_file, indent=2)
    
if __name__ == '__main__':

    options = main()

    executor = ThreadPoolExecutor(options.threads)

    nbins = int(options.nbins)

    pvparam = options.pv
    
    param = [pvparam]
    param += options.param.split(',')
    
    print 'Parameterisation(s): '+', '.join(param)
    if options.method == 'c': print 'Method: Constant bin width'
    else: print 'Method: Variable bin width'

    flush('Open file: ')
    trk = uproot.open(options.input)[options.tree]
    f = ROOT.TFile(options.input, "READ")
    nev = f.Get(options.tree).GetEntries()
    f.Close()
    ROOT.gROOT.Reset()
    
    crop = int(float(options.crop)*nev)

    nmin = float(options.nmin)
    nminsc = nmin*float(options.crop)
    
    flush('read '+str(crop)+' events (crop factor = '+str(options.crop)+') .. ')

    dtrk = trk.arrays(param, executor=executor, entrystop=crop)

    done()

    flush('Apply quantiles .. ')

    dfc, qup, qcutup, qdown, qcutdown = {}, {}, {}, {}, {}

    for p in param:
        
        if p == pvparam: df = pd.Series(dtrk[p])
        else: df = pd.Series(np.concatenate(dtrk[p]))

        nminopt = nminsc
        if options.meas == 'pv':
            if p == pvparam:
                nminopt = nminsc*float(nbins)
        
        total = df.count()
        qt = nminopt/total
        qup[p] = 1.-qt
        qdown[p] = qt

        if qup[p] < 0.5 or qdown[p] > 0.5:
            print ''
            print '\033[1m'+p+'\033[0m: the requested number of events per fit ('+str(int(nminopt))+') is too high'
            sys.exit()

        qcutup[p] = df.quantile(qup[p])
        qcutdown[p] = df.quantile(qdown[p])
 
        dfc[p] = df[(df < qcutup[p]) & (df > qcutdown[p])]
        
        del df

    done()
    
    for p in param:
        print '->   \033[1m'+p+'\033[0m: events='+str(dfc[p].count()), 'q=['+str(qdown[p])+','+str(qup[p])+']', 'range=['+str(qcutdown[p])+','+str(qcutup[p])+']'

    flush('Optimise binning: ')

    opt = dict.fromkeys(param,False)

    res, bl, bh = {}, {}, {}
        
    for p in param:

        flush(p+' ')
        
        nminopt = nminsc
        if options.meas == 'pv':
            if p == pvparam:
                nminopt = nminsc*float(nbins)
        
        opt[p], res[p] = optimise(dfc[p], nbins, nminopt)
        bl[p] = res[p][0]
        bh[p] = res[p][1]
                
    done()
                
    for p in param:
            
        if not opt[p]:
            print '\033[1m'+p+'\033[0m: bin optimisation failed due to insufficient stats'
            sys.exit()
        else:
            print '->   \033[1m'+p+'\033[0m: nbins='+str(len(bl[p]))
        
    if options.meas == 'bs':
        
        flush('Save bins to file .. ')
        save(res)
        done()
        
    else:
            
        nbpv = len(bl[pvparam])
        
        param.remove(pvparam) 

        flush('Split track data on primary vertex bins .. ')
        
        trk = {}
        for p in param:
            trk[p] = {}
            for i in range(nbpv):
                trk[p][i] = []
        
        for iev, pv in enumerate(dtrk[pvparam]):
            for i in range(nbpv):
                if pv >= bl[pvparam][i] and pv < bh[pvparam][i]:
                    for p in param:
                        for t in dtrk[p][iev]:
                            trk[p][i].append(t)
                            
        del dtrk
        
        done()

        flush('Optimise binning for tracks: ')
        
        opttrk, restrk, bltrk, bhtrk = (dict.fromkeys(param,{}) for _ in range(4))
        
        for k, r in opttrk.iteritems():
            for i in range(nbpv):
                r[i] = False

        for p in param:
            
            flush(p+' ')
            
            restrk[p], bltrk[p], bhtrk[p] = ({} for _ in range(3))
            
            for i in range(nbpv):
                
                dftrk = pd.Series(trk[p][i])
                
                opttrk[p][i], restrk[p][i] = optimise(dftrk, nbins, nminsc)
                bltrk[p][i] = restrk[p][i][0]
                bhtrk[p][i] = restrk[p][i][1]
                
                del dftrk

        done()
        
        minbin = {}
        for p in param:
            
            min = nbins
            for i in range(nbpv):
            
                nbinsc = len(bltrk[p][i])
                if nbinsc <= min:
                    min = nbinsc
                    minbin[p] = i
                    
        for p in param:

            for i in range(nbpv):
        
                if not opttrk[p][i]:
                    print '\033[1m'+p+'\033[0m: bin optimisation failed due to insufficient stats'
                    sys.exit()
                    
            print '->   \033[1m'+p+'\033[0m: pvmin=['+str(bl[pvparam][minbin[p]])+','+str(bh[pvparam][minbin[p]])+']: nbins='+str(len(bltrk[p][minbin[p]]))
        
        flush('Save bins to file .. ')
        save(res)        
        done()
        
