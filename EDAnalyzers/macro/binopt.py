#!/bin/env python

from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, sys
import json
import uproot
import ROOT

from optparse import OptionParser

def cmsstyle():
    
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['axes.labelsize'] = 17.0
    plt.rcParams['axes.unicode_minus'] = False
    plt.rcParams['axes.labelpad'] = 10.0
    plt.rcParams['xtick.labelsize'] = 16.0
    plt.rcParams['ytick.labelsize'] = 16.0
    plt.rcParams['legend.fontsize'] = 'small'
    plt.rcParams['legend.handlelength'] = 1.5
    plt.rcParams['legend.borderpad'] = 0.5
    plt.rcParams['legend.frameon'] = True
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['grid.alpha'] = 0.8
    plt.rcParams['grid.linestyle'] = ':'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['savefig.transparent'] = False
    plt.rcParams['figure.subplot.left'] = 0.25
    plt.rcParams['figure.subplot.bottom'] = 0.2
    plt.rcParams['figure.subplot.right'] = 0.96
    plt.rcParams['figure.subplot.top'] = 0.95

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]
        
    usage = "usage: %prog [options]\n Bin optimisation tool"
    
    parser = OptionParser(usage)
    
    parser.add_option("--method", default='v', help="Method for optimisation (variable or constant) [default: %default]")
    parser.add_option("--input", default='JetHT.root', help="Input file name [default: %default]")
    parser.add_option("--output", default='data/bins/qcd_bs.json', help="Output file name [default: %default]")
    parser.add_option("--tree", default='trackTree', help="Input tree name [default: %default]")
    parser.add_option("--nmin", type=int, default=40000, help="Minimum number of events per bin in the first-level parameterisation [default: %default]")
    parser.add_option("--nbins", type=int, default=50, help="Maximum number of bins [default: %default]")
    parser.add_option("--threads", default=8, help="Number of threads [default: %default]")
    parser.add_option("--crop", type=float, default=0.1, help="Crop factor [default: %default]")
    parser.add_option("--quantile", type=float, default=0.01, help="Fraction of events to cut from the sides [default: %default]")
    parser.add_option("--meas", default='bs', help="Measurement type (bs or pv) [default: %default]")
    parser.add_option("--pv", default='sumTrackPtSq', help="PV parameterisation [default: %default]")
    parser.add_option("--param", default='pt,eta,phi,npv,dr', help="List of track parameterisations [default: %default]")
    parser.add_option("--plot", action='store_true', help="Draw validation plots [default: %default]")
    parser.add_option("--validation", action='store_true', help="Fill validation histograms for multidimensional binning [default: %default]")
    parser.add_option("--symmetric", action='store_true', help="Use symmetric binning for eta and phi parameterisations [default: %default]")

    (options, args) = parser.parse_args(sys.argv[1:])

    return options

def done():

    sys.stdout.write('\x1b[\33[32mdone\x1b[0m\n'); sys.stdout.flush()

def flush(s, flush=True):

    sys.stdout.write(s)
    if flush: sys.stdout.flush()
    
def plot(x, bins, p, pref=''):

    pdir = 'pics/'
    
    if not os.path.isdir(pdir):
        os.system("mkdir "+pdir)
    
    fig = plt.hist(x, bins, facecolor='royalblue', alpha=0.5)
    
    for b in bins:
        plt.axvline(x=b, color='r', linestyle='dashed', linewidth=1)            
        
    plt.xlabel(p)
    plt.ylabel('Events')
    
#    if 'pt' in p.lower() or 'dr' in p.lower(): plt.yscale('log')
#    if 'pt' in p.lower(): plt.yscale('log')
    if 'pt' in p: plt.yscale('log')

    plt.savefig('pics/'+p+pref+'.eps')
    plt.close()

def lastbin(data, nbins, nmin, increase = True, rebin = False):

    nbinmax = 10000

    opt = False
    
    res, bl, bh, bli, bhi = [], [], [], [], []
    
    x = np.array(data.values.tolist())
    
    fig, hbins, hpatches = plt.hist(x, nbinmax)
    
    lbin, ibin = 0, -1
    
    for ib in reversed(range(len(hbins)-1)):
        
        lbin += fig[ib]
        ibin = ib
        
        dltlb = len(hbins)-ibin-1
        nw = len(hbins)/dltlb
        
        if (lbin > nmin) and (nw <= nbins):
            
            dltlbmod = dltlb
            ibinlb = ibin
            
            mergelb = False
            while True:
                
                r = float(len(hbins)-dltlbmod) % float(dltlb)
                if r == 0:
                    if ibin > ibinlb+dltlb/2:
                        mergelb = True
                    break
                else:
                    dltlbmod -= 1
                    ibin += 1
            
            nbinsmain = (len(hbins)-dltlbmod)/dltlb

            for ibb in range(nbinsmain):
                
                bl.append(hbins[ibb*dltlb])
                bli.append(ibb*dltlb)
                if ibb == nbinsmain-1 and mergelb: break
                bh.append(hbins[(ibb+1)*dltlb])
                bhi.append((ibb+1)*dltlb)

            if not mergelb:
                bl.append(hbins[ibin])
                bli.append(ibin)
            bh.append(hbins[-1])
            bhi.append(nbinmax)
            
            break

    lbsum = 0
    for ib in range(bli[-1], bhi[-1]): lbsum += fig[ib]
    lpsum = 0
    for ib in range(bli[-2], bhi[-2]): lpsum += fig[ib]

    blf = bl
    bhf = bh

    if rebin:
        
        if lpsum < lbsum:
        
            blf, bhf = [], []
        
            for ib in range(len(bl)):
                r = ib % 2
                if r == 0:
                    blf.append(bl[ib])
                else:
                    bhf.append(bh[ib])
                if ib == len(bl)-1:
                    bhf.append(bh[ib])
                
#    for i in range(len(blf)):
#        print '[', blf[i], bhf[i], ']', 'size=', bhf[i]-blf[i]
        
    if (ibin < nbins) or (lbin < nmin):
        print ''
        print 'Requested number of bins or last bin stats is too big'
        sys.exit()

    opt = True
    res = [blf, bhf]
            
    if not opt:
        print ''
        print 'Optimisation failed'
        sys.exit()

    plt.close()
        
    return opt, res
    
def optimise(data, nbins, nmin, var = ''):
    
    nbinsc = nbins
    
    opt = False; res = []
            
    while not opt:
             
        bl = []
        bh = []
        
        if options.method == 'c' or (var in ['pt']):
            print 'constant for pt'
            nmin = 1.
            nbinsc = 40
            splt = pd.cut(data, nbinsc, duplicates='drop') # fixme
        else: splt = pd.qcut(data, nbinsc, duplicates='drop')
        bins = splt.value_counts(sort=False)
#        print bins
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

def save(res, dfc, dfcpv = None):
    
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
                
                d[p][bname]['pullx'] = [600, -5.0, 5.0]
                d[p][bname]['pully'] = [600, -5.0, 5.0]
                d[p][bname]['pullz'] = [600, -5.0, 5.0]
                d[p][bname]['resox'] = [600, -300.0, 300.0]
                d[p][bname]['resoy'] = [600, -300.0, 300.0]
                d[p][bname]['resoz'] = [600, -500.0, 500.0]
                
            else:

                if options.meas == 'pv':
                    d[p][bname]['d0'] = [800, -1800.0, 1800.0]
                else:
                    d[p][bname]['d0'] = [1800, -1800.0, 1800.0]
                
                if p in ['eta']:
                    if abs(bl[i]) < 1.5 and abs(bh[i]) < 1.5:
                        d[p][bname]['dz'] = [600, -1500.0, 1500.0]
                    elif abs(bl[i]) < 2.3 and abs(bh[i]) < 2.3:
                        d[p][bname]['dz'] = [600, -3000.0, 3000.0]
                    elif abs(bl[i]) < 2.6 and abs(bh[i]) < 2.6:
                        d[p][bname]['dz'] = [1200, -6000.0, 6000.0]
                    else:
                        d[p][bname]['dz'] = [1200, -12000.0, 12000.0]
                else:
                    d[p][bname]['dz'] = [600, -3000.0, 3000.0]
    
            d[p]['allbins'] = np.append(d[p]['allbins'], bl[i])
            
            if i == nb-1:
        
                d[p][''] = {}
                d[p]['']['bins'] = [0, bl[0], bh[nb-1]]
                
                if p in [options.pv]:
                    
                    d[p]['']['pullx'] = [600, -5.0, 5.0]
                    d[p]['']['pully'] = [600, -5.0, 5.0]
                    d[p]['']['pullz'] = [600, -5.0, 5.0]
                    d[p]['']['resox'] = [600, -300.0, 300.0]
                    d[p]['']['resoy'] = [600, -300.0, 300.0]
                    d[p]['']['resoz'] = [600, -500.0, 500.0]
                    
                else:
                    
                    if options.meas == 'pv':
                        d[p][bname]['d0'] = [800, -1800.0, 1800.0]
                    else:
                        d[p][bname]['d0'] = [1800, -1800.0, 1800.0]
                
                    if p in ['eta']:
                        if abs(bl[i]) < 1.5 and abs(bh[i]) < 1.5:
                            d[p][bname]['dz'] = [600, -1500.0, 1500.0]
                        elif abs(bl[i]) < 2.3 and abs(bh[i]) < 2.3:
                            d[p][bname]['dz'] = [600, -3000.0, 3000.0]
                        elif abs(bl[i]) < 2.6 and abs(bh[i]) < 2.6:
                            d[p][bname]['dz'] = [1200, -6000.0, 6000.0]
                        else:
                            d[p][bname]['dz'] = [1200, -12000.0, 12000.0]
                    else:
                        d[p][bname]['dz'] = [600, -3000.0, 3000.0]
                    
                d[p]['allbins'] = np.append(d[p]['allbins'], bh[i])
                
        bins = d[p]['allbins']
        d[p]['allbins'] = pd.Series(d[p]['allbins']).to_json(orient='values')

        if dfc != None:
            
            if options.plot and options.meas == 'bs':
                x = np.array(dfc[p].values.tolist())
                plot(x, bins, p)
            elif options.plot and options.meas == 'pv':
                if p not in options.pv:
                    for k, v in dfc[p].iteritems():
                        x = np.array(v.values.tolist())
                        plot(x, bins, p, str(k))
                else:
                    x = np.array(dfcpv[p].values.tolist())
                    plot(x, bins, p)

    with open(options.output, 'w') as write_file:
        json.dump(d, write_file, indent=2)
    
if __name__ == '__main__':

    options = main()
    
    cmsstyle()

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
    
    flush('use '+str(crop)+' events (crop factor = '+str(options.crop)+') .. ')

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
        qt = nminopt/(total*float(options.crop))
#        if 'pt' in p.lower() or 'dr' in p.lower(): qt = max(qt, options.quantile)
#        if 'pt' in p.lower(): qt = max(qt, options.quantile)
##        if 'pt' in p: qt = max(qt, options.quantile)
        
        qup[p] = 1.-qt
        qdown[p] = qt

        if qup[p] < 0.5 or qdown[p] > 0.5:
            print ''
            print '\033[1m'+p+'\033[0m: the requested number of events per fit ('+str(int(nminopt))+') is too high'
            sys.exit()
            
        # do not apply quantiles in case of angular parameterisations
        if p in ['eta', 'phi']:
            qup[p] = 1.
            qdown[p] = 0.

        qcutup[p] = df.quantile(qup[p])
        qcutdown[p] = df.quantile(qdown[p])
        
        if p in ['pt']: 
            qcutup[p] = 10.0
            qcutdown[p] = 0.4
 
        dfc[p] = df[(df < qcutup[p]) & (df > qcutdown[p])]

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

#        if (options.method in ['c', 'v']) and ('pt' not in p.lower() and 'dr' not in p.lower()):
#        if (options.method in ['c', 'v']) and ('pt' not in p.lower()):
        if (options.method in ['c', 'v']) and ('pt' not in p):
            opt[p], res[p] = optimise(dfc[p], nbins, nminopt, p)
        else:
            opt[p], res[p] = optimise(dfc[p], nbins, nminopt, p)
#            opt[p], res[p] = lastbin(dfc[p], nbins, nminopt)
            
        if options.symmetric and p in ['eta', 'phi']:
                         
            sbins = []

            for b in range(len(res[p][0])):
                sbins.append(res[p][0][b])
            sbins.append(res[p][1][len(res[p][0])-1])
                
            spos = -1
            sneg = -1
            for b in range(len(sbins)):
                if sbins[b] > 0:
                    spos = b
                    sneg = b-1
                    break

            binssym = []
            
            if sbins[spos] > abs(sbins[sneg]):
                del sbins[0:sneg+1]
                neg = [-b for b in reversed(sbins)]
                binssym += neg
                binssym.append(0.)
                binssym += sbins
            else: 
                del sbins[spos-1:]
                pos = [-b for b in reversed(sbins)]
                binssym += sbins
                binssym.append(0.)
                binssym += pos
                
            res[p][0] = []
            res[p][1] = []

            b = 0
            while b < len(binssym)-1:
                res[p][0].append(binssym[b])
                res[p][1].append(binssym[b+1])
                b += 1

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
        
        flush('Save results .. ')
        save(res, dfc)
        done()
        
    else:
            
        nbpv = len(bl[pvparam])
        
        param.remove(pvparam)
        
        dftrk = None
        
        if options.validation:

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
        
            flush('Create data sets with tracks using final binning: ')
        
            opttrk, restrk, bltrk, bhtrk = (dict.fromkeys(param,{}) for _ in range(4))
        
            for k, r in opttrk.iteritems():
                for i in range(nbpv):
                    r[i] = False

            dftrk = {}
            for p in param:
            
                flush(p+' ')
            
                restrk[p], bltrk[p], bhtrk[p] = ({} for _ in range(3))
            
                dftrk[p] = {}
                for i in range(nbpv):
                
                    dftrk[p][i] = pd.Series(trk[p][i])
                    if p in ['dr']: dftrk[p][i] = dftrk[p][i].drop_duplicates()

            done()
        
        flush('Save results .. ')
        save(res, dftrk, dfc)
        done()
        
