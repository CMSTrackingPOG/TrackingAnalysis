#!/bin/env python

from concurrent.futures import ThreadPoolExecutor
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import OptimizeWarning
import warnings
import os, sys
import json
import uproot
import ROOT

from optparse import OptionParser

def cmsstyle():
    
#    plt.rcParams['font.family'] = 'sans-serif'
#    plt.rcParams['font.sans-serif'] = 'Helvetica'
#    plt.rcParams['mathtext.fontset'] = 'custom'
#    plt.rcParams['mathtext.rm'] = 'Helvetica'
#    plt.rcParams['mathtext.bf'] = 'Helvetica:bold'
#    plt.rcParams['mathtext.sf'] = 'Helvetica'
#    plt.rcParams['mathtext.it'] = 'Helvetica:italic'
#    plt.rcParams['mathtext.tt'] = 'Helvetica'
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
#    plt.rcParams['xtick.major.size'] = 12
#    plt.rcParams['xtick.minor.size'] = 6
#    plt.rcParams['xtick.major.pad'] = 6
#    plt.rcParams['xtick.top'] = False
#    plt.rcParams['xtick.major.top'] = False
#    plt.rcParams['xtick.major.bottom'] = True
#    plt.rcParams['xtick.minor.top'] = False
#    plt.rcParams['xtick.minor.bottom'] = True
#    plt.rcParams['xtick.minor.visible'] = True
#    plt.rcParams['ytick.direction'] = 'in'
#    plt.rcParams['ytick.major.size'] = 12
#    plt.rcParams['ytick.minor.size'] = 6.0
#    plt.rcParams['ytick.right'] = True
#    plt.rcParams['ytick.major.left'] = False
#    plt.rcParams['ytick.major.right'] = False
#    plt.rcParams['ytick.minor.left'] = False
#    plt.rcParams['ytick.minor.right'] = False
#    plt.rcParams['ytick.minor.visible'] = False
    plt.rcParams['grid.alpha'] = 0.8
    plt.rcParams['grid.linestyle'] = ':'
    plt.rcParams['axes.linewidth'] = 1
    plt.rcParams['savefig.transparent'] = False
    plt.rcParams['figure.subplot.left'] = 0.1
    plt.rcParams['figure.subplot.bottom'] = 0.2
    plt.rcParams['figure.subplot.right'] = 0.96
    plt.rcParams['figure.subplot.top'] = 0.95

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]
        
    usage = "usage: %prog [options]\n Beam spot monitoring"
    
    parser = OptionParser(usage)
    
    parser.add_option("--input",default='JetHT.root',help="Input file name [default: %default]")
    parser.add_option("--output",default='data/bins/qcd_bsw.json',help="Output file name [default: %default]")
    parser.add_option("--tree",default='trackTree',help="Input tree name [default: %default]")
    parser.add_option("--lumi",default='../crab/lumi.csv',help="Luminosity data file name [default: %default]")
    parser.add_option("--threads",default=8,help="Number of threads [default: %default]")
    parser.add_option("--crop",default=1.,help="Crop factor [default: %default]")
    parser.add_option("--cut",action='store_true',help="Apply selection [default: %default]")
    parser.add_option("--weight",action='store_true',help="Apply lumi weights [default: %default]")
    parser.add_option("--qcd",action='store_true',help="Use QCD events [default: %default]")
    parser.add_option("--draw",action='store_true',help="Draw bin regions [default: %default]")

    (options, args) = parser.parse_args(sys.argv[1:])

    return options

def done():
    sys.stdout.write('\x1b[\33[32mdone\x1b[0m\n'); sys.stdout.flush()

def flush(s, flush=True):
    sys.stdout.write(s)
    if flush: sys.stdout.flush()

def lumi(fname):    
    lm = pd.read_csv(fname, low_memory=True, engine='python', sep=',|:', comment='#', names=['run','fill','ls1','ls2','time1','time2','time3','beamstatus','E','delivered','recorded','avgpu','source'])
    return lm
    
def fitf(x, a, b):
    return b
    
if __name__ == '__main__':

    options = main()
    
    cmsstyle()

    warnings.simplefilter("ignore", OptimizeWarning)
    pd.options.mode.chained_assignment = None
#    pd.set_option('display.max_rows', -1)
    
    executor = ThreadPoolExecutor(options.threads)

    if options.weight:        
        flush('Read lumi summary: ')    
        lm = lumi(options.lumi)        
        done()
    
    flush('Open data file: ')
    
    tbs = uproot.open(options.input)[options.tree]
    f = ROOT.TFile(options.input, "READ")
    nev = f.Get(options.tree).GetEntries()
    f.Close()
    ROOT.gROOT.Reset()
    
    flush('read '+str(nev)+' .. ')

    var = ['run', 'lumi', 'beamWidthX', 'beamWidthY', 'beamSigmaZ', 'beamWidthXError', 'beamWidthYError', 'beamSigmaZError']
    titl = ['run', 'lumi', r'Beam width in x [$\mu m$]', r'Beam width in y [$\mu m$]', r'Beam $\sigma_{z}$ [mm]', \
    r'Beam width error in x [$\mu m$]', r'Beam width error in y [$\mu m$]', r'Beam $\sigma_{z}$ error [mm]']

    dbs = tbs.pandas.df(var, executor=executor, entrystop=nev*float(options.crop))
    dbs = dbs.drop_duplicates().sort_values(['run', 'lumi'], ascending = (True, True))
    
    bvar = ['beamWidthX', 'beamWidthY', 'beamWidthXError', 'beamWidthYError']
    dbs[bvar] = dbs[bvar].apply(lambda x: x*10000, axis=1)
    
#    if options.cut: dbs = dbs[(dbs.beamWidthXError < 0.1) & (dbs.beamWidthYError < 0.1) & (dbs.beamSigmaZError < 0.1)]
    if options.cut: dbs = dbs[dbs.beamWidthXError < 0.1]
    
    # remove lb dependence for the final parameterisation
    dbspl = dbs.drop_duplicates(subset=['run', 'beamWidthX']).sort_values(['run', 'lumi'], ascending = (True, True)).reset_index(drop=True)
    dbspl['bin'] = dbspl.index

    if options.weight:
        # keep lb splitting for luminosity weighting
        dbs = dbs.reset_index(drop=True)
        dbs['bin'] = dbs.index

    done()
    
    if options.weight:
        
        flush('Assign lumi weights: ')
    
        run = dbs['run'].values.tolist()
        lumi = dbs['lumi'].values.tolist()
        width = dbs['beamWidthX'].values.tolist()
        
        il = {} # as function of width
        lumib, lumie = ([] for _ in range(2))
        
        wpred = str(run[0])+'_'+str(width[0])
        lpred = lumi[0]
    
        for ir, r in enumerate(run):
        
            w = width[ir]
            l = lumi[ir]
        
            lmrun = lm['run']
            lmls1 = lm['ls1']
        
            k = str(r)+'_'+str(w)
        
            if k != wpred: lumie.append(lpred)
            if ir == len(run)-1: lumie.append(l)
        
            if k not in il:
                il[k] = []
                lumib.append(l)
        
            res = lm.loc[(lmrun == r) & (lmls1 == l)]['recorded']
        
            if len(res) != 1:
                print ''
                print 'Multiple matches found'
                sys.exit()            
            
            il[k].append(float(res))
        
            wpred = k
            lpred = l
        
        ilrun = []
        run = dbspl['run'].values.tolist()
        width = dbspl['beamWidthX'].values.tolist()
        for i, w in enumerate(width):
            tot = sum(il[str(run[i])+'_'+str(w)])
            ilrun.append(1./tot)
        iltot = sum(ilrun)
        ilrunnorm = [x/iltot for x in ilrun]
    
        dbspl['weight'] = pd.Series(ilrunnorm, index=dbspl.index)
        dbspl['lumib'] = pd.Series(lumib, index=dbspl.index)
        dbspl['lumie'] = pd.Series(lumie, index=dbspl.index)
            
        done()
    
    flush('Produce plots .. ')

    figsize = [15.0, 7.0]
    ylim = {'beamWidth':[0., 35.], 'beamWidthX':[0., 35.], 'beamWidthY':[0., 35.], \
    'beamWidthError':[0., 2.0], 'beamWidthXError':[0., 2.0], 'beamWidthYError':[0., 2.0], \
    'beamSigmaZ':[30., 50.], 'beamSigmaZError':[0., 1.]}

    # final binning
    if options.qcd:
        roi = [0, 150, 385, 920, dbspl['bin'].iloc[-1]+1]
    else:
        roi = [0, 140, 870, 1310, 1830, dbspl['bin'].iloc[-1]+1]

    plots = {}

    for ix, x in enumerate(['run', 'lumi']):
        var.remove(x)
        titl.remove(x)
    
    for i, v in enumerate(var):
        plots[v] = dbspl.plot(kind='scatter', x='bin', y=v, color='black', figsize=figsize)
        plots[v].set_xlabel(r'Time [Arbitrary units]')
        plots[v].set_ylabel(titl[i])
        if v in ylim: plots[v].set_ylim(ylim[v][0], ylim[v][1])
        if options.draw:
            for r in roi:
                plots[v].axvline(x=r-0.5, ymin=ylim[v][0], ymax=ylim[v][1], linestyle='--', color='gray', linewidth=2.0)
        fig = plots[v].get_figure()
        fig.savefig('pics/'+v+'.pdf')

    cols = ['lightskyblue', 'salmon']
    fcols = ['blue', 'red']
    labs = ['x', 'y']
    name = ['beamWidth', 'beamWidthError']
    titl =  [r'Beam width [$\mu m$]', r'Beam width error [$\mu m$]']
    same = [['beamWidthX', 'beamWidthY'], ['beamWidthXError', 'beamWidthYError']]
    
    param = {}
    columns = ['runstart', 'runend', 'lumistart', 'lumiend', 'beamwidthx', 'beamwidthy']
    for c in columns: param[c] = []
    
    for i, s in enumerate(same):
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        for j, v in enumerate(s):
            dbspl.plot(kind='scatter', x='bin', y=v, color=cols[j], label=labs[j], figsize=figsize, ax=ax)
            if v in ['beamWidthX', 'beamWidthY']:
                
                y = dbspl[v].values.tolist()
                x = dbspl['bin'].values.tolist()
                if options.weight: w = dbspl['weight'].values.tolist()
                
                for r in range(len(roi)-1):
                    xsel, ysel, wsel = ([] for _ in range(3))
                    for e in range(len(x)):
                        if roi[r] <= x[e] < roi[r+1]:
                            xsel.append(x[e])
                            ysel.append(y[e])
                            if options.weight: wsel.append(w[e])
                    if options.weight: par = curve_fit(fitf, xsel, ysel, sigma=wsel, absolute_sigma=False)[0]
                    else: par = curve_fit(fitf, xsel, ysel)[0]
                    yf = [fitf(d, par[0], par[1]) for d in xsel]
                    ax.plot(xsel, yf, 'r', color=fcols[j], linewidth=3)
                    ytext = 2. if v == 'beamWidthX' else 4.
                    vname = 'x' if v == 'beamWidthX' else 'y'
                    if options.draw:
                        ax.text((xsel[0]+xsel[-1])/2.-100., ytext, r''+vname+' = '+format(par[1], '.2f')+' $\mu m$', fontsize=12, color=fcols[j])
                    
                    istart = roi[r]
                    iend =  roi[r+1]-1
                    
                    if iend > len(dbspl.index):
                        print ''
                        print 'Requested bin edge exceeds the derived parametrisation'
                        sys.exit()

                    if v == 'beamWidthX':
                        param['runstart'].append(dbspl[dbspl['bin'] == istart].run.values.tolist()[0])
                        param['runend'].append(dbspl[dbspl['bin'] == iend].run.values.tolist()[0])
                        param['lumistart'].append(dbspl[dbspl['bin'] == istart].lumi.values.tolist()[0])
                        param['lumiend'].append(dbspl[dbspl['bin'] == iend].lumi.values.tolist()[0])
                        param['beamwidthx'].append(par[1])
                    else:
                        param['beamwidthy'].append(par[1])
                    
        ax.set_xlabel(r'Time [Arbitrary units]')
        ax.set_ylabel(titl[i])
        if v in ylim: ax.set_ylim(ylim[v][0], ylim[v][1])
        if options.draw:
            for r in roi:
                ax.axvline(x=r-0.5, ymin=ylim[v][0], ymax=ylim[v][1], linestyle='--', color='gray', linewidth=2.0)
        ax.legend(loc='upper center', prop={'size': 14})
        fig.savefig('pics/'+name[i]+'.pdf')
    
    done()

    flush('Save to file .. ')
    
    with open(options.output, 'w') as json_file:
        json.dump(param, json_file, indent=2)

    done()
