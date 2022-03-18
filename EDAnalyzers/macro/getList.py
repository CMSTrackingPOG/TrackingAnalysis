#!/usr/bin/env python3

import os
import sys
import json

ntupleProd = 'v20220314'
year = 'UL17'

sys.tracebacklimit = 0

if __name__ == '__main__':

    lsc = 'xrdfs maite.iihe.ac.be ls ' 
    
    loc = '/pnfs/iihe/cms/store/user/kskovpen/Track/Ntuple/Track-'+ntupleProd+'/'

    files = {}

    dlist = os.popen(lsc+loc).read().splitlines()

    for i in dlist:

        ds = i.split('/')[-1]
    
        flist = os.popen(lsc+i).read().splitlines()

        for j in range(len(flist)):

            if year not in flist[j] and year.replace('17', '2017') not in flist[j]: continue

            dss = flist[j].split('/')[-1]

            dname = ds+'_'+dss
            print(dname)

            # do not run on samples with BS constraint
            if 'withBS' in dname: continue
            
            files[ds+'_'+dss] = []

            flistt = os.popen(lsc+flist[j]).read().splitlines()

            for jj in range(len(flistt)):
        
                flisttt = os.popen(lsc+flistt[jj]).read().splitlines()
                
                for jjj in range(len(flisttt)):
                    
                    flistttt = os.popen(lsc+flisttt[jjj]).read().splitlines()
                    
                    for jjjj in range(len(flistttt)):
                        
                        if 'log' not in flistttt[jjjj]:
                
                            floc = 'root://maite.iihe.ac.be/'+flistttt[jjjj]
                            files[ds+'_'+dss].append(floc)

    with open("samples"+year+".json", "w") as write_file:
        json.dump(files, write_file, indent=2)
