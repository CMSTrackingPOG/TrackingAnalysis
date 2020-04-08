#!/usr/bin/env python

import os
import sys
import json

ntupleProd = 'v20200224'

sys.tracebacklimit = 0

if __name__ == '__main__':

    loc = 'srm://maite.iihe.ac.be:/pnfs/iihe/cms/store/user/kskovpen/Track/Ntuple/Track-'+ntupleProd+'/'

#    fout = open("list.txt","w+")

    files = {}

    dlist = os.popen('gfal-ls '+loc).read().splitlines()

    for i in dlist:

        ds = i
    
        ds = ds.replace('/','')
#        print ds
        ds = ds.replace('-','')
        ds = ds.replace('_','')

#        files[ds] = []
        
#        fout.write(ds+' = [\n')

        flist = os.popen('gfal-ls '+loc+i).read().splitlines()
    
        for j in range(len(flist)):
            
            dss = flist[j]
            dss = dss.replace('/','')
            dss = dss.replace('-','')
            dss = dss.replace('_','')
            print ds+'_'+dss
            
            files[ds+'_'+dss] = []
            
            flistt = os.popen('gfal-ls '+loc+i+'/'+flist[j]).read().splitlines()

            for jj in range(len(flistt)):
        
                flisttt = os.popen('gfal-ls '+loc+i+'/'+flist[j]+'/'+flistt[jj]).read().splitlines()
                
                for jjj in range(len(flisttt)):
                    
                    flistttt = os.popen('gfal-ls '+loc+i+'/'+flist[j]+'/'+flistt[jj]+'/'+flisttt[jjj]).read().splitlines()
                    
                    for jjjj in range(len(flistttt)):
                        
                        if 'log' not in flistttt[jjjj]:
                
#                            fout.write("'"+loc+i+'/'+flist[j]+'/'+flistt[jj]+'/'+flisttt[jjj]+'/'+flistttt[jjjj]+"'")
                            floc = loc+i+'/'+flist[j]+'/'+flistt[jj]+'/'+flisttt[jjj]+'/'+flistttt[jjjj]
                            floc = floc.replace('srm','root')
                            files[ds+'_'+dss].append(floc)

#                        if jjjj != len(flistttt)-1:
#                            fout.write(',\n')
#                        else:
#                            fout.write(']\n\n')

    with open("list.json", "w") as write_file:
        json.dump(files, write_file, indent=2)
