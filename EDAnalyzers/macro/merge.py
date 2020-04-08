#!/usr/bin/env python

import os, sys, glob

from optparse import OptionParser

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]
        
    usage = "usage: %prog [options]\n Merge the output of processed jobs"
    
    parser = OptionParser(usage)
    parser.add_option("-d","--dir",default="jobs/",help="directory with the plot results [default: %default]")

    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

if __name__ == '__main__':
    
    options = main()

    cwd = os.getcwd()
    
    dir = options.dir
    if not os.path.isdir(dir):
        print 'Directory '+dir+'/ does not exist'
        sys.exit()

    jobs = glob.glob('jobs/list_*.root')
    
    fd = {}
    
    for j in jobs:
 
        jname = j.replace('.root','').replace(dir+'list_','')
        samp = jname.split('_')[0]+'_'+jname.split('_')[1]
        
        if samp in fd.keys():
            fd[samp].append(j)
        else:
            fd[samp] = []

    for k,v in fd.items():
        
        flist = " ".join(fd[k])
        os.system('hadd -f jobs/'+k+'.root '+flist)
        os.system('rm jobs/list_'+k+'_*.root')
