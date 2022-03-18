#!/usr/bin/env python3

import os, sys, glob
import htcondor
import json

from optparse import OptionParser

def main(argv = None):
    
    if argv == None:
        argv = sys.argv[1:]
        
    usage = "usage: %prog [options]\n Submit jobs"

    parser = OptionParser(usage)
    parser.add_option("--output", default="jobs", help="Output directory [default: %default]")
    parser.add_option("--splitdata", type=int, default=10, help="Number of files per job [default: %default]")
    parser.add_option("--splitmc", type=int, default=10, help="Number of files per job [default: %default]")
    parser.add_option("--param", default="sumTrackPtSq", help="Parameterisation for PV resolution measurement [default: %default]")
    parser.add_option("--data", action='store_true', help="Only run on data [default: %default]")
    parser.add_option("--mc", action='store_true', help="Only run on mc [default: %default]")
    parser.add_option("--bs", action='store_true', help="Produce trees with only BS information [default: %default]")
    parser.add_option("--pileup", action='store_true', help="do pileup reweighting [default: %default]")
    parser.add_option("--reweight", action='store_true', help="do variable reweighting [default: %default]")
    parser.add_option("--year", default="UL17", help="Year of data taking [default: %default]")
    
    (options, args) = parser.parse_args(sys.argv[1:])
    
    return options

def job(js, output, home, wdir, fsh):

    j = "#!/bin/bash\n\n"

    j += "export X509_USER_PROXY=/user/kskovpen/proxy/x509up_u20657\n"

    j += "source /cvmfs/cms.cern.ch/cmsset_default.sh\n"
    j += "cd "+home+"\n"
    j += "export SCRAM_ARCH=slc7_amd64_gcc820\n"
    j += "eval `scramv1 runtime -sh`;cd -\n"

    j += "cd "+wdir+"\n"

    j += "python "+home+"/./plot.py --input=" + home+"/"+js + " --output=" + home+"/"+output + " --param=" + options.param + " --year=" + options.year + (" --bs" if options.bs else "") + (" --pileup" if options.pileup else "") + (" --reweight" if options.reweight else "") + "\n"

    with open(fsh, 'w') as f:
        f.write(j)
    
    os.system('chmod u+x '+fsh)

if __name__ == '__main__':
    
    options = main()

    home = os.getcwd()
    
    outdir = options.output+options.year
    if os.path.isdir(outdir):
        os.system("rm -rf "+outdir)
    os.system("mkdir "+outdir)

    lszdata = options.splitdata
    lszmc = options.splitmc

    with open("samples"+options.year+".json", "r") as read_file:
        flist = json.load(read_file)
    
    for k, v in flist.items():
        
        if options.data and ('Run20' not in k): continue
        if options.mc and ('Run20' in k): continue
        
        lsz = lszdata if ('Run20' in k) else lszmc
        
        jdir = outdir+'/'+k
        os.system('mkdir '+jdir)
        
        njob = 0
        datasplit = []
        for i in range(len(v)):
            datasplit.append(v[i])
            if len(datasplit) >= lsz or i == len(v)-1:
                with open(jdir+'/'+str(njob) + '.json', 'w') as outfile:
                    data = {}; data[k] = datasplit
                    json.dump(data, outfile, indent=2)
                datasplit = []
                njob += 1
                
    procs = glob.glob(outdir+'/*/')
    
    schedd = htcondor.Schedd()
    njob = 0

    for p in procs:
        
        procName = p.split('/')[-2]
        print(procName)
        jobs = glob.glob(outdir+'/'+procName+'/*.json')
        
        for jidx, j in enumerate(jobs):
            
            jname = j.replace('.json','')
            outname = jname
            outd = '/'.join(j.split('/')[:-1])
            outlog = outname+'.log'
            output = outname+'.root'
            
            job(j, output, home, outd, outname+'.sh')
            
            js = htcondor.Submit({\
            "executable": outname+'.sh', \
            "requirements": '(Machine != \"node25-2.wn.iihe.ac.be\" && Machine != \"node19-8.wn.iihe.ac.be\" && Machine != \"node32-4.wn.iihe.ac.be\" && Machine != \"node22-10.wn.iihe.ac.be\" && Machine != \"node34-10.wn.iihe.ac.be\" && Machine != \"node35-4.wn.iihe.ac.be\" && Machine != \"node19-13.wn.iihe.ac.be\" && Machine != \"node27-3.wn.iihe.ac.be\" && Machine != \"node23-1.wn.iihe.ac.be\")', \
            "output": outname+'.out', \
            "error": outname+'.err', \
            "log": outname+'.log' \
            })
            
            with schedd.transaction() as shd:
                cluster_id = js.queue(shd)
                
            njob += 1
        
    print('Njobs = ', njob)
