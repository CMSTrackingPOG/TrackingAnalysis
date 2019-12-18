#!/usr/bin/env python

import os
import json

with open('list.json', "r") as read_file:
    flist = json.load(read_file)

if os.path.isdir('jobs'):
    os.system("rm -rf jobs")
os.system("mkdir jobs")

for k, v in flist.iteritems():
    print k
    ijob = 0
    for f in v:
        
        jname = k+'_'+str(ijob)

        fjob = {}
        fjob[k] = f
        
        with open('jobs/'+jname+'.json', "w") as write_file:
            json.dump(fjob, write_file)
        
        jsub = "executable            = jobs/"+jname+".sh\n"
        jsub += "arguments             = --input=jobs/"+jname+".json --output=jobs/"+jname+".root\n"
        jsub += "output                = jobs/"+jname+".out\n"
        jsub += "error                 = jobs/"+jname+".err\n"
        jsub += "log                   = log/"+jname+".log\n"
        jsub += "+JobFlavour           = \"espresso\"\n" # espresso (20min), microcentury (1h), longlunch (2h), workday (8h)
        jsub += "queue"
        
        with open('jobs/'+jname+'.sub', "w") as sub_file:
            sub_file.write(jsub)
            
        jsh = "echo 'hello'"
        
        with open('jobs/'+jname+'.sh', "w") as sh_file:
            sh_file.write(jsh)

        print 'jobs/'+jname+'.sub'
        os.system('condor_submit jobs/'+jname+'.sub')
        
        ijob += 1
