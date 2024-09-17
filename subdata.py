#!/usr/bin/env python3
import os
import re
import math
import time
import sys
print()
print('START')
print()
########   YOU ONLY NEED TO FILL THE AREA BELOW   #########
########   customization area #########
interval = 1  # number files to be processed in a single job, take care to split your file so that you run on all files. The last job might be with smaller number of files (the ones that remain).
queue = "workday"
tag = str(sys.argv[1]) #Muon or Electron
NumberOfJobs = -1
doSubmit = False
isdata = sys.argv[2]  # mc or data
eta = sys.argv[3] # Eta or noEta
era = sys.argv[4] # 2016, 2016APV, 2017, 2018
file_list = sys.argv[6]
output = sys.argv[5]  # "{tag}_efficiencies"
proxy_path = "/afs/cern.ch/user/r/rresendi/private/x509up_u168502"
files = []

with open(file_list, 'r') as f:
    for line in f:
        files.append(line.strip())
        
########   customization end   #########
path = os.getcwd()
print()
print('do not worry about folder creation:')
os.system(f"rm -rf tmp{tag}")
os.system(f"rm -rf exec{tag}")
os.system(f"rm -rf batchlogs{tag}")
os.system(f"mkdir -p tmp{tag}")
os.system(f"mkdir -p exec{tag}")
print()

if NumberOfJobs == -1:
    NumberOfJobs = (len(files) + interval) // interval
##### loop for creating and sending jobs #####
for x in range(1, int(NumberOfJobs) + 1):
    ##### creates directory and file list for job #######
    jobFiles = files[max(0, (x - 1) * interval):min(x * interval, len(files))]
    with open(f'exec{tag}/job_{x}.sh', 'w') as fout:
        fout.write("#!/bin/sh\n")
        fout.write("echo\n")
        fout.write("echo\n")
        fout.write("echo 'START---------------'\n")
        fout.write("echo 'WORKDIR ' ${PWD}\n")
        fout.write("export HOME=$PWD\n")
        fout.write("export X509_USER_PROXY=${proxy_path}\n")
        fout.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
        fout.write("cmsrel CMSSW_13_3_3\n")
        fout.write("cd CMSSW_13_3_3/src\n")
        fout.write("cmsenv\n")
        fout.write("cd -\n")
        fout.write(f"python3 /afs/cern.ch/user/r/rresendi/data.py {tag} {isdata} {eta} {era} output_{x}.root {' '.join(jobFiles)}\n")
        fout.write(f"echo 'cp output_{x}.root /afs/cern.ch/user/r/rresendi/{output}/output_{x}.root'\n")
        fout.write(f"cp output_{x}.root /afs/cern.ch/user/r/rresendi/{output}/output_{x}.root\n")
        fout.write("echo 'STOP---------------'\n")
        fout.write("echo\n")
        fout.write("echo\n")
    os.system(f"chmod 755 exec{tag}/job_{x}.sh")

###### create submit.sub file ####
os.mkdir(f"batchlogs{tag}")
with open('submit.sub', 'w') as fout:
    fout.write("executable              = $(filename)\n")
    fout.write("arguments               = $(Proxy_path) $(ClusterId)$(ProcId)\n")
    fout.write(f"output                  = batchlogs{tag}/$(ClusterId).$(ProcId).out\n")
    fout.write(f"error                   = batchlogs{tag}/$(ClusterId).$(ProcId).err\n")
    fout.write(f"log                     = batchlogs{tag}/$(ClusterId).log\n")
    fout.write(f'+JobFlavour = "{queue}"\n')
    fout.write(f"x509userproxy = {proxy_path}\n")
    fout.write("\n")
    fout.write(f"queue filename matching (exec{tag}/job_*sh)\n")

###### sends bjobs ######
os.system("echo submit.sub")
if doSubmit:
    os.system("condor_submit -spool submit.sub")

print()
print("your jobs:")
os.system("condor_q")
print()
print('END')
