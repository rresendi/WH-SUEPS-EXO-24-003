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
interval = 1  # number of files to be processed in a single job
queue = "workday"
tag = str(sys.argv[1])  # Muon or Electron
NumberOfJobs = -1
doSubmit = True
isdata = sys.argv[2]  # mc or data
eta = sys.argv[3]  # Eta or noEta
era = sys.argv[4]  # 2016, 2016APV, 2017, 2018
mainfolder = sys.argv[6]  # Placeholder, will be set based on the era
output = sys.argv[5]  # "{tag}_efficiencies"
proxy_path = "/tmp/x509up_u168502"
files = []

# Set the appropriate path based on the era
if era == "2016APV":
    mainfolder = "root://cms-xrd-global.cern.ch//store/data/Run2016B/MET/MINIAOD/UL2016_MiniAODv2-v1/"
elif era == "2016":
    mainfolder = "root://cms-xrd-global.cern.ch//store/data/Run2016G/MET/MINIAOD/UL2016_MiniAODv2-v1/"
elif era == "2017":
    mainfolder = "root://cms-xrd-global.cern.ch//store/data/Run2017B/MET/MINIAOD/UL2017_MiniAODv2-v1/"
elif era == "2018":
    mainfolder = "root://cms-xrd-global.cern.ch//store/data/Run2018A/MET/MINIAOD/UL2018_MiniAODv2-v1/"

# Populate the list of files from the xrootd server
# (In practice, you would get these from a DAS query or similar)
# For demonstration, I'll assume a function that lists remote files:
for fil in os.listdir(mainfolder):  # This would be replaced with an appropriate xrootd listing method
    if "root" in fil:
        files.append(os.path.join(mainfolder, fil))

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
