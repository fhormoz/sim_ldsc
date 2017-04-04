import os
import glob
import pandas as pd
import numpy as npy
import subprocess
from subprocess import Popen, PIPE

def maxSubmitReached(max):
        p1 = Popen(["qstat", "-u", "fhormoz"], stdout=PIPE)
        p2 = Popen(["wc", "-l"], stdin=p1.stdout, stdout=PIPE)
        output = p2.communicate()[0]
        if(int(output) < max) :
                return False;
        else:
                return True;


TEMPLATE_SERIAL = """
#####################################
#!/bin/bash
#BSUB -n 1                #each  job run on 1 core
#BSUB -W 06:00            #job run 2 hour
#BSUB -J {name}
#BSUB -o {logfile}        #lsf output file
#BSUB -e {errfile}        #lsf error file
#BSUB -q short            #submit to "short" queue
#####################################
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
{script}
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
"""

currentPath = os.getcwd();
logFolder   = currentPath + "/log/";
LDSCPath    = "/home/fh80/Code/LDSC/ldsc/" 

for simFiles in glob.glob(currentPath + "/simulated_data/sim2??"):
	script = "rm -rf " + simFiles + "/log/"
	print script;
	script = "rm -rf " + simFiles + "/*log"
	print script;
	script = "rm -rf " + simFiles + "/*nosex";
	print script;