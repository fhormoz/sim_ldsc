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

for simFiles in glob.glob(currentPath + "/annots/TopciseQTL/*"):
	simFileName = simFiles.replace(currentPath + "/annots/TopciseQTL/", "");
	annotFileName=simFiles + "/" + "TopciseQTL.1.annot"	

	if (os.path.exists(simFiles + "/TopciseQTL.1.l2.ldscore.gz")):
		continue;
	else :
		print simFiles;
	
	scriptfile = logFolder + "/qsub/GeneStat_" + str(simFileName);
        logfile    = logFolder + "/log/GeneStat_"  + str(simFileName) + ".log";
        errfile    = logFolder + "/err/GeneStat_"  + str(simFileName) + ".err";

	script = "python " + LDSCPath+ "ldsc.py" +\
		" --l2 --bfile /n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000_HapMapCommon/UKBioBank_chr_500ind_eQTL.CM " +\
		" --ld-wind-cm 1 --annot " + annotFileName +\
		" --out " + annotFileName.replace(".annot","") +\
		" --print-snps /home/fh80/Code/RunLDSC/list.txt";
	scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="LDSC", logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True)
