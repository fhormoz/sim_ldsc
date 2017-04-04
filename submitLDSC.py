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
outFolder   = "final_ldsc_result";

if not os.path.exists(outFolder):
        os.makedirs(outFolder);
if not os.path.exists(outFolder + "/TopciseQTL"):
        os.makedirs(outFolder + "/TopciseQTL");

for simFiles in glob.glob(currentPath + "/annots/TopciseQTL/*"):
	print simFiles;
	simFileName = simFiles.replace(currentPath + "/annots/TopciseQTL/", "");
	annotFileName=simFiles + "/" + "TopciseQTL.1.annot"	

	if os.path.exists(currentPath + "/" + outFolder +  "/TopciseQTL/" + simFileName + ".results"):
		continue;
	sumFiles = currentPath + "/simulated_data/SummaryStatistics/"+simFileName + ".sumstats";

	scriptfile = logFolder + "/qsub/LDSC_" + str(simFileName);
        logfile    = logFolder + "/log/LDSC_"  + str(simFileName) + ".log";
        errfile    = logFolder + "/err/LDSC_"  + str(simFileName) + ".err";

	script =  "python " + LDSCPath + "ldsc.py" +\
                        " --h2 " + sumFiles +\
                        " --ref-ld " + annotFileName.replace(".annot","") +\
                        " --chisq-max 9999 --frqfile /groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.1" +\
                        " --w-ld /groups/price/ldsc/reference_files/1000G_EUR_Phase3/weights/weights.hm3_noMHC.1" +\
                        " --overlap-annot --print-cov --print-coefficients --print-delete-vals " +\
                        " --out " + currentPath + "/" + outFolder +  "/TopciseQTL/" + simFileName + "\n";
	print script;	
	scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="LDSC", logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True)
