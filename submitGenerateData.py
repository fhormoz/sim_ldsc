import os
import glob
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
#BSUB -W 06:00            #job run 5 hour
#BSUB -J {name}
#BSUB -o {logfile}        #lsf output file
#BSUB -e {errfile}        #lsf error file
#BSUB -q short         	  #submit to "short" queue
#####################################
echo "------------------------------------------------------------------------"
echo "Job started on" `date`
echo "------------------------------------------------------------------------"
{script}
echo "------------------------------------------------------------------------"
echo "Job ended on" `date`
echo "------------------------------------------------------------------------"
"""

Nsim = 1000;
currentPath = os.getcwd();
outFolder   = currentPath + "/simulated_data";
logFolder   = currentPath + "/log/";

if not os.path.exists(outFolder):
	os.makedirs(outFolder);
if not os.path.exists(logFolder):
	os.makedirs(logFolder);

if not os.path.exists(logFolder+"/qsub/"):
        os.makedirs(logFolder+"/qsub");
if not os.path.exists(logFolder+"/log/"):
        os.makedirs(logFolder+"/log/");
if not os.path.exists(logFolder+"/err/"):
        os.makedirs(logFolder+"/err/");

for sim in range(Nsim):
	if os.path.exists(outFolder +  "/sim" + str(sim)):
		continue;
		
	scriptfile = logFolder + "/qsub/sim_" + str(sim);
	logfile    = logFolder + "/log/sim_"  + str(sim) + ".log";
	errfile    = logFolder + "/err/sim_"  + str(sim) + ".err";

	script =  "python " + currentPath + "/generateData.py -s " + str(sim) + " -o simulated_data";
	print script;
	scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="simGenData", logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True);
