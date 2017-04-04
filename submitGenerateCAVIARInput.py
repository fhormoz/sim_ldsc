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
#BSUB -W 05:00            #job run 5 hour
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

currentPath = os.getcwd();
logFolder   = currentPath + "/log/";

if not os.path.exists(logFolder):
	os.makedirs(logFolder);

for sim in glob.glob(currentPath + "/simulated_data/sim*"):
	simFilName = sim.replace(currentPath + "/simulated_data/", "");
	print simFilName;
	#python generateSimStatAnnot.py -f sim92 -p /groups/price/farhad/Simulation/simulated_data/
	
	scriptfile = logFolder + "/qsub/GeneCAVIAR_" + str(simFilName);
	logfile    = logFolder + "/log/GeneCAVIAR_"  + str(simFilName) + ".log";
	errfile    = logFolder + "/err/GeneCAVIAR_"  + str(simFilName) + ".err";

	script =  "python " + currentPath + "/generateCAVIARInput.py -f " + simFilName + " -p " + currentPath + "/simulated_data/";
	print script;
	scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="GeneStatCAVIAR", logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True);
