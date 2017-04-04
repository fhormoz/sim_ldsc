import os
import glob
import subprocess
from subprocess import Popen, PIPE


def maxSubmitReached(max):
        p1 = Popen(["bjobs", "-u", "fh80"], stdout=PIPE)
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

script = "";
currentPath = os.getcwd();
logFolder   = currentPath + "/log/";

for sim in glob.glob(currentPath + "/simulated_data/sim*"):
	print sim.replace(currentPath + "/simulated_data/", "");
	scriptfile = logFolder + "/qsub/sim_" + sim.replace(currentPath + "/simulated_data/", "");
	logfile    = logFolder + "/log/sim_"  + sim.replace(currentPath + "/simulated_data/", "") + ".log";
	errfile    = logFolder + "/err/sim_"  + sim.replace(currentPath + "/simulated_data/", "") + ".err";
	outFolder  = sim + "/caviar_output/";
	for caviarInput in glob.glob(sim + "/caviar_input/*.Z"):
		zFile = caviarInput;
		ldFile = caviarInput.replace(".Z", ".ld");
		outputFileName = caviarInput.replace(sim + "/caviar_input/", "").replace(".Z","")
		outputfile = outFolder + outputFileName ;
		if not os.path.exists(outFolder):
			os.makedirs(outFolder); 
		if(os.path.exists(outputfile+"_post")):
			continue;		
		script = script + "/home/fh80/Code/CAVIAR/caviar/CAVIAR-C++/CAVIAR " +\
			" -o " + outputfile +\
			" -z " + zFile +\
			" -l " + ldFile +\
			" -g 0.001 " +\
			" -c 3 " +\
			" -f 1 -r 0.95\n";
	scriptFILEHandler = open(scriptfile+'.qsub', 'wb');
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="CAVIAR", logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('bsub < ' + scriptfile + '.qsub', shell=True);
	script = "";
