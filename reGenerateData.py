import os
import glob
import optparse
import subprocess
import numpy as np
from bisect import bisect

TEMPLATE_SERIAL = """
#####################################
#!/bin/bash
#BSUB -n 1                #each  job run on 1 core
#BSUB -W 03:00            #job run 5 hour
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
class Bed(object):
	start = 0;
	end   = 0;
	rsIDStart = "";
	rsIDEnd   = "";

	def __init__(self, start, end, rsIDStart, rsIDEnd):
		self.start = start;
		self.end   = end;
		self.rsIDStart = rsIDStart;
		self.rsIDEnd   = rsIDEnd;
	def __str__(self):
		return str(self.start) + "_" + str(self.end);
	def __repr__(self):
		return str(self.start) + "_" + str(self.end);

def resubmitScript2Orchestra(scriptfile,jobId):
	scriptfile = "log/" + str(jobId) + '.qsub';
        logfile = "log/" + str(jobId) + '.log';
        errfile = "log/" + str(jobId) + '.err';
        subprocess.call('bsub < ' + scriptfile, shell=True);
	
def submitScript2Orchestra(script, jobId):
	scriptfile = "log/" + str(jobId) + '.qsub';
	scriptFILEHandler = open(scriptfile, 'wb');
	logfile = "log/" + str(jobId) + '.log';
	errfile = "log/" + str(jobId) + '.err';
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="Re-RUN", logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('bsub < ' + scriptfile, shell=True);

def main(parser):
	gwashg2 = 0.5;
	eqtlhg2 = 0.1;
	(options, args) = parser.parse_args();
	simID = options.simID;
	outFolder = options.outFolder;
	geneBedFile = "/groups/price/farhad/Data/Genes/Wright_chr1_genePos.bed";
	snpBedFileFam = "/n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000/UKBioBank_chr_40000ind_GWASQC.bim";
	
	currentPath = os.getcwd() + "/"+ outFolder;
	os.chdir(currentPath);
	if(not os.path.exists("sim" + simID)):
		os.makedirs("sim" + simID);
		os.makedirs(currentPath + "/sim" + simID + "/log");
	os.chdir(currentPath + "/sim" + simID);	
	print currentPath + "/sim" + simID + "/log";
	
	PLINK_PATH = "/opt/plink2-1.90b3/bin/plink" 
	GCTA_PATH  = "/groups/price/farhad/Software/gcta/binary/gcta64";
	eqtlPlinkFile = "/n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000/UKBioBank_chr_500ind_eQTLQC"

	##RUN GCTA to generate the GWAS phenotype
	COMMAND="";
	command_index = 10000000;
	for logFile in glob.glob("log/*.log"):
		failedCOMMAND = "";
		erroFlag = False;
		for data in open(logFile):
			if ("Error" in data):
				erroFlag = True;
		if(erroFlag):
			inputLogFile = open(logFile);
			index = 0;
			print logFile;
			inputLogFile.readline();
			inputLogFile.readline();
			eqtlFileName = logFile.replace("_out.log", "_causalsnp.eqtl");
			phenFileName = "";
			for index in range(0,8):
				dataLine = inputLogFile.readline().rstrip();
				if ("_phen" in dataLine):
					indexFind = dataLine.find("test");
					phenFileName = dataLine[indexFind:].replace(".phen","");
				failedCOMMAND += dataLine;
				
			COMMAND += GCTA_PATH  + " --simu-qt --bfile " + eqtlPlinkFile +\
                         " --out "+ phenFileName +\
                         " --simu-causal-loci " +  eqtlFileName + "\n";

			COMMAND += PLINK_PATH + failedCOMMAND + "\n";
	
	submitScript2Orchestra(COMMAND, command_index);

	for qsubFile in glob.glob("log/*.qsub"):
		if(not os.path.exists(qsubFile.replace(".qsub",".err"))):
			jobId = qsubFile.replace(".qsub","").replace("log/","");
			print jobId;
			resubmitScript2Orchestra(qsubFile, jobId);	

if __name__ == "__main__":
        parser = optparse.OptionParser("usage: %prog [options] ")
        parser.add_option("-s", "--sim", dest="simID",
                default="1", type="string",
                help="specify the simulation ID");
	parser.add_option("-o", "--out", dest="outFolder",
                default="simulation", type="string",
                help="specify the folder to keep the simulation files");
	main(parser);
