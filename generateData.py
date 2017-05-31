import os
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

def submitScript2Orchestra(script, jobId):
	scriptfile = "log/" + str(jobId) + '.qsub';
	scriptFILEHandler = open(scriptfile, 'wb');
	logfile = "log/" + str(jobId) + '.log';
	errfile = "log/" + str(jobId) + '.err';
	scriptFILEHandler.write(TEMPLATE_SERIAL.format(script=script, name="LDSC", logfile=logfile, errfile=errfile, slots=1))
	scriptFILEHandler.close();
	subprocess.call('bsub < ' + scriptfile, shell=True);

def main(parser):
	gwashg2 = 0.5;
	eqtlhg2 = 0.1;
	(options, args) = parser.parse_args();
	simID = options.simID;
	outFolder = options.outFolder;
	geneBedFile = "/groups/price/farhad/Data/Genes/Wright_chr1_genePos.bed";
	snpBedFileFam = "/n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000_HapMapCommon/UKBioBank_chr_40000ind_GWAS.CM.bim";
	snpBedFamData = np.genfromtxt(snpBedFileFam, dtype='str');
	snpPos  = map(int, snpBedFamData[:,3]);
	snprsID = snpBedFamData[:,1]; 

	##Generate the Causal variants for GWAS and eQTL
	causalIndexList  = [];
	causalGeneList   = [];
	allGeneListStart = [];
	for data in open(geneBedFile):
		fileds = data.split();
		if (int(fileds[3]) in allGeneListStart):
			continue;
		indexStart = bisect(snpPos, int(fileds[3]) );
		indexEnd   = bisect(snpPos, int(fileds[4]) );
		if(indexStart < indexEnd):
			#causalIndex = np.random.randint(indexStart, indexEnd);
			causalIndex = indexStart + 1;
			if(not (causalIndex in causalIndexList) ):
				causalIndexList.append(causalIndex);
				causalGeneList.append(Bed(fileds[3], fileds[4], snprsID[indexStart], snprsID[indexEnd]));
				allGeneListStart.append(int(fileds[3]));

	##MAKE THE SIMULATION DIRECTORY	
	currentPath = os.getcwd() + "/"+ outFolder;
	os.chdir(currentPath);
	if(not os.path.exists("sim" + simID)):
		os.makedirs("sim" + simID);
		os.makedirs(currentPath + "/sim" + simID + "/log");
	os.chdir(currentPath + "/sim" + simID);	
	print currentPath + "/sim" + simID;
	ouputFileName = "test"
	gwasPlinkFile = "/n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000_HapMapCommon/UKBioBank_chr_40000ind_GWAS.CM";
	eqtlPlinkFile = "/n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000_HapMapCommon/UKBioBank_chr_500ind_eQTL.CM"
	##Generate the effect size of GWAS
	gwasCausalFile = open(ouputFileName + "_causalsnp.gwas", "w");
	for index in causalIndexList:
		gwasCausalFile.write(snpBedFamData[index,1] + "\t" + str(np.random.normal(0,scale=np.sqrt(gwashg2/len(causalIndexList)))) + "\n" );
	gwasCausalFile.close();

	GCTA_PATH  = "/groups/price/farhad/Software/gcta/binary/gcta64";
	PLINK_PATH = "/opt/plink2-1.90b3/bin/plink" 
	##RUN GCTA to generate the GWAS phenotype
	COMMAND = GCTA_PATH + " --simu-qt --bfile " + gwasPlinkFile + " --out " + ouputFileName + "_phen" + " --simu-causal-loci " +  ouputFileName + "_causalsnp.gwas\n";
	COMMAND += PLINK_PATH + " --noweb --bfile " + gwasPlinkFile + " --out " + ouputFileName + "_out"  + " --pheno " + ouputFileName + "_phen.phen --allow-no-sex --assoc\n";
	submitScript2Orchestra(COMMAND, 1);
	##Generate the effect size of eQTL
	COMMAND="";
	command_index = 0;
	for index, bedData in zip(causalIndexList, causalGeneList):
		print index, bedData, command_index;
		eqtlFileName = ouputFileName + "_" + str(bedData) + "_causalsnp.eqtl";
		eqtlCausalFile = open(eqtlFileName, "w");
		eqtlCausalFile.write(snpBedFamData[index,1] + "\t" + str(np.random.normal(0,scale=np.sqrt(eqtlhg2))) + "\n" );
		eqtlCausalFile.close();
		##RUN GCTA to generate the eQTL phenotype	
		COMMAND += GCTA_PATH  + " --simu-qt --bfile " + eqtlPlinkFile +\
			 " --out "+ ouputFileName + "_" + str(bedData) + "_phen" +\
			 " --simu-causal-loci " +  eqtlFileName + "\n";
		COMMAND += PLINK_PATH + " --noweb --bfile " + eqtlPlinkFile +\
			" --from " + str(bedData.rsIDStart) + " --to " + str(bedData.rsIDEnd) +\
			" --out " + ouputFileName + "_" + str(bedData) +  "_out"  +\
			" --pheno " + ouputFileName + "_" + str(bedData) + "_phen.phen --allow-no-sex --assoc\n";
		command_index = command_index + 1;
		if (command_index % 50 == 0):
			submitScript2Orchestra(COMMAND, command_index);
			COMMAND = "";
	if(COMMAND != ""):
		submitScript2Orchestra(COMMAND, command_index);

if __name__ == "__main__":
        parser = optparse.OptionParser("usage: %prog [options] ")
        parser.add_option("-s", "--sim", dest="simID",
                default="1", type="string",
                help="specify the simulation ID");
	parser.add_option("-o", "--out", dest="outFolder",
                default="simulation", type="string",
                help="specify the folder to keep the simulation files");
	main(parser);
