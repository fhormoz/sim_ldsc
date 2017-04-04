import os
import glob
import optparse
import subprocess
import numpy as np
from scipy.stats import norm

def main(parser):
        (options, args) = parser.parse_args();
        simPath = options.simPath;
	simFileName = options.simFileName

	currentPath = os.getcwd();
	simPath = currentPath + "/simulated_data";
	
	PLINK_PATH = "/opt/plink2-1.90b3/bin/plink"
	
	for eQTLFiles in glob.glob(currentPath + "/simulated_data/" + simFileName + "/*_out.qassoc"):
		if ("test_out.qassoc" in eQTLFiles):
			continue;
		else:
			if not os.path.exists("simulated_data/"+simFileName + "/caviar_input/"):
                		os.makedirs("simulated_data/"+simFileName + "/caviar_input/")
			data  = np.genfromtxt(eQTLFiles, dtype="str", usecols=(1,4,8), skip_header=1);
			eQTLFileName = eQTLFiles.replace(currentPath + "/simulated_data/" + simFileName + "/", "");
			eQTLFileNameHandler = open("simulated_data/"+simFileName + "/caviar_input/"+eQTLFileName+".Z", "w")
			eQTLFileSNPHandler  = open("simulated_data/"+simFileName + "/caviar_input/"+eQTLFileName+"_snp", "w")

			Ptmp = map(float, data[:,2]);
			sigCutoff = 1;
			maxSNP = 50;
			if (len(Ptmp) > maxSNP):
				Ptmp = np.sort(Ptmp);
				sigCutoff = np.abs(Ptmp[maxSNP]); 
			for snp,beta,P in data:
				Z = np.abs(norm.ppf(float(P)/2));
				if(float(beta) < 0):
					Z = -Z;
				if(float(P) < sigCutoff):
					eQTLFileNameHandler.write( snp + "\t" + str(Z) + "\n" );
					eQTLFileSNPHandler.write( snp + "\n");
			eQTLFileNameHandler.close();
			eQTLFileSNPHandler.close();
			
			LDFile  = "simulated_data/"+simFileName + "/caviar_input/"+eQTLFileName;
			snpFile = "simulated_data/"+simFileName + "/caviar_input/"+eQTLFileName + "_snp";
			COMMAND = PLINK_PATH + " --bfile /n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000_HapMapCommon/UKBioBank_chr_500ind_eQTL.CM --r --matrix --out " + LDFile  + "  --extract " + snpFile; 
			subprocess.call(COMMAND, shell=True);
			os.remove(snpFile);
			os.remove("simulated_data/"+simFileName + "/caviar_input/"+eQTLFileName+".log");		
			os.remove("simulated_data/"+simFileName + "/caviar_input/"+eQTLFileName+".nosex");

if __name__ == "__main__":
        parser = optparse.OptionParser("usage: %prog [options] ")
        parser.add_option("-f", "--simFILE", dest="simFileName",
                default="1", type="string",
                help="specify the simulation FILE");
	parser.add_option("-p", "--simPATH", dest="simPath",
                default="1", type="string",
                help="specify the simulation Path")

	main(parser);
