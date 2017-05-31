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
	
	TopciseQTL = {};
	for eQTLFiles in glob.glob(currentPath + "/simulated_data/" + simFileName + "/*_out.qassoc"):
		if ("test_out.qassoc" in eQTLFiles):
			continue;
		else:
			data  = np.genfromtxt(eQTLFiles, dtype="str", usecols=(1,8), skip_header=1, missing_values='NA', filling_values='1');
			index = np.argmin(map(float, data[:,1]));
			TopciseQTL[data[index,0]] = 1;

	if not os.path.exists("annots"):
                os.makedirs("annots")
	if not os.path.exists("annots/TopciseQTL"):
                os.makedirs("annots/TopciseQTL");
	if not os.path.exists("annots/TopciseQTL/"+simFileName):
                os.makedirs("annots/TopciseQTL/"+simFileName);

	bimFileHandler = open("/n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000_HapMapCommon/UKBioBank_chr_500ind_eQTL.CM.bim", "r");
	annotFileHandler = open(currentPath + "/" + "annots/TopciseQTL/"+simFileName + "/TopciseQTL.1.annot", "w");
	annotFileHandler.write("CHR\tBP\tSNP\tCM\tANNOT1\tANNOT2\n");
	for bimData in bimFileHandler:
		data = bimData.split();
		Annot = 0;
		if (data[1] in TopciseQTL):
			Annot = 1;
		else:
			Annot = 0;	
		annotFileHandler.write(data[0] + "\t" + data[3] + "\t" + data[1] + "\t" + data[2] + "\t1\t" +  str(Annot) + "\n");
	annotFileHandler.close();
	bimFileHandler.close();

if __name__ == "__main__":
        parser = optparse.OptionParser("usage: %prog [options] ")
        parser.add_option("-f", "--simFILE", dest="simFileName",
                default="1", type="string",
                help="specify the simulation FILE");
	parser.add_option("-p", "--simPATH", dest="simPath",
                default="1", type="string",
                help="specify the simulation Path")

	main(parser);
