import os
import glob
import optparse
import subprocess
import numpy as np
import statsmodels.sandbox.stats.multicomp as mp
from scipy.stats import norm

def main(parser):
        (options, args) = parser.parse_args();
        simPath = options.simPath;
	simFileName = options.simFileName

	currentPath = os.getcwd();
	simPath = currentPath + "/simulated_data";
	
	outputPath = simPath + "/SummaryStatistics";
	if not os.path.exists(outputPath):
		os.makedirs(outputPath)

	allCausalSetCAVIAR = {};
	for eQTLFiles in glob.glob(currentPath + "/simulated_data/" + simFileName + "/caviar_output/*_post"):
		for data in open(eQTLFiles):
			splitData = data.split();
			snp  = splitData[0];
			post = float(splitData[2]);
			if (snp in allCausalSetCAVIAR):
				tmpPost = allCausalSetCAVIAR[snp];
				maxPPC  = np.maximum(tmpPost, post);
				allCausalSetCAVIAR[snp] = maxPPC;
			else:
				allCausalSetCAVIAR[snp] = post;	

	if not os.path.exists("annots"):
                os.makedirs("annots")
	if not os.path.exists("annots/MaxPPC_CAVIAR"):
                os.makedirs("annots/MaxPPC_CAVIAR");
	if not os.path.exists("annots/MaxPPC_CAVIAR/"+simFileName):
                os.makedirs("annots/MaxPPC_CAVIAR/"+simFileName);

	bimFileHandler = open("/n/scratch2/fh80/UKBioBank_SimulatedData/UKBiobank/500_40000_HapMapCommon/UKBioBank_chr_500ind_eQTL.CM.bim", "r");
	annotFileHandler = open(currentPath + "/" + "annots/MaxPPC_CAVIAR/"+simFileName + "/MaxPPC_CAVIAR.1.annot", "w");
	annotFileHandler.write("CHR\tBP\tSNP\tCM\tANNOT1\tANNOT2\n");
	for bimData in bimFileHandler:
		data = bimData.split();
		Annot = 0;
		if (data[1] in allCausalSetCAVIAR):
			Annot = allCausalSetCAVIAR[data[1]];
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
