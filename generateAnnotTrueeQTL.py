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

	allciseQTLSNP = [];
	allciseQTLP   = [];
	trueCausalSNP = [];
	
	for data in open(currentPath + "/simulated_data/" + simFileName +"/test_causalsnp.gwas"):
		trueCausalSNP.append(data.split("\t")[0]);
	print trueCausalSNP;
	
	for eQTLFiles in glob.glob(currentPath + "/simulated_data/" + simFileName + "/*_out.qassoc"):
		if ("test_out.qassoc" in eQTLFiles):
			continue;
		else:
			data  = np.genfromtxt(eQTLFiles, dtype="str", usecols=(1,8), skip_header=1);
			allciseQTLSNP.extend(data[:,0]);
			allciseQTLP.extend(map(float, data[:,1]));

	allSigciseQTLSNP ={};
	allciseQTLPadj = mp.multipletests(allciseQTLP, method="fdr_bh", alpha=0.05);
	for snp, p, prej, padj in zip(allciseQTLSNP, allciseQTLP, allciseQTLPadj[0], allciseQTLPadj[1]):
		if(prej):
			allSigciseQTLSNP[snp] = p;


	if not os.path.exists("annots"):
                os.makedirs("annots")
	if not os.path.exists("annots/TrueQTL"):
                os.makedirs("annots/TrueQTL");
	if not os.path.exists("annots/TrueQTL/"+simFileName):
                os.makedirs("annots/TrueQTL/"+simFileName);

	bimFileHandler = open("/n/scratch2/fh80/UKBioBank_SimulatedData/1000G/1000G.EUR.QC.1.bim", "r");
	annotFileHandler = open(currentPath + "/" + "annots/TrueQTL/"+simFileName + "/TrueQTL.1.annot", "w");
	annotFileHandler.write("CHR\tBP\tSNP\tCM\tANNOT1\tANNOT2\n");
	for bimData in bimFileHandler:
		data = bimData.split();
		AnnotTrue = 0;
		if (data[1] in trueCausalSNP):
                        AnnotTrue = 1;
                else:
                        AnnotTrue = 0;
		annotFileHandler.write(data[0] + "\t" + data[3] + "\t" + data[1] + "\t" + data[2] + "\t1\t" + str(AnnotTrue) + "\n");
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
