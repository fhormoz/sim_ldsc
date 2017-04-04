import os
import glob
import optparse
import numpy as np

def main(parser):
	(options, args) = parser.parse_args();
	annot_name = options.annot_name;
	print annot_name;
	ldsc_enr_est_list = [];
	true_enr_list_common = [];
	true_enr_list = [];
	true_per_hg2 = [];
	true_per_hg2_common = [];
	true_per_snp = [];
	true_per_snp_common = [];

	print "final_ldsc_result/" + annot_name + "/*.results";
	for simFiles in glob.glob("final_ldsc_result/" + annot_name + "/*.results"):
		enr_est = np.loadtxt(simFiles,skiprows=2, usecols=(4,));
		#if(float(enr_est)>0):
		ldsc_enr_est_list.append(float(enr_est));
	print ldsc_enr_est_list;
	rsID2beta = {};
	rsID2frq  = {};
	hapMapfrq = "/groups/price/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.1.frq"
	hapMapfrqHandler = open(hapMapfrq);
	hapMapfrqHandler.readline();
	for snpdataHapMap in hapMapfrqHandler:
		data = snpdataHapMap.split();
		snp = data[1];
		frq = data[4];
		rsID2frq[snp] = float(frq);	

	Nsim = hg2 = est_var_snp = num_snp_genome = num_snp_annot = 0;
	hg2_common = est_var_snp_common = num_snp_genome_common = num_snp_annot_common = 0;
	for simPath in glob.glob("simulated_data/sim*"):
		rsID2beta = {};
		hg2 = est_var_snp = num_snp_genome = num_snp_annot = 0;	
		hg2_common = est_var_snp_common = num_snp_genome_common = num_snp_annot_common = 0;
		simFileName = simPath.replace("simulated_data/", "");
		simFileGWAS = simPath + "/test_causalsnp.gwas";
		annotFile   = "annots/" + annot_name + "/" +  simFileName + "/" + annot_name + ".1.annot";
		print annotFile;
		if(not os.path.exists(annotFile)):
			continue;
		for gwasCausal in open(simFileGWAS):
			data = gwasCausal.split();
			rsID2beta[data[0]] = float(data[1]);
		annotFileHandler = open(annotFile);
		annotFileHandler.readline();
		for annotData in annotFileHandler:
			data = annotData.split();
			if (data[2] in rsID2beta):
				est_var_snp = est_var_snp + (float(data[5]) * rsID2beta[data[2]] * rsID2beta[data[2]]);
				hg2 = hg2 + (rsID2beta[data[2]] * rsID2beta[data[2]]);
			num_snp_annot = num_snp_annot + float(data[5]);
			num_snp_genome = num_snp_genome + 1;

			if (data[2] in rsID2beta and rsID2frq[data[2]] >= 0.05):
				est_var_snp_common = est_var_snp_common + (float(data[5]) * rsID2beta[data[2]] * rsID2beta[data[2]]);
				hg2_common = hg2_common + (rsID2beta[data[2]] * rsID2beta[data[2]]);	  
			if (rsID2frq[data[2]] >= 0.05):
				num_snp_annot_common = num_snp_annot_common + float(data[5]);
				num_snp_genome_common = num_snp_genome_common + 1;
		if (num_snp_annot_common == 0):
			continue;
		hg2 = 0.5;
		hg2_common = 0.5;
		true_enr_list.append((est_var_snp/hg2)/(float(num_snp_annot)/float(num_snp_genome)));
		true_enr_list_common.append((est_var_snp_common/hg2_common)/(float(num_snp_annot_common)/float(num_snp_genome_common)));
		true_per_hg2.append(est_var_snp/hg2);
		true_per_hg2_common.append(est_var_snp_common/hg2_common);

		true_per_snp.append(float(num_snp_annot)/float(num_snp_genome));
		true_per_snp_common.append(float(num_snp_annot_common)/float(num_snp_genome_common));
		
		Nsim = Nsim + 1;
	print ldsc_enr_est_list;
	print annot_name + "_all.out";
	outfileHandler = open(annot_name + "_all.out", "w");
	for ldsc_enr, true_enr in zip(ldsc_enr_est_list, true_enr_list):
		print ldsc_enr, true_enr;
		outfileHandler.write(str(ldsc_enr) + "\t" + str(true_enr)+"\n");
	outfileHandler.close();
	
	outfileHandler = open(annot_name + "_summary.out", "w");
	outfileHandler.write("LDSC Median ENR=\t" + str(np.median(ldsc_enr_est_list))+"\n");
	outfileHandler.write("LDSC Mean ENR=\t"   + str(np.mean(ldsc_enr_est_list))+"\n");
	outfileHandler.write("LDSC ENR(se)=\t"    + str(np.std(ldsc_enr_est_list)/np.sqrt(Nsim))+"\n");
	outfileHandler.write("TRUE ENR=\t"        + str(np.mean(true_enr_list))+"\n");
	outfileHandler.write("TRUE ENR(se)=\t"    + str(np.std(true_enr_list)/np.sqrt(Nsim))+"\n");
	outfileHandler.write("TRUE ENR COMMON=\t"        + str(np.mean(true_enr_list_common))+"\n");
        outfileHandler.write("TRUE ENR(se) COMMON=\t"    + str(np.std(true_enr_list_common)/np.sqrt(Nsim))+"\n");
	outfileHandler.write("TRUE %hg2=" + str(np.mean(true_per_hg2)) + "\n");
	outfileHandler.write("TRUE %hg2 COMMON =" + str(np.mean(true_per_hg2_common)) + "\n");
	outfileHandler.write("TRUE %SNP=" + str(np.mean(true_per_snp)) + "\n");
        outfileHandler.write("TRUE %SNP COMMON =" + str(np.mean(true_per_snp_common)) + "\n");
	outfileHandler.close();


if __name__ == "__main__":
        parser = optparse.OptionParser("usage: %prog [options] ")
        parser.add_option("-a", "--annot", dest="annot_name",
                default="1", type="string",
                help="specify the name of Annotation");
	main(parser);
