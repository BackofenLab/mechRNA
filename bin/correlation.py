import sys, getopt, sets, os

from scipy import stats
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests, _ecdf as ecdf, fdrcorrection as fdrcorrection0, fdrcorrection_twostage
import numpy as np

#######################################################
#
# Main Code
#
#######################################################

def main(argv):

	# try:
	# 	opts, args = getopt.getopt(argv,"hi:r:m:c:p:x:",["ifile=","mfile="])#,"ofile="])
	# except getopt.GetoptError:
	# 	print 'infer.py -i <inputfile> -r <reference> -m <models> -c <cancer> -p <peak-cutoff>'
	# 	sys.exit(2)

	# for opt, arg in opts:
	# 	if opt == '-h':
	# 		print 'infer.py -i <inputfile> -r <reference> -m <models> -c <cancer> -p <peak-cutoff> -x <top-percent>'
	# 		sys.exit()
	# 	elif opt in ("-i", "--ifile"):
	# 		inputfile = arg
	# 	elif opt in ("-r", "--reference"):
	# 		ref = arg
	# 	elif opt in ("-m", "--models"):
	# 		modelsfile = arg
	# 	elif opt in ("-c", "--cancer"):
	# 		cancer_type = arg
	# 	elif opt in ("-p", "--peak-cutoff"):
	# 		p_cutoff = arg
	# 	elif opt in ("-x", "--topX"):
	# 		topx = arg

	expression_file = "/home/agawrons/data/expression/tcga/prostate.fpkm.norm.c03.test" 
	genes = list()
	data = list()
	pvalues = list()
	writer = open(expression_file + ".correl", 'w')

	with open(expression_file, 'r') as exp_in:
		samples = exp_in.readline()
		for line in exp_in:
			tokens = line[:-1].split()
			genes.append(tokens[0])
			data.append((np.array(tokens[1:])).astype(np.float))

	for x in xrange(0, len(data)):
   		for y in xrange(x+1, len(data)):
   			r_row, p_value = pearsonr(data[x], data[y])
			#if(p_value < 0.05):
			writer.write(genes[x] + "\t" + genes[y] + "\t" + str(r_row) + "\t" + str(p_value) + "\n")
			pvalues.append(p_value)

	writer.close()

	writer = open(expression_file + ".fdr", 'w')

	reject, pvalues_corr = fdrcorrection0(pvalues)

	for fdr in pvalues_corr:
		writer.write(str(fdr) + "\n")

	writer.close()
	

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))


#"Gene_Name,Accessible_Region,Target_Start,Target_End,lncRNA,lncRNA_Start,lncRNA_End,Free_Energy,LT_P-Value,
#Protein_lncRNA,Protein_Start,Protein_End,PL_P-Value,Protein_Target,Protein_Start,Protein_End,PT_P-Value"
#LT_Correl, LT_Correl_P-Value, PL_Correl, PL_Correl_P-Value, PT_Correl, PT_Correl_P-Value, Inferred_Mechanism, Combined_P-value
