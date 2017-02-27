import sys, getopt, sets

from sets import Set
from scipy import stats

# mech_table = dict()

# #No overlap

# mech_table["-11-11-10"] = "complex_formation_downregulation"
# mech_table["-1001-10"] = "localization_downregulation"
# mech_table["-1111-10"] = "nonsense"
# mech_table["-11-1000"] = "direct_downregulation"
# mech_table["-100000"] = "direct_downregulation"
# mech_table["-111000"] = "de-stabilization"
# mech_table["-11-1110"] = "nonsense"
# mech_table["-100110"] = "nonsense"
# mech_table["-111110"] = "nonsense"
# mech_table["01-11-10"] = "complex_formation_downregulation"
# mech_table["0001-10"] = "localization_downregulation"
# mech_table["0111-10"] = "nonsense"
# mech_table["01-1000"] = "protein_only_downregulation"
# mech_table["000000"] = "nonsense"
# mech_table["011000"] = "protein_only_upregulation"
# mech_table["01-1110"] = "nonsense"
# mech_table["000110"] = "localization_upregulation"
# mech_table["011110"] = "complex_formation_upregulation"
# mech_table["11-11-10"] = "nonsense"
# mech_table["1001-10"] = "nonsense"
# mech_table["1111-10"] = "nonsense"
# mech_table["11-1000"] = "nonsense"
# mech_table["100000"] = "direct_upregulation"
# mech_table["111000"] = "stabilization"
# mech_table["11-1110"] = "nonsense"
# mech_table["100110"] = "localization_upregulation"
# mech_table["111110"] = "complex_formation_upregulation"

# #Overlap

# mech_table["-11-11-11"] = "nonsense"
# mech_table["-1111-11"] = "competitive_downregulation"
# mech_table["-11-1001"] = "nonsense"
# mech_table["-111001"] = "competitive_downregulation"
# mech_table["-11-1111"] = "nonsense"
# mech_table["-111111"] = "nonsense"
# mech_table["01-11-11"] = "nonsense"
# mech_table["0111-11"] = "competitive_downregulation"
# mech_table["01-1001"] = "competitive_upregulation"
# mech_table["011001"] = "competitive_downregulation"
# mech_table["01-1111"] = "competitive_upregulation"
# mech_table["011111"] = "nonsense"
# mech_table["11-11-11"] = "nonsense"
# mech_table["1111-11"] = "nonsense"
# mech_table["11-1001"] = "competitive_upregulation"
# mech_table["111001"] = "nonsense"
# mech_table["11-1111"] = "competitive_upregulation"
# mech_table["111111"] = "nonsense"

# ENST00000484846	ENSG00000110047	EHD1	614	804	ENST00000441257	ENSG00000237036	ZEB1-AS1	263	462	92.7431	-172.9	1.86429	-80.1569	4.51460091532851e-10	1.84811344742948e-07	
# 0	1	1	NA	
# NA	NA	NA	NA	0	1	NA	
# 0	1	1	NA	
# NA	NA	NA	NA	0	1	NA	
# 0	1	1	NA	

def get_mechanism(candidate):

	#print candidate
	protein_target = False
	protein_lncRNA = False

	if((candidate[20] != "NA")):
		protein_target = True

	if((candidate[31] != "NA")):
		protein_lncRNA = True

	if(protein_target):
		if( ((int(candidate[3]) >= int(candidate[20])) and  (int(candidate[3]) <= int(candidate[21]))) or ((int(candidate[4]) >= int(candidate[20])) and  (int(candidate[4]) <= int(candidate[21]))) ):
			if(protein_lncRNA == True):
				if(float(candidate[16]) <= 0 and float(candidate[27]) > 0 and float(candidate[38]) < 0):
					return "competitive_downregulation"
				elif(float(candidate[16]) >= 0 and float(candidate[27]) < 0 and float(candidate[38]) > 0):
					return "competitive_upregulation"
				else:
					return "nonsense"
			else:
				if(float(candidate[16]) <= 0 and float(candidate[27]) > 0):
					return "competitive_downregulation"
				elif(float(candidate[16]) >= 0 and float(candidate[27]) < 0):
					return "competitive_upregulation"
				else:
					return "nonsense"
		else:
			if(protein_lncRNA == True):
				if(float(candidate[16]) <= 0 and float(candidate[27]) < 0 and float(candidate[38]) < 0):
					return "complex_formation_downregulation"
				elif(float(candidate[16]) >= 0 and float(candidate[27]) < 0 and float(candidate[38]) > 0):
					return "complex_formation_upregulation"
				else:
					return "nonsense"
			else:
				if(float(candidate[16]) < 0):
					return "de-stabilization"
				elif(float(candidate[16]) > 0):
					return "stabilization"
				else:
					return "nonsense"
	else:
		if(protein_lncRNA == True):
			if(float(candidate[16]) <= 0 and float(candidate[38]) < 0):
				return "localization_downregulation"
			elif(float(candidate[16]) >= 0  and float(candidate[38]) > 0):
				return "localization_upregulation"
			else:
				return "nonsense"
		else:
			if(float(candidate[16]) < 0):
				return "direct_downregulation"
			elif(float(candidate[16]) > 0):
				return "direct_upregulation"
			else:
				return "nonsense"


def main(argv):

	try:
		opts, args = getopt.getopt(argv,"hi:r:m:c:p:",["ifile=","mfile="])#,"ofile="])
	except getopt.GetoptError:
		print 'infer.py -i <inputfile> -r <reference> -m <models> -c <cancer> -p <peak-cutoff>'
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print 'infer.py -i <inputfile> -r <reference> -m <models> -c <cancer> -p <peak-cutoff>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-r", "--reference"):
			ref = arg
		elif opt in ("-m", "--models"):
			modelsfile = arg
		elif opt in ("-c", "--cancer"):
			cancer_type = arg
		elif opt in ("-p", "--peak-cutoff"):
			p_cutoff = arg
	
	#correlationfile = "/home/agawrons/data/net.edges.sig.genes"
	correlationfile = "./data/{0}.net.edges.sig.genes".format(cancer_type)
	proteinfile = "./data/peaks/graphprot.peaks.{0}.sig".format(ref) 
	#correlationfile = "../MechRNA-v0.1/data/test.net"
	#proteinfile = "./test/test.peaks" 

	intarna_data = list()
	models = list()
	lncRNA_peak_keys = dict()
	graphprot_data = dict()
	correl_lncRNA_data = dict()
	correl_proteins_data = dict()

	with open(modelsfile, 'r') as model_in:
		for model in model_in:
			model = model[:-1]
			models.append(model)
			correl_proteins_data[model] = dict()
			graphprot_data[model] = dict()
		model = "NA;NA"
		models.append(model)
		correl_proteins_data[model] = dict()
		graphprot_data[model] = dict()
	
	with open(inputfile, 'r') as lncRNA_in:
		line = lncRNA_in.readline()
		tokens = line.split();
		lncRNA = tokens[3].split(";")
		lncRNA[0] = lncRNA[0][1:] 

	#print "Loading correlation file..."		
	with open(correlationfile, 'r') as correl_in:
		for line in correl_in:
			tokens = line.split()
			if(tokens[0] == lncRNA[1]):
				correl_lncRNA_data[tokens[1]] = tokens[2:6]
			if(tokens[1] == lncRNA[1]):
				correl_lncRNA_data[tokens[0]] = tokens[2:6]
			for model in models:
				protein_id = model.split(";")[0]
				protein_gene = model.split(";")[1]
				if(tokens[0] == protein_id):
					correl_proteins_data[model][tokens[1]] = tokens[2:6]
				elif(tokens[1] == protein_id):
					correl_proteins_data[model][tokens[0]] = tokens[2:6]
				elif(protein_id == "NA"):
					correl_proteins_data[model][tokens[0]] = ["0","1","1","NA"]
					correl_proteins_data[model][tokens[1]] = ["0","1","1","NA"]
	
	#print "Loading protein file..."
	with open(proteinfile, 'r') as graphprot_in:
		for line in graphprot_in:
			tokens = line.split()
			ids = tokens[0].split(";")
			trans = ids[0]
			gene = ids[1]
			if (trans == lncRNA[0]):
				for protein in tokens[3].split(","):
					if (float(tokens[5]) < float(p_cutoff)):
						if protein not in lncRNA_peak_keys:
							lncRNA_peak_keys[protein] = list()
						lncRNA_peak_keys[protein].append(tokens[0])
			for model in models:
				if(model == "NA;NA"):
					peak = ["NA","NA","NA","NA","0","1","NA","0","1","1","NA"]
					if tokens[0] not in graphprot_data[model]:
						graphprot_data[model][tokens[0]] = list()
						graphprot_data[model][tokens[0]].append(peak)
				else:
					for protein in tokens[3].split(","):
						protein_id = protein.split(";")[0]
						protein_gene = protein.split(";")[1]
						if(protein == model):
							peak = tokens[1:]
							peak[2] = protein_id
							peak.insert(3, protein_gene)
							if gene in correl_proteins_data[model]:
								peak.extend(correl_proteins_data[model][gene])
							else:
								peak.extend(["0","1","1","NA"])
							if tokens[0] not in graphprot_data[model]:
								graphprot_data[model][tokens[0]] = list()
							graphprot_data[model][tokens[0]].append(peak)
					


	#print "Reading RNA-RNA interactions file..."
	with open(inputfile, 'r') as intarna_in:

		cur_transcript = None
		RR_interactions = list()
		candidates = list()

		for line in intarna_in:
			tokens = line.split()
			ids = tokens[0][1:]
			RR_interaction = ids.split(';')
			RR_interaction.extend(tokens[1:3])
			RR_interaction.extend(tokens[3][1:].split(';'))
			RR_interaction.extend(tokens[4:])
			transcript = RR_interaction[0]
			if(cur_transcript == None): cur_transcript = transcript

			if(transcript != cur_transcript):
				gene = RR_interactions[0][1]
				target_ids = ';'.join(RR_interactions[0][0:3])
				lncRNA_ids = ';'.join(RR_interactions[0][5:8])
				top_candidate = None
				ties = 0
				for RR_inter in RR_interactions:
					for protein_target in models:
						pt_correl = correl_proteins_data[protein_target].get(gene)
						pt_peaks = graphprot_data[protein_target].get(target_ids)
						if(protein_target == "NA;NA"):
							pt_correl = ["0","1","1","NA"]
							pt_peaks = list()
							pt_peaks.append(["NA","NA","NA","NA","0","1","NA","0","1","1","NA"])
						if(pt_correl == None or pt_peaks == None): continue
						for peak_pt in pt_peaks:
							for protein_lncRNA in models:
								pl_correl = correl_proteins_data[protein_lncRNA].get(gene)
								pl_peaks = graphprot_data[protein_lncRNA].get(target_ids)
								if(protein_lncRNA == "NA;NA"): 
									pl_correl = ["0","1","1","NA"]
									pl_peaks = list()
									pl_peaks.append(["NA","NA","NA","NA","0","1","NA","0","1","1","NA"])
								if(pl_correl == None or pl_peaks == None or (protein_target == protein_lncRNA and protein_target != "NA;NA")): continue
								for peak_pl in pl_peaks: 
									if(float(peak_pl[5]) > float(p_cutoff) and peak_pl[0] != "NA"): continue
									candidate = list()
									candidate.extend(RR_inter)
									candidate.extend(correl_lncRNA_data.get(gene, ["0","1","1","NA"]))#[gene])
									candidate.extend(peak_pt)
									candidate.extend(peak_pl)

									mechanism = get_mechanism(candidate)
									
									if(mechanism != "nonsense"):

										candidate.extend([mechanism])
										pvalues = list()

										for index in [15,17,25,28,36,39]:
											if(float(candidate[index]) != 1): 
												pvalues.append(float(candidate[index]))

										if(len(pvalues) > 1): 
											xsq, p = stats.combine_pvalues(pvalues, method='fisher', weights=None)
											candidate.extend([p])
										else:
											candidate.extend([pvalues[0]])

										if(top_candidate == None or float(candidate[43]) < float(top_candidate[43])):
											top_candidate = candidate
											ties = 0
										elif(top_candidate != None and float(candidate[43]) == float(top_candidate[43])):
											ties += 1
									#print candidate
				if(top_candidate != None):
					print '\t'.join(map(str, top_candidate))
				#print ties

				cur_transcript = transcript
				RR_interactions = list()
				candidates = list()

			RR_interactions.append(RR_interaction)

	

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))


#"Gene_Name,Accessible_Region,Target_Start,Target_End,lncRNA,lncRNA_Start,lncRNA_End,Free_Energy,LT_P-Value,
#Protein_lncRNA,Protein_Start,Protein_End,PL_P-Value,Protein_Target,Protein_Start,Protein_End,PT_P-Value"
#LT_Correl, LT_Correl_P-Value, PL_Correl, PL_Correl_P-Value, PT_Correl, PT_Correl_P-Value, Inferred_Mechanism, Combined_P-value
