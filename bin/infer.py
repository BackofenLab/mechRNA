import sys, getopt, sets

from sets import Set

def parse_accessible_region(region):

	ids = list()
	region = region[1:]

	for token in region.split(';'):
		trans = token.split(',',1)
		ids.append(trans[0])

	return ids


def main(argv):

	try:
		opts, args = getopt.getopt(argv,"hi:r:m:",["ifile=","mfile="])#,"ofile="])
	except getopt.GetoptError:
		print 'infer.py -i <inputfile> -r <reference> -m <models>'
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print 'infer.py -i <inputfile> -r <reference> -m <models>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-r", "--reference"):
			ref = arg
		elif opt in ("-m", "--models"):
			modelsfile = arg

	correlationfile = "../data/net.edges.sig.genes"
	proteinfile = "../data/peaks/graphprot.peaks.{0}.sig".format(ref) 
	#correlationfile = "../data/test.net"
	#proteinfile = "../data/test.peaks.pvalues.genes" 

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

	with open(inputfile, 'r') as lncRNA_in:
		line = lncRNA_in.readline()
		tokens = line.split();
		lncRNA = tokens[3].split(";")
		lncRNA[0] = lncRNA[0][1:] 
	
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
				if(tokens[1] == protein_id):
					correl_proteins_data[model][tokens[0]] = tokens[2:6]

	with open(proteinfile, 'r') as graphprot_in:
		for line in graphprot_in:
			tokens = line.split()
			ids = tokens[0].split(";")
			trans = ids[0]
			gene = ids[1]
			if (trans == lncRNA[0]):
				for protein in tokens[3].split(","):
					if protein not in lncRNA_peak_keys:
						lncRNA_peak_keys[protein] = list()
					lncRNA_peak_keys[protein].append(tokens[0])
			for model in models:
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

	with open(inputfile, 'r') as intarna_in:
		for line in intarna_in:
			tokens = line.split()
			ids = tokens[0][1:]
			interaction = ids.split(';')
			interaction.extend(tokens[1:3])
			interaction.extend(tokens[3][1:].split(';'))
			interaction.extend(tokens[4:])
			if interaction[1] in correl_lncRNA_data:
				interaction.extend(correl_lncRNA_data[interaction[1]])
			else:
				interaction.extend(["0","1","1","NA"])
			for model in models:
				annotated_interaction = list(interaction)
				if ids in graphprot_data[model]:
					for peak in graphprot_data[model][ids]:
						if((int(annotated_interaction[3]) >= int(peak[0]) and int(annotated_interaction[3]) <= int(peak[1])) or (int(annotated_interaction[4]) >= int(peak[0]) and int(annotated_interaction[4]) <= int(peak[1]))):# or (peak[5] != 0)):
							annotated_interaction_peak = list(annotated_interaction)
							annotated_interaction_peak.extend(peak)
							intarna_data.append(annotated_interaction_peak)
				else:
					if model in lncRNA_peak_keys:
						protein_id = model.split(";")[0]
						protein_gene = model.split(";")[1]
						if annotated_interaction[1] in correl_proteins_data[model]:
							annotated_interaction.extend(["NA","NA",protein_id,protein_gene,"0","1"])
							annotated_interaction.extend(correl_proteins_data[model][annotated_interaction[1]])
							intarna_data.append(annotated_interaction)
					#	else:
					#		annotated_interaction.extend(["NA","NA","NA","0","1","0","1","1","NA"])
					#else:
						#annotated_interaction.extend(["NA","NA","NA","0","1","0","1","1","NA"]) 
					

		# print intarna_data[0]
		# print intarna_data[0][6]
		# print intarna_data[0][7]
		# print intarna_data[0][23]
		# print intarna_data[0][24]
		# print graphprot_data["ENSG00000066044"]["ENST00000555375,1,102,18,80477"]

	########################################################################################################
	# Inference 6-7 15-18 19-24 25-28
	# IntaRNA 0-15 correl 16-19
	# GraphProt 20-25 correl 26-29
	#print intarna_data[0]

	for mech in intarna_data:
		protein_target = False
		protein_lncRNA = False
		overlap = False
		protein = ""
		protein_symbol = ""
		
		#Determine features of interaction
		if((mech[20] != "NA") and (float(mech[27]) != 1)):
			protein_target = True
			protein = mech[22]
			protein_symbol = mech[23]
		if(any(lncRNA_peak_keys)):
			protein_lncRNA = True
		if(protein_target):
			if( ((int(mech[3]) >= int(mech[20])) and  (int(mech[3]) <= int(mech[21]))) or ((int(mech[4]) >= int(mech[20])) and  (int(mech[4]) <= int(mech[21]))) ):
				overlap = True

		#Mechanism inference
		if((protein_target == False) and (protein_lncRNA == False) and (overlap == False)):
			mech.append("direct")
		elif((protein_target == True) and (protein_lncRNA == False or (protein not in lncRNA_peak_keys)) and (overlap == True)):
			if(float(mech[27]) > 0):
				mech.append("competitive_de-stabilizing")
			elif(float(mech[27]) < 0):
				mech.append("competitive_stabilizing")
			else:
				mech.append("competitive")
		elif((protein_target == True) and (protein_lncRNA == True) and (overlap == False)):
			if(protein in lncRNA_peak_keys):
				mech.append("localization_single")
			else:
				mech.append("localization_multiple")
		elif((protein_target == False) and (protein_lncRNA == True) and (overlap == False)):
			mech.append("localization")
		elif((protein_target == True) and (protein_lncRNA == False) and (overlap == False)):
			if(float(mech[27]) > 0):
				mech.append("stabilizing")
			elif(float(mech[27]) < 0):
				mech.append("de-stabilizing")
			else:
				mech.append("direct")
		elif((protein_target == True) and (protein_lncRNA == True) and (overlap == True)):
			if(protein_symbol == "AGO1"):
				mech.append("miRNA-related")
			elif(protein_symbol == "STAU"):
				mech.append("STAU-mediated_decay")
			else:
				mech.append("Unknown")
		else:
			mech.append("Unknown")
		print '\t'.join(mech)

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))


#"Gene_Name,Accessible_Region,Target_Start,Target_End,lncRNA,lncRNA_Start,lncRNA_End,Free_Energy,LT_P-Value,
#Protein_lncRNA,Protein_Start,Protein_End,PL_P-Value,Protein_Target,Protein_Start,Protein_End,PT_P-Value"
#LT_Correl, LT_Correl_P-Value, PL_Correl, PL_Correl_P-Value, PT_Correl, PT_Correl_P-Value, Inferred_Mechanism, Combined_P-value