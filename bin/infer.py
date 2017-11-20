import sys, getopt, sets, os
import statsmodels.api as sm
import subprocess

from sets import Set
from scipy import stats
import numpy as np
from statsmodels.stats.multitest import multipletests, _ecdf as ecdf, fdrcorrection as fdrcorrection0, fdrcorrection_twostage
from heapq import heapify, heappop, heappush
from itertools import islice, cycle
from tempfile import gettempdir

#from subprocess import call

indels = dict()
upregulators = dict()
downregulators = dict()
complexes = dict()
mechs = set()
all_genes = dict()
all_trans = dict()
LOCAL_DIST = 1000
OVERLAP_DIST = 50

#rna_inter_members = ["t_trans_id","t_gene_id","t_gene_symbol","l_trans_id","l_gene_id","l_gene_symbol","t_start","t_end","l_start","l_end","E_init","E_loops","E_dangleL","E_dangleR","E_endL","E_endR","ED1","ED2","E","pvalue","fdr","mutations","anno1","anno2"]
rna_inter_members = ["t_trans_id","t_gene_id","t_gene_symbol","l_trans_id","l_gene_id","l_gene_symbol","t_start","t_end","l_start","l_end","E","pvalue","fdr","mutations","anno1","anno2"]
#protein_inter_members = ["start", "end", "gene_id", "gene_symbol", "peak", "pvalue", "pid", "mutations"]
protein_inter_members = ["start", "end", "gene_id", "gene_symbol", "pvalue", "mutations"]
#correlation_members = ["pcor", "pvalue", "dir_pvalue", "dir"]
correlation_members = ["pcor", "pvalue"]

#######################################################
#
# Data Structures
#
#######################################################

class Candidate:
	def __init__(self, rna_inter):

		self.rna = rna_inter

		self.correl = Correlation([0,1,1,"NA"])

		self.t_peak = Protein_Interaction(["NA","NA","NA","NA","0","1","NA"], Correlation([0,1,1,"NA"]))
		self.l_peak = Protein_Interaction(["NA","NA","NA","NA","0","1","NA"], Correlation([0,1,1,"NA"]))

		self.mechanism = "NA"
		self.joint_pvalue = 1

	def get_pvalues(self):
		pvalues = list()
		pvalues.append(float(self.rna.pvalue))
		if(float(self.correl.pvalue) != 1): pvalues.append(float(self.correl.pvalue))
		if(float(self.t_peak.pvalue) != 1): pvalues.append(float(self.t_peak.pvalue))
		if(float(self.t_peak.correl.pvalue) != 1): pvalues.append(float(self.t_peak.correl.pvalue))
		if(float(self.l_peak.pvalue) != 1): pvalues.append(float(self.l_peak.pvalue))
		if(float(self.l_peak.correl.pvalue) != 1): pvalues.append(float(self.l_peak.correl.pvalue))
		for i in xrange(len(pvalues)):
			if(pvalues[i] == 0):
				pvalues[i] = float(1e-22)
		return pvalues

	def toString(self):
		return "{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(self.rna.toString(), self.correl.toString(), self.t_peak.toString(), self.l_peak.toString(), self.mechanism, self.joint_pvalue) 

class RNA_Interaction:
	def __init__(self, data):

		tokens = data.split(';')#'\t')

		self.t_trans_id = tokens[0][0:15]
		self.t_gene_id = tokens[0][16:31]
		self.t_gene_symbol = tokens[0][32:]
		self.l_trans_id = tokens[1][0:15]
		self.l_gene_id = tokens[1][16:31]
		self.l_gene_symbol = tokens[1][32:]
		self.t_start = tokens[2]
		self.t_end = tokens[3]
		self.l_start = tokens[4]
		self.l_end = tokens[5]
		#self.E_init = tokens[6]
		#self.E_loops = tokens[7]
		#self.E_dangleL = tokens[8]
		#self.E_dangleR = tokens[9]
		#self.E_endL = tokens[10]
		#self.E_endR = tokens[11]
		#self.ED1 = tokens[12]
		#self.ED2 = tokens[13]
		self.E = tokens[14]
		self.pvalue = 0 #tokens[15]
		self.fdr = 0 #tokens[16]
		self.mutations = 0
		self.anno1 = ""
		self.anno2 = ""

	def get_target_id(self):
		return "{0};{1};{2}".format(self.t_trans_id, self.t_gene_id, self.t_gene_symbol)

	def get_lncRNA_id(self):
		return "{0};{1};{2}".format(self.l_trans_id, self.l_gene_id, self.l_gene_symbol)

	def toString(self):
		inter_string = ""
		for member in rna_inter_members:
			inter_string += str(getattr(self, member)) + "\t"
		return inter_string[:-1]

class Protein_Interaction:
	def __init__(self, data, correl):
		self.start = data[0]
		self.end = data[1]
		self.gene_id = data[2]
		self.gene_symbol = data[3]
		#self.peak = data[4]
		self.pvalue = data[5]
		#self.pid = data[6]
		self.correl = correl
		self.mutations = 0

	def toString(self):
		inter_string = ""
		for member in protein_inter_members:
			inter_string += str(getattr(self, member)) + "\t"
		return "{0}\t{1}".format(inter_string[:-1], self.correl.toString()) 

class Correlation:
	def __init__(self, data):
		self.pcor = data[0]
		self.pvalue = data[1]
		#self.dir_pvalue = data[2]
		#self.dir = data[3]

	def toString(self):
		correl_string = ""
		for member in correlation_members:
			correl_string += str(getattr(self, member)) + "\t"
		return correl_string[:-1]

#######################################################
#
# Utilities
#
#######################################################

def sort_key(item):
	return float(item.split(';')[14][:-1])

def sort_key_p(item):
	return float(item.split('\t')[35][:-1]) #(len(item)-1)

def sort_key_id(item):
	return float(item.split('\t')[0])

def merge(chunks,key=None):
	if key is None:
		key = lambda x : x

	values = []

	for index, chunk in enumerate(chunks):
		try:
			iterator = iter(chunk)
			value = iterator.next()
		except StopIteration:
			try:
				chunk.close()
				os.remove(chunk.name)
				chunks.remove(chunk)
			except:
				pass
		else:
			heappush(values,((key(value),index,value,iterator,chunk)))

	while values:
		k, index, value, iterator, chunk = heappop(values)
		yield value
		try:
			value = iterator.next()
		except StopIteration:
			try:
				chunk.close()
				os.remove(chunk.name)
				chunks.remove(chunk)
			except:
				pass
		else:
			heappush(values,(key(value),index,value,iterator,chunk))

def filter(repeatsfile, inputfile, outputfile):

	writer = open(outputfile, 'w')
	repeats = list()
	l = 0

	with open(inputfile, 'r') as data_in:
		line = data_in.readline()
		tokens = line.split(';')
		trans_id = tokens[1][0:15]

		with open(repeatsfile, 'r') as repeats_in:
			for line in repeats_in:

				tokens = line.split()
				l = l + 1

				if(tokens[0] == trans_id):
					print "Masking: " + tokens[1] + " " + tokens[2]
 					repeats.append(tokens[1:3])
				elif(len(repeats) > 0):
					break

	if(len(repeats) == 0): return False

	with open(inputfile, 'r') as data_in:
		for line in data_in:
			repeat = False;
			tokens = line.split(';')

			start = tokens[4]
			end = tokens[5]

			for reps in repeats:
				if((int(reps[0]) <= int(start) and int(reps[1]) >= int(start)) or (int(reps[0]) <= int(end) and int(reps[1]) >= int(end))):
					repeat = True

			if(not repeat):
				writer.write(line)

	return True


def batch_sort(input,output,key=None,buffer_size=32000,tempdirs=[]):
	if not tempdirs:
		tempdirs.append(gettempdir())
	
	input_file = file(input,'rb',64*1024)
	try:
		input_iterator = iter(input_file)
		
		chunks = []
		try:
			for tempdir in cycle(tempdirs):
				current_chunk = list(islice(input_iterator,buffer_size))
				if current_chunk:
					current_chunk.sort(key=key)
					output_chunk = file(os.path.join(tempdir,'%06i'%len(chunks)),'w+b',64*1024)
					output_chunk.writelines(current_chunk)
					output_chunk.flush()
					output_chunk.seek(0)
					chunks.append(output_chunk)
				else:
					break
		except:
			for chunk in chunks:
				try:
					chunk.close()
					os.remove(chunk.name)
				except:
					pass
			if output_chunk not in chunks:
				try:
					output_chunk.close()
					os.remove(output_chunk.name)
				except:
					pass
			return
	finally:
		input_file.close()
	
	output_file = file(output,'wb',64*1024)
	try:
		output_file.writelines(merge(chunks,key))
	finally:
		for chunk in chunks:
			try:
				chunk.close()
				os.remove(chunk.name)
			except:
				pass
		output_file.close()

def load_data(inputfile, top, sample):

	lookup = dict()
	data = list()
	energies = list()
	extra_energies = list()
	num_lines = 0
	count = 0

	with open(inputfile, 'r') as data_in:
		for line in data_in:
			num_lines = num_lines+1

	index = int(float(num_lines)*(float(top)/100))
	
	with open(inputfile, 'r') as data_in:
		for line in data_in:

			if(count < index):

				candidate = RNA_Interaction(line[:-1])
				candidate.mutations = get_mutation_count(candidate.t_trans_id, candidate.t_start, candidate.t_end)
				candidate.mutations = candidate.mutations + get_mutation_count(candidate.l_trans_id, candidate.l_start, candidate.l_end)

				if(lookup.get(candidate.t_trans_id) == None):
					lookup[candidate.t_trans_id] = list()
				lookup[candidate.t_trans_id].append(len(data))

				data.append(candidate)

				all_genes[candidate.t_gene_id] = True
				all_trans[candidate.t_trans_id] = True

				energies.append(float(candidate.E))

				if(count % sample == 0):
					extra_energies.append(float(candidate.E))

			elif(count == (index+1)):
				
				for energy in extra_energies:
					energies.append(energy)

			elif(count % sample == 0):

			 	energies.append(float(line.split(';')[14]))

			count = count+1

	return (lookup, data, energies, index)

def compute_pvalues(data, energies, index):

	eng_obs = [float(x) * -1 for x in energies[0:index]]

	if(index == len(energies)): index = 0 #Use whole data as background

	eng_back = [float(x) * -1 for x in energies[index:]]

	fit_alpha, fit_loc, fit_beta=stats.gamma.fit(eng_back)

	pvalues = stats.gamma.sf(eng_obs, fit_alpha, loc=fit_loc, scale=fit_beta)
	
	reject, pvalues_corr = fdrcorrection0(pvalues, is_sorted=True)

	for i in xrange(len(data)):
		data[i].pvalue = pvalues[i]
		data[i].fdr = pvalues_corr[i]

	# for i in xrange(0, len(pvalues_corr)):
	# 	print "{0}\t{1}\t{2}".format(str(energies[i]), pvalues[i], pvalues_corr[i]) 

def get_mutation_count(trans_id, start, end):
	
	count = 0

	if trans_id in indels:
		for indel in indels[trans_id]:
			if(int(indel[0]) > int(start) and int(indel[1]) < int(end)):
				count = count + 1

	return count

def get_mechanism(candidate, apriori):

	#print candidate
	protein_target = False
	protein_lncRNA = False

	if((candidate.t_peak.start != "NA")):
		protein_target = True

	if((candidate.l_peak.start != "NA")):
		protein_lncRNA = True

	if(apriori):
		if(protein_target):
			if( ((int(candidate.rna.t_start) >= (int(candidate.t_peak.start)-int(OVERLAP_DIST))) and  (int(candidate.rna.t_start) <= (int(candidate.t_peak.end)+int(OVERLAP_DIST)))) or ((int(candidate.rna.t_end) >= (int(candidate.t_peak.start)-int(OVERLAP_DIST))) and  (int(candidate.rna.t_end) <= (int(candidate.t_peak.end)+int(OVERLAP_DIST)))) ):
				if(protein_lncRNA == True):
					if(upregulators.get(candidate.t_peak.gene_id) != None and downregulators.get(candidate.l_peak.gene_id) != None):
						return "competitive_downregulation"
					elif(downregulators.get(candidate.t_peak.gene_id) != None and upregulators.get(candidate.l_peak.gene_id) != None):
						return "competitive_upregulation"
					else:
						return "nonsense"
				else:
					if(upregulators.get(candidate.t_peak.gene_id) != None):
						return "competitive_downregulation"
					elif(downregulators.get(candidate.t_peak.gene_id) != None):
						return "competitive_upregulation"
					else:
						return "nonsense"

			elif(  abs((int(candidate.t_peak.start) + (int(candidate.t_peak.end)-int(candidate.t_peak.start))/2) - (int(candidate.rna.t_start) + (int(candidate.rna.t_end)-int(candidate.rna.t_start))/2)) <= int(LOCAL_DIST)):
				if(protein_lncRNA == True):
					if(complexes[candidate.t_peak.gene_id].get(candidate.l_peak.gene_id) != None and downregulators.get(candidate.t_peak.gene_id) != None and downregulators.get(candidate.l_peak.gene_id) != None):
						return "complex_formation_downregulation"
					elif(complexes[candidate.t_peak.gene_id].get(candidate.l_peak.gene_id) != None and upregulators.get(candidate.t_peak.gene_id) != None and upregulators.get(candidate.l_peak.gene_id) != None):
						return "complex_formation_upregulation"
					else:
						return "nonsense"
				else:
					return "(de-)stabilization"
			else:
				return "nonsense"
		else:
			if(protein_lncRNA == True):
				if(downregulators.get(candidate.l_peak.gene_id) != None):
					return "localization_downregulation"
				elif(upregulators.get(candidate.l_peak.gene_id) != None):
					return "localization_upregulation"
				else:
					return "nonsense"
			else:
				return "direct_regulation"
	else:
		if(protein_target):
			#if( ((int(candidate.rna.t_start) >= int(candidate.t_peak.start)) and  (int(candidate.rna.t_start) <= int(candidate.t_peak.end))) or ((int(candidate.rna.t_end) >= int(candidate.t_peak.start)) and  (int(candidate.rna.t_end) <= int(candidate.t_peak.end))) ):
			if( ((int(candidate.rna.t_start) >= (int(candidate.t_peak.start)-int(OVERLAP_DIST))) and  (int(candidate.rna.t_start) <= (int(candidate.t_peak.end)+int(OVERLAP_DIST)))) or ((int(candidate.rna.t_end) >= (int(candidate.t_peak.start)-int(OVERLAP_DIST))) and  (int(candidate.rna.t_end) <= (int(candidate.t_peak.end)+int(OVERLAP_DIST)))) ):
				if(protein_lncRNA == True):
					if(float(candidate.correl.pcor) < 0 and float(candidate.t_peak.correl.pcor) > 0 and float(candidate.l_peak.correl.pcor) < 0 and upregulators.get(candidate.t_peak.gene_id) != None): #I am here
						return "competitive_downregulation"
					elif(float(candidate.correl.pcor) > 0 and float(candidate.t_peak.correl.pcor) < 0 and float(candidate.l_peak.correl.pcor) > 0 and downregulators.get(candidate.t_peak.gene_id) != None):
						return "competitive_upregulation"
					else:
						return "nonsense"
				else:
					if(float(candidate.correl.pcor) < 0 and float(candidate.t_peak.correl.pcor) > 0 and upregulators.get(candidate.t_peak.gene_id) != None):
						return "competitive_downregulation"
					elif(float(candidate.correl.pcor) > 0 and float(candidate.t_peak.correl.pcor) < 0 and downregulators.get(candidate.t_peak.gene_id) != None):
						return "competitive_upregulation"
					else:
						return "nonsense"
			else:
				if(protein_lncRNA == True):
					if(float(candidate.correl.pcor) < 0 and float(candidate.t_peak.correl.pcor) < 0 and float(candidate.l_peak.correl.pcor) < 0 and complexes[candidate.t_peak.gene_id].get(candidate.l_peak.gene_id) != None):
						return "complex_formation_downregulation"
					elif(float(candidate.correl.pcor) > 0 and float(candidate.t_peak.correl.pcor) > 0 and float(candidate.l_peak.correl.pcor) > 0 and complexes[candidate.t_peak.gene_id].get(candidate.l_peak.gene_id) != None):
						return "complex_formation_upregulation"
					else:
						return "nonsense"
				elif((downregulators.get(candidate.t_peak.gene_id) != None and float(candidate.t_peak.correl.pcor) < 0) or (upregulators.get(candidate.t_peak.gene_id) != None and float(candidate.t_peak.correl.pcor) > 0)):
					if(float(candidate.correl.pcor) < 0):# and downregulators.get(candidate.t_peak.gene_id) != None):
						return "de-stabilization"
					elif(float(candidate.correl.pcor) > 0):# and upregulators.get(candidate.t_peak.gene_id) != None):
						return "stabilization"
					else:
						return "nonsense"
				else:
					return "nonsense"
		else:
			if(protein_lncRNA == True):
				if(float(candidate.correl.pcor) <= 0 and float(candidate.l_peak.correl.pcor) < 0 and downregulators.get(candidate.l_peak.gene_id) != None):
					return "localization_downregulation"
				elif(float(candidate.correl.pcor) >= 0  and float(candidate.l_peak.correl.pcor) > 0 and upregulators.get(candidate.l_peak.gene_id) != None):
					return "localization_upregulation"
				else:
					return "nonsense"
			else:
				if(float(candidate.correl.pcor) < 0):
					return "direct_downregulation"
				elif(float(candidate.correl.pcor) > 0):
					return "direct_upregulation"
				else:
					return "nonsense"

def find_functional_domains(cancerfile, interactions, coor1, coor2, writer, type, bin_size, min_len):

	counts = dict()
	genes_back = set()
	cancer_genes = dict()
	c_count_back = 0
	n_count_back = 0
	m_count_back = 0
	w_count_back = 0
	max_bin = 0
	start = -1

	# with open(cancerfile, 'r') as data_in:
	# 	for gene in data_in:
	# 		cancer_genes[gene[:-1]] = True

	for inter in interactions:

		bin1 = int(float(getattr(inter, coor1))/float(bin_size))
		bin2 = int(float(getattr(inter, coor2))/float(bin_size))
		

		if(bin2 > max_bin): max_bin = bin2

		for i in xrange(bin1, bin2+1):
			if(i in counts):
				counts[i] = counts[i] + 1
			else:
				counts[i] = 1

		if(type != "Protein"): 
			genes_back.add(inter.t_gene_symbol)
			if(get_mutation_count(inter.t_trans_id, inter.t_start, inter.t_end) > 0):
				m_count_back = m_count_back + 1
			else:
				w_count_back = w_count_back + 1

	# for gene in genes_back:
	# 	if(gene in cancer_genes):
	# 		c_count_back = c_count_back + 1
	# 	else:
	# 		n_count_back = n_count_back + 1

	for i in xrange(0, len(counts)):
		if(i not in counts):
			counts[i] = 0

	counts_list = counts.values()
	counts_list.sort()
	arr = np.array(counts_list)
	mean = np.mean(arr)
	stdev = np.std(arr)
	exp = int(mean + 2*stdev)

	for i in xrange(0, max_bin+1):
		if(i in counts):
			counts[i] = counts[i] + 1
		else:
			counts[i] = 0

		if(int(counts[i]) >= int(exp)):   #TODO replace with stat test
			if(start == -1): 
				start = i
		else:
			d_start = (start)*bin_size
			d_end = (i)*bin_size

			if(start != -1 and int(d_end-d_start) >= int(min_len)): 
				if(type == "Protein"):
					proteins = set()

					for inter in interactions:
						if((int(getattr(inter, coor1)) >= int(d_start) and int(getattr(inter, coor1)) <= int(d_end)) or (int(getattr(inter, coor2)) >= int(d_start) and int(getattr(inter, coor2)) <= int(d_end))):
							proteins.add(inter.gene_id)

					writer.write(type + " DOMAIN:\t" + str(d_start) + "\t" + str(d_end) + "\t" + ",".join(proteins) + "\n")
				else:
					genes = set()
					m_count = 0
					w_count = 0

					for inter in interactions:
						if((int(getattr(inter, coor1)) >= int(d_start) and int(getattr(inter, coor1)) <= int(d_end)) or (int(getattr(inter, coor2)) >= int(d_start) and int(getattr(inter, coor2)) <= int(d_end))):
							genes.add(inter.t_gene_symbol)
							if(get_mutation_count(inter.t_trans_id, inter.t_start, inter.t_end) > 0):
								m_count = m_count + 1
							else:
								w_count = w_count + 1

					max_count = 0
					for i in xrange(d_start, d_end+1):
						if(int(counts[i]) > int(max_count)): max_count = counts[i]

					c_count = 0
					n_count = 0

					# for gene in genes:
					# 	if(gene in cancer_genes):
					# 		c_count = c_count + 1
					# 	else:
					# 		n_count = n_count + 1

					#print "Cancer: " + str(c_count) + "\t" + str(n_count) + "\t" + str(c_count_back) + "\t" + str(n_count_back)
					#print "Mutation: " + str(m_count) + "\t" + str(w_count) + "\t" + str(m_count_back) + "\t" + str(w_count_back)
					#oddsratio, c_pvalue = stats.fisher_exact([[c_count, n_count], [c_count_back, n_count_back]], alternative='greater')
					oddsratio, m_pvalue = stats.fisher_exact([[m_count, w_count], [m_count_back, w_count_back]], alternative='greater')

					#writer.write(type + " DOMAIN:\t" + str(d_start) + "\t" + str(d_end) + "\t" + str(max_count) + "\t" + str(c_pvalue) + "\t" + str(m_pvalue) + "\n")
					writer.write(type + " DOMAIN:\t" + str(d_start) + "-" + str(d_end) + "\tCount:\t" + str(max_count) + "\tP-Value:\t" + str(m_pvalue) + "\n")
			start = -1


#######################################################
#
# Main Code
#
#######################################################

def main(argv):

	try:
		opts, args = getopt.getopt(argv,"hi:r:m:n:c:p:x:a:",["ifile=","mfile="])#,"ofile="])
	except getopt.GetoptError:
		print 'infer.py -i <inputfile> -r <reference> -m <models> -n <mechanisms> -c <correlation> -p <peak-cutoff> -a <apriori>'
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			print 'infer.py -i <inputfile> -r <reference> -m <models> -n <mechanisms> -c <cancer> -p <peak-cutoff> -x <top-percent> -a <apriori>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-r", "--reference"):
			ref = arg
		elif opt in ("-m", "--models"):
			modelsfile = arg
		elif opt in ("-n", "--mechanisms"):
			mechsfile = arg
		elif opt in ("-c", "--correlation"):
			correlationfile = arg
		elif opt in ("-p", "--peak-cutoff"):
			p_cutoff = arg
		elif opt in ("-x", "--topX"):
			topx = arg
		elif opt in ("-a", "--apriori"):
			if(arg in ['True', 'true', 'T', 't', '1']):
				apriori = True
			else:
				apriori = False

	workdir = os.path.dirname(os.path.abspath(inputfile))
	proteinfile = workdir + "/data/peaks/graphprot.peaks.{0}.sig.pp25".format(ref) 
	annotationfile = workdir + "/data/anno.{0}".format(ref)
	indelsfile = workdir + "/data/indels.{0}".format(ref)
	repeatsfile = workdir + "/data/repeats.{0}".format(ref)
	upregulatorsfile = workdir + "/data/upregulators"
	downregulatorsfile = workdir + "/data/downregulators"
	complexesfile = workdir + "/data/complexes"
	cancerfile = workdir + "/data/cancer.genes"
	outputfile = inputfile + ".mechs"
	outputfiledomains = inputfile + ".domains"

	intarna_data = list()
	models = list()
	graphprot_data = dict()
	anno_data = dict()
	correl_lncRNA_data = dict()
	correl_proteins_data = dict()
	cur_transcript = None
	sampling = 10

	if(apriori):
		topx = 100
		sampling = 1

	with open(indelsfile, 'r') as indels_in:
		for line in indels_in:

			tokens = line.split()

			if tokens[0] not in indels:
				indels[tokens[0]] = list()
			indels[tokens[0]].append(tokens[1:3])

	print "Loading and sorting input file..."

	if(filter(repeatsfile, inputfile, inputfile+".tmp.filter")):
		batch_sort(inputfile+".tmp.filter", inputfile+".tmp.sort", key=sort_key)
	else:
		batch_sort(inputfile, inputfile+".tmp.sort", key=sort_key)

	lookup, data, energies, index = load_data(inputfile + ".tmp.sort", topx, sampling)

	os.remove(inputfile+".tmp.filter")
	os.remove(inputfile+".tmp.sort")
	
	lncRNA_ids = data[0].get_lncRNA_id()
	lncRNA = [data[0].l_trans_id, data[0].l_gene_id, data[0].l_gene_symbol]
	all_genes[data[0].l_gene_id] = True
	all_trans[data[0].l_trans_id] = True

	print "Computing RNA-RNA interaction p-values..."

	compute_pvalues(data, energies, index)

	print "Annotating RNA-RNA interactions..."

	with open(annotationfile, 'r') as anno_in:
		for line in anno_in:
			tokens = line.split('\t')
			tokens[2] = tokens[2][:-1]
			anno_data[tokens[0]] = tokens[1:3]

	for RR_inter in data:

		if(anno_data[RR_inter.t_trans_id][0] == 0):
			RR_inter.anno1 = "non-coding"
			RR_inter.anno2 = "non-coding"
		else:
			if(RR_inter.t_start <= anno_data[RR_inter.t_trans_id][0]):
				RR_inter.anno1 = "3'UTR"
			elif(RR_inter.t_start >= anno_data[RR_inter.t_trans_id][1]):
				RR_inter.anno1 = "5'UTR"
			else:
				RR_inter.anno1 = "CDS"

			if(RR_inter.t_end <= anno_data[RR_inter.t_trans_id][0]):
				RR_inter.anno2 = "3'UTR"
			elif(RR_inter.t_end >= anno_data[RR_inter.t_trans_id][1]):
				RR_inter.anno2 = "5'UTR"
			else:
				RR_inter.anno2 = "CDS"

	print "Loading auxiliary files..."

	with open(upregulatorsfile, 'r') as up_in:
		for line in up_in:
			upregulators[line[:-1]] = True

	with open(downregulatorsfile, 'r') as down_in:
		for line in down_in:
			downregulators[line[:-1]] = True

	with open(complexesfile, 'r') as complex_in:
		for line in complex_in:
			if(line[0] == "#"):continue

			tokens = line.split()

			if(len(tokens) == 1):
				complexes[line[:-1]] = dict()
			else:
				complexes[tokens[0]][tokens[1]] = True

	with open(modelsfile, 'r') as model_in:
		for model in model_in:
			model = model[:-1]
			models.append(model)
			all_genes[model.split(";")[0]] = True
			correl_proteins_data[model] = dict()
			graphprot_data[model] = dict()
		model = "NA;NA"
		models.append(model)
		correl_proteins_data[model] = dict()
		graphprot_data[model] = dict()

	with open(mechsfile, 'r') as mechs_in:
		for mech in mechs_in:
			mech = mech[:-1]
			mechs.add(mech)

	if(not apriori):
		print "Loading correlation file..."		
		with open(correlationfile, 'r') as correl_in:
			for line in correl_in:
				tokens = line.split()
				if(tokens[0] not in all_genes or tokens[1] not in all_genes): continue

				if(tokens[0] == lncRNA[1]):
					correl_lncRNA_data[tokens[1]] = Correlation(tokens[2:6])
				if(tokens[1] == lncRNA[1]):
					correl_lncRNA_data[tokens[0]] = Correlation(tokens[2:6])
				for model in models:
					protein_id = model.split(";")[0]
					protein_gene = model.split(";")[1]
					if(tokens[0] == protein_id):
						correl_proteins_data[model][tokens[1]] = Correlation(tokens[2:6])
					elif(tokens[1] == protein_id):
						correl_proteins_data[model][tokens[0]] = Correlation(tokens[2:6])
					elif(protein_id == "NA"):
						correl_proteins_data[model][tokens[0]] = Correlation(["0","1","1","NA"])
						correl_proteins_data[model][tokens[1]] = Correlation(["0","1","1","NA"])

	print "Loading protein file..."
	with open(proteinfile, 'r') as graphprot_in:
		for line in graphprot_in:
			tokens = line.split()
			ids = tokens[0].split(";")
			trans = ids[0]
			gene = ids[1]

			if(trans not in all_trans): continue

			for model in models:
				if(model == "NA;NA"):
					if tokens[0] not in graphprot_data[model]:
						graphprot_data[model][tokens[0]] = list()
						graphprot_data[model][tokens[0]].append(Protein_Interaction(["NA","NA","NA","NA","0","1","NA"], Correlation(["0","1","1","NA"])))
				else:
					for protein in tokens[3].split(","):
						protein_id = protein.split(";")[0]
						protein_gene = protein.split(";")[1]
						if(protein == model):
							peak_data = tokens[1:]
							peak_data[2] = protein_id
							peak_data.insert(3, protein_gene)
							peak = Protein_Interaction(peak_data, Correlation(["0","1","1","NA"]))
							peak.mutations = get_mutation_count(trans, peak.start, peak.end)

							if gene in correl_proteins_data[model]:
								peak.correl = correl_proteins_data[model][gene]

							if tokens[0] not in graphprot_data[model]:
								graphprot_data[model][tokens[0]] = list()
							graphprot_data[model][tokens[0]].append(peak)


	print "Inferring mechanisms..."

	writer = open(outputfile, 'w')

	for key in lookup.keys():

		indices = lookup[key]

		gene = data[indices[0]].t_gene_id
		target_ids = data[indices[0]].get_target_id()
		top_candidate = None

		for index in indices:

			RR_inter = data[index]

			for protein_target in models:
				pt_correl = correl_proteins_data[protein_target].get(gene)
				pt_peaks = graphprot_data[protein_target].get(target_ids)
				if(protein_target == "NA;NA"):
					pt_correl = Correlation(["0","1","1","NA"])
					pt_peaks = list()
					pt_peaks.append(Protein_Interaction(["NA","NA","NA","NA","0","1","NA"], Correlation(["0","1","1","NA"])))
				if((not apriori and (pt_correl == None)) or pt_peaks == None): continue

				for peak_pt in pt_peaks:
					for protein_lncRNA in models:
						pl_correl = correl_proteins_data[protein_lncRNA].get(gene)
						pl_peaks = graphprot_data[protein_lncRNA].get(lncRNA_ids)
						if(protein_lncRNA == "NA;NA"): 
							pl_correl = Correlation(["0","1","1","NA"])
							pl_peaks = list()
							pl_peaks.append(Protein_Interaction(["NA","NA","NA","NA","0","1","NA"], Correlation(["0","1","1","NA"])))
						if((not apriori and (pl_correl == None or (protein_target == protein_lncRNA and protein_target != "NA;NA"))) or pl_peaks == None): continue

						for peak_pl in pl_peaks: 						
							candidate = Candidate(RR_inter)
							candidate.correl = (correl_lncRNA_data.get(gene, Correlation(["0","1","1","NA"])))
							candidate.t_peak = peak_pt
							candidate.l_peak = peak_pl

							mechanism = get_mechanism(candidate, apriori)
							
							if(mechanism != "nonsense" and mechanism in mechs):

								candidate.mechanism = mechanism
								pvalues = candidate.get_pvalues()

								if(len(pvalues) > 1): 
									xsq, p = stats.combine_pvalues(np.array(pvalues), method='stouffer', weights=None)
									candidate.joint_pvalue = p
								else:
									candidate.joint_pvalue = pvalues[0]

								if(top_candidate == None or float(candidate.joint_pvalue) < float(top_candidate.joint_pvalue)):
									top_candidate = candidate
									
		if(top_candidate != None):
			writer.write(top_candidate.toString() + "\n")

	writer.close()
	
	batch_sort(outputfile, (outputfile + ".sort.tsv"), key=sort_key_p)

	header = "target_transcript_ID\ttarget_gene_ID\ttarget_gene_symbol\tlncRNA_transcript_id\tlncRNA_gene_ID\tlncRNA_gene_symbol\ttarget_start\ttarget_end\tlncRNA_start\tlncRNA_end\tfree_energy\tP_value\tFDR\tmutation_count\tleft_context\tright_context\tpartial_correl\tpartial_correl_fdr\trbp_target_peak_start\trbp_target_peak_end\trbp_gene_ID\trbp_gene_symbol\trbp_p_value\tmutation_count\tpartial_correl\tpartial_correl_fdr\trbp_lncRNA_peak_start\trbp_lncRNA_peak_end\trbp_gene_ID\trbp_gene_symbol\trbp_p_value\tmutation_count\tpartial_correl\tpartial_correl_fdr\tmechanism\tjoint_p_value"

	with open((outputfile + ".sort.tsv"), 'r+') as f:
		content = f.read()
		f.seek(0, 0)
		f.write(header + '\n' + content)


	#==================================
	#= Domain Analysis
	#==================================

	print "Predicting functional domains..."

	writer = open(outputfiledomains, 'w')

	find_functional_domains(cancerfile, data, "l_start", "l_end", writer, "RNA-RNA", 1, 20)

	peaks = list()

	for model in models:
		if(model != "NA;NA"):
			if (graphprot_data[model].get(lncRNA_ids) != None):
				peaks.extend(graphprot_data[model].get(lncRNA_ids)) 
	
	find_functional_domains(cancerfile, peaks, "start", "end", writer, "Protein", 1, 20)
	
	writer.close()

	print "Done"

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))


#"Gene_Name,Accessible_Region,Target_Start,Target_End,lncRNA,lncRNA_Start,lncRNA_End,Free_Energy,LT_P-Value,
#Protein_lncRNA,Protein_Start,Protein_End,PL_P-Value,Protein_Target,Protein_Start,Protein_End,PT_P-Value"
#LT_Correl, LT_Correl_P-Value, PL_Correl, PL_Correl_P-Value, PT_Correl, PT_Correl_P-Value, Inferred_Mechanism, Combined_P-value
