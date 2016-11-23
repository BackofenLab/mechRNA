import sys

def main(argv):

	reference_file = argv[0]
	protein_file = argv[1]
	ensembl_id = argv[2]
	motif_size = int(argv[3])


	graphprot_data = dict()

	with open(protein_file, 'r') as protein_in:
		for line in protein_in:
			tokens = line.split()
			trans = tokens[0]

			if trans not in graphprot_data:
				graphprot_data[trans] = list()
				graphprot_data[trans].append(tokens)
			else:
				graphprot_data[trans].append(tokens)


	with open(reference_file, 'r') as reference_in:
		for line in reference_in:
			if(line[0] == ">"):
				end = len(next(reference_in))
				header = line[1:]
				tokens = header.split(";")
				trans = tokens[0]
				gene = tokens[1]
				gene_sym = tokens[2]
				
				if trans in graphprot_data:
					for peak in graphprot_data[trans]:
						loc = int(peak[1])
						motif_start = loc-motif_size
						motif_end = loc+motif_size

						if(motif_start < 1):
							motif_start = 1

						if(motif_end > end):
							motif_end = end

						print header[:-1] + "\t{0}\t{1}\t".format(motif_start, motif_end) + ensembl_id + "\t{0}\t{1}".format(peak[2], peak[5])

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))

