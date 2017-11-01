#!/usr/bin/env python
#-*- encoding:utf-8 -*-
import os, sys, operator, math, time


def usage():
	print '\nUsage: python get_contig.py fasta output line_length'
	sys.exit(-1)

####################
def LoadFasta( fasta_file ):
	ref_list = [] # keep order of reference in original input file
	ref_dict = {}
	tmp_list = []
	ref_id   = ""
	sr = open(fasta_file, "r")
	for line in sr:
		if ( ">" == line[0]):
			if ( ref_id != "" ):
				ref_dict[ ref_id] ="".join(tmp_list)
			ref_id   = line.strip()[1:]#.split()[0][1:]
			if ref_id not in ref_dict:
				ref_list.append( ref_id )
			del tmp_list[:]
		else:
			tmp_list.append( line.strip() )
	
	# adding records for the last chromosome
	if ( ref_id != "" ):
		ref_dict[ ref_id ] = "".join(tmp_list)

	sr.close()

	return ref_dict, ref_list

####################
def main():
	args = sys.argv[1:]
	if len(args) == 2:
		max_line = 100
	elif len(args) == 3:
		max_line = int(sys.argv[3])
	else:
		usage()


	start_time = time.time()

	#ref_dict, ref_list = LoadFasta0( sys.argv[1] )
	ref_dict, ref_list = LoadFasta( sys.argv[1] )
	print("--- Loading %s seconds ---" % (time.time() - start_time))
	start_time = time.time()
	sw = open( sys.argv[2], 'w')
	#for x in ref_list:
	for x,y in ref_dict.iteritems():
		y = ref_dict[x]
		#sw.write(">%s\n%s\n" % (x, y))
		sw.write(">%s\n" % (x))
		for i in range(0, len(y), max_line):
			sw.write("%s\n" % (y[i:i+max_line]))
	sw.close()
	print("--- Writing %s seconds ---" % (time.time() - start_time))
	
		
#############################################################################################
if __name__ == "__main__":
    sys.exit(main())

