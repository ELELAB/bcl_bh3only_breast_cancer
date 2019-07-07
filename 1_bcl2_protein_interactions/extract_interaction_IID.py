#!/usr/bin/env python3

# This program: (1) Extract protein interaction partners of the Bcl-2 family members from IID database, (2) retrieve [fasta](http://www.uniprot.org), (3) filter out PPI's containing the BH3 motif and (4) output tab-delimited file for [Cytoscape](http://www.cytoscape.org) containing source and target interactors.
# sys.argv[1] file containing the query ID's
# sys.argv[2] database file from IID 
# sys.arvv[3] output file 

import sys
import re
from urllib.request import urlopen
import argparse
from collections import defaultdict

usage = ''' This program: (1) Extract protein interaction partners of the Bcl-2 family members from a IID database file, (2) retrieve [fasta](http://www.uniprot.org), (3) filter out PPI's containing the BH3 motif, 
(4) outputs tab-delimited file containing source and target interactors (gene names),(5) outputs a file of all interaction partners (gene names) and (6) outputs a fasta file of interactions partners.
 '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 1.0'
    )
parser.add_argument(
	'query_file',
	help='file of uniprot IDs',
    type=argparse.FileType('r'),
    )
parser.add_argument(
    'ppi_file',
    help='IID database file',
    type=argparse.FileType('r'),
    )
parser.add_argument(
    'outfile_1',
    help='BC2-2 family interactions',
    type=argparse.FileType('w'),
    default=sys.stdout
    )
parser.add_argument(
    'outfile_2',
    help='All interantion partners - candidates',
    type=argparse.FileType('w'),
    default=sys.stdout
    )
parser.add_argument(
    'outfile_3',
    help='Fasta file of protein sequences for candidates',
    type=argparse.FileType('w'),
    default=sys.stdout
    )

args = parser.parse_args()
parser.parse_args()




with open(sys.argv[1], 'r') as BCL_2_file, open(sys.argv[2], 'r') as iid_database, open(sys.argv[3], 'w') as f_out, open(sys.argv[4], 'w') as f_out2, open(sys.argv[5], 'w') as f_out3:

	BCL_2 = [line.rstrip() for line in BCL_2_file]
	unique_ids = set()
	fasta = defaultdict(set)
	for line in iid_database:
		if not line.startswith('Query'):	
			query = line.rstrip().split()[1]
			partner = line.rstrip().split()[2]
			query_sym = line.rstrip().split()[3]
			partner_sym = line.rstrip().split()[4]
			if query in BCL_2:


				# Fetch the protein's fasta file and get rid of newlines. Check if all accessions are found on uniprot
				uniprot = 'http://www.uniprot.org/uniprot/'+partner+'.fasta'
				try:
					file = urlopen(uniprot).read().decode('utf-8')
					
					header = ''.join(file.split('\n')[:1])
					sequence = ''.join(file.split('\n')[1:])

					#output fasta formatted file before motif search if desirable
					#print(header, file=f_out) 
					#print(sequence, file = f_out)

					motif_match = re.search(r'((A|V|I|L|M|F|Y|W|C|P).{3,4}L.{2,3}(A|V|I|L|M|F|Y|W|C|P)(G|A|S|C).?(D|E|Q|N))', sequence) #  motifs L.{3}(G|A)D L.{4}D
					
					if motif_match:
						print('{}\t{}'.format(query_sym, partner_sym), file = f_out)
						unique_ids.add(partner_sym)
						unique_ids.add(query_sym)
						fasta[header].add(sequence)
						
						
						
						#output fasta formatted file:
						#print(header, file=f_out3)
						#print(sequence, file=f_out3)
						
				# Give message if the entry is not found at uniprot
				except:
					print('{}{}'.format('entry is obsolete in the uniprot database: ',partner))
					pass

	# output unique candidates
	for i in unique_ids:
		print(i, file = f_out2)
	for head, seq in fasta.items():
		print(head, file = f_out3)
		print(''.join(seq), file = f_out3)

	

				


