#!/usr/bin/env python3
### Quick and dirty program to extract mutations from the cosmic database for a list of genes specifying a tissue of interest
### Usage: Python3 retrieve_mutations_cosmic.py tissue gene_list.txt outfile.csv

import sys, os, re, csv, argparse
from collections import defaultdict

usage = '''This program retrieves mutation from the COSMIC database for a specified tissue and set of genes '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 1.0'
    )
parser.add_argument(
	'database',
	help='path to cosmic database',
    )
parser.add_argument(
    'tissue',
    help='tissue to query',
    )
parser.add_argument(
    'gene_file',
    help='file with genes to query',
    type=argparse.FileType('r'),
    )
parser.add_argument(
    'outfile',
    help='file(csv) to print results to',
    type=argparse.FileType('w'),
    default=sys.stdout
    )

args = parser.parse_args()
parser.parse_args()



#path = '/data/databases/cosmic-v84/'
path = str(sys.argv[1])

# Search all variants of COSMIC mutation databases
databases = ['CosmicMutantExport.tsv',
         'CosmicCompleteTargetedScreensMutantExport.tsv',
         'CosmicGenomeScreensMutantExport.tsv']

tissue = str(sys.argv[2]).lower()

mutations = defaultdict(set)
with open(sys.argv[3], 'r') as f_in, open(sys.argv[4], 'w') as f_out:
	print('{},{},{},{},{}'.format('Gene name', 'Primary Site', 'Mutation AA', 'Mutation CDS', 'Mutation Description'), file=f_out) 
	gene_names = [line.rstrip().upper() for line in f_in]

	for gene in gene_names:
		gene_reg = re.compile(gene+'\s')
		for file in databases:
			with open(os.path.join(path, file), 'r') as cosmic_db:
				for line in cosmic_db:
					line = line.rstrip()
					if gene_reg.match(line):
						gene_name = line.split('\t')[0]
						tis = line.split('\t')[7]
						mut_aa = line.split('\t')[18]
						mut_cds = line.split('\t')[17]
						mut_des = line.split('\t')[19]
		
						if tis == tissue and mut_des == 'Substitution - Missense':
							mutations[gene_name].add(mut_aa)
							print('{},{},{},{},{}'.format(gene_name, tissue, mut_aa, mut_cds, mut_des), file=f_out)
							

'''
	for key, value in mutations.items():
		print('{}\t{}'.format(key, ','.join(value)), file = f_out)
'''


