#!/usr/bin/env python3
import re, sys, argparse



usage = ''' This program: (1) Takes a fasta file of protein sequences, (2) extracts the region(s) matching the BH3 motif and (3) outputs to file .  
 '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s 1.0'
    )
parser.add_argument(
	'in_file',
	help='fasta file of protein sequences',
    type=argparse.FileType('r'),
    )
parser.add_argument(
    'outfile_file',
    help='name of file to write to',
    type=argparse.FileType('w'),
    )

args = parser.parse_args()
parser.parse_args()

with open(sys.argv[1], 'r') as f_in, open(sys.argv[2], 'w') as f_out:
	for line in f_in:
		line = line.rstrip()
		if line.startswith('>'):
			print(line, file = f_out)
		else:
			sequence = line
			matches = re.finditer(r'((A|V|I|L|M|F|Y|W|C|P).{3,4}L.{2,3}(A|V|I|L|M|F|Y|W|C|P)(G|A|S|C).?(D|E|Q|N))', sequence) 
			for m in matches: 
			    base = m.group() 
			    start  = m.start() 
			    stop = m.end()
			    print('BH3 motif: ' + base + " found at position " + str(start) + '-' + str(stop), file = f_out)
			    print('Length of sequence: ' + str(len(sequence)), file = f_out)
		

