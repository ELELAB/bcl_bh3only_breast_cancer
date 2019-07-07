# Identification of BCL-2 interaction partners containing the BH-3 motif

We exploited the [Integrated Interactions Database](http://iid.ophid.utoronto.ca/iid/About/) of tissue – and organism specific interactions to retrieve known interactions partners of the globular Bcl-2 family members. Subsequently we filtered the interaction partners to retain only those encompassing the BH3 motif. Finally we extracted region(s) of the protein sequences matching the BH3 motif. The interaction partners encompassing the BH3 motif we used as candidate genes in the downstream analysis.

## To retrieve protein-protein interaction (PPI) partners

### Prerequisites

```
Python version 3.0 or higher
```

The [Integrated interactions database] needs to be downloaded before running the scripts. We downloaded the database file of PPI for the 13 globular BCL2 family members with the following criteria:   

* Query genes: P10415, Q07817, Q92843, Q07820, Q9HD36, Q16548, Q9UMX3, Q07812, Q16611, Q9HB09, Q9BXK5, Q9BZR8, Q5TBC7
* Selected species: Human
* Selected tissues: Any 
* Find interaction partners supported by: Experimental evidence
* Required evidence: gene OR protein expression

To download the database file used by ```extract_interaction_IID.py``` see [Data Downloads](http://cancer.sanger.ac.uk/cosmic/download).

The database was downloaded 6 February 2018 - version 2017-04 (latest)

### Running the analysis
	
*To run ```extract_interaction_IID.py```:*

```
This program: (1) Extract protein interaction partners of the Bcl-2 family
members from a IID database file, (2) retrieve
[fasta](http://www.uniprot.org), (3) filter out PPI's containing the BH3
motif, (4) outputs tab-delimited file containing source and target interactors
(gene names),(5) outputs a file of all interaction partners (gene names) and
(6) outputs a fasta file of interactions partners.

positional arguments:
  query_file     file of uniprot IDs
  ppi_file       IID database file
  outfile_1      BC2-2 family interactions
  outfile_2      All interantion partners - candidates
  outfile_3      Fasta file of protein sequences for candidates

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit

```

```
$ cd bcl2_bh3_breast_cancer/bcl2_protein_interactions/
$ chmod +x extract_interaction_IID.py
$ ./extract_interaction_IID.py data/BCL2.queries.txt data/PPIs.txt output/bcl2.interactions.txt output/candidate.genes.txt output/candidate.genes.faa 
```

3. *To run ```find_bh3.py```:*
```
This program: (1) Takes a fasta file of protein sequences, (2) extracts the
region(s) matching the BH3 motif and (3) outputs to file .

positional arguments:
  in_file        fasta file of protein sequences
  outfile_file   name of file to write to

optional arguments:
  -h, --help     show this help message and exit
  -v, --version  show program's version number and exit
```


```
$ cd bcl2_bh3_breast_cancer/bcl2_protein_interactions/
$ chmod +x find_bh3.py 
$ ./find_bh3.py output/candidate.genes.faa output/bh3.motif.candidate.genes.txt
```


