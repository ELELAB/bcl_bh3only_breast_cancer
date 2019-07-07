#!/bin/bash

./retrieve_mutations_cosmic.py /data/databases/cosmic-v84/ breast genes_of_int.txt mutations_cosmic.csv

./retrieve_mutations_cbio.R genes_of_int.txt



