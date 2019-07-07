#!/bin/bash

traj="multimodel.pdb"
#pdb="../../../../md/c22star/9-md/Mol_An/update.gro"
pdb="model_0.pdb"

#as a default it will use all hbonds
filter_graph -d contact.dat -o contact_filter.dat -t 20
graph_analysis -r $pdb -a contact_filter.dat -c -u -k 3 -cb ccomp_contact.pdb -ub hubs_contact.pdb > log_graph_analysis                                                                             
