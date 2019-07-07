#!/bin/bash

source /usr/local/gromacs-5.1.5/bin/GMXRC.bash
#gmx_mpi editconf -f ../../../../md/c27/9-md/Mol_An/update.gro -o update.pdb
traj="multimodel.pdb"
#pdb="../../../../md/c22star/9-md/Mol_An/update.gro"
pdb="model_0.pdb"

#as a default it will use all hbonds
#all
pyinteraph -v -s $pdb -t $traj -r $pdb -f --hc-graph contact.dat --hc-co 5.0 --ff-masses encads --hc-residues ALA,CYS,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,HIS,LYS,ARG,GLU,ASP,ASN,GLN,SER,THR
