#!/bin/bash

traj="multimodel.pdb"
#pdb="../../../../md/c22star/9-md/Mol_An/update.gro"
pdb="model_0.pdb"
mutsite="M75"
r2="A-75MET"
#prepare folders
mkdir V48 L52 V74 R88 T91 F95

#mutation paths
cd V48 
mkdir $mutsite
cd $mutsite
graph_analysis -r ../../../$pdb -a ../../../contact_filter.dat  -p -r1 A-48VAL -r2 $r2 -l 10  > log_path_48V_$mutsite 
cd ../../L52
mkdir $mutsite
cd $mutsite
graph_analysis -r ../../../$pdb -a ../../../contact_filter.dat  -p -r1 A-52LEU -r2 $r2 -l 10  > log_path_52LEU_$mutsite 
cd ../../V74
mkdir $mutsite
cd $mutsite
graph_analysis -r ../../../$pdb -a ../../../contact_filter.dat  -p -r1 A-74VAL -r2 $r2 -l 15  > log_path_74V_$mutsite   
cd ../../R88
mkdir $mutsite
cd $mutsite
graph_analysis -r ../../../$pdb -a ../../../contact_filter.dat  -p -r1 A-88ARG -r2 $r2 -l 15  > log_path_88R_$mutsite 
cd ../../T91
mkdir $mutsite
cd $mutsite
graph_analysis -r ../../../$pdb -a ../../../contact_filter.dat  -p -r1 A-91THR -r2 $r2 -l 15  > log_path_91T_$mutsite 
cd ../../F95
mkdir $mutsite
cd $mutsite
graph_analysis -r ../../../$pdb -a ../../../contact_filter.dat  -p -r1 A-95PHE -r2 $r2 -l 15  > log_path_95PHE_$mutsite 
