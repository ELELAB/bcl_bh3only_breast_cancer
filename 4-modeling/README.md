# Modelling

Protein-peptide interactions will be modelled with the scope of: i) predicting the binding interface and the three-dimensional structure of their complex, ii) identifying the mutations that are mapping at the interface or can communicate long-range with them.

To model protein-peptide interaction we applied comparative modelling by satisfaction of spacial restraints, implemented in the program MODELLER  [@Sali1993]. 


[MODELLER](https://salilab.org/modeller/) carries out comparative protein structure modelling by statisfaction of spacial restraints including (1) homology-derived restraints on the distances and dihedral angles in the target sequence, extracted from its alignment with the template structures [@Sali1993], (2) stereochemical restraints on bond length and bond angle preferences, obtained from the CHARMM22 molecular mechanics forcefield [@MacKerell1998], (3) statistical preferences for dihedral angles and non-bonded inter-atomic distances, obtained from a representative set of known protein structures [@Sali1994, @Sali2006], and (4) optional restraints. The model is subsequently obtained by minimizing the violations of all restraints, archieved by real-space optimization [@Marti-Renom2000].
We restrained the CA-CA distance of our model to 7 (st.dev.=0.1) angstroms between Val74 in chain A (Bcl2a1) and the highly concerved Leu in the BH3 motif in chain B (Hrk). Val74 is located in the hydrophopic groove acting as the interaction site on the template complex.


## Modelling Hrk

As a template structure we used the protein sequence of [Bcl2a1: Bcl-2-related protein A1](http://www.uniprot.org/uniprot/Q16548) in complex with the BH3 motif of [Puma: Bcl-2-binding component 3](http://www.uniprot.org/uniprot/Q9BXH1), experimentally derrived from X-RAY DIFFRACTION [PDB 5uul](https://www.rcsb.org/structure/5uul).
This complex was alligned to the target sequence complex between Bcl2a1 and the predicted BH3 regions from the [Hrk: Activator of apoptosis harakiri](http://www.uniprot.org/uniprot/O00198) (extended to correspond to the length of the Puma BH3 region in the template).


We computed 10 models for each alignments. To infer reliability of the models we implemented the SOAP-Peptide potential, trained specifically for protein-peptide complexes. 

## To reproduce this analysis

### Prerequisites

MODELLER is available free of charge to academic non-profit institutions; you will, however, need to register for a license in order to use the software. See [Download and Installation](https://salilab.org/modeller/download_installation.html).

### Running the analysis
```
$ cd bcl2_bh3_breast_cancer/modelling/output/
$ chmod +x run_models.sh
$ ./run_models.sh
```


