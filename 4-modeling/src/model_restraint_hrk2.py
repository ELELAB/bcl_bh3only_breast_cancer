# Addition of restraints to the default ones
from modeller import *
from modeller.automodel import *    # Load the automodel class
from modeller import soap_peptide

log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../data']

class MyModel(automodel):
     def special_patches(self, aln):
         self.rename_segments(segment_ids=['A','B'], renumber_residues=[1, 28])

     def special_restraints(self, aln):
         rsr = self.restraints
         at = self.atoms
#
#       Restrain the specified CA-CA distance to 7 angstroms (st. dev.=0.1)
#       Use a harmonic potential and X-Y distance group.
         rsr.add(forms.gaussian(group=physical.xy_distance,
                                 feature=features.distance(at['CA:74:A'],		# Val74 on Bcl2a1 (PDB numbering)
                                                          at['CA:37:B']),	# Leu37 on Hrk - Leu79 on PUMA (PDB numbering)
                                mean=7.0, stdev=0.1))

a = MyModel(env,
            alnfile  = '../data/5uulpuma_bcl2a1hrk_2.ali',    # alignment filename
            knowns   = '5uul_temp.pdb',    # known template structure
            sequence = 'bcl2a1_hrk_2',     # target sequence
	    assess_methods=soap_peptide.Scorer())  # assess with SOAP

a.starting_model= 1                 # index of the first model
a.ending_model  = 10                # index of the last model
                                    # (determines how many models to calculate)
a.make()                            # do comparative modeling
