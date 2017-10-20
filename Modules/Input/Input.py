#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# 
# www.genesilico.pl 
#

#creates ranked 3D models of macromoleular complexes 
#based on experimental restraints and a whole complex shape.


__author__ = "Joanna M. Kasprzak"
__copyright__ = "Copyright 2010, The PyRy3D Project"
__credits__ = ["Janusz Bujnicki"]
__license__ = "GPL"
__version__ = "0.1.0"
__maintainer__ = "Joanna Kasprzak"
__email__ = "jkasp@amu.edu.pl"
__status__ = "Prototype"

from InputChecker           import DataChecker, DisorderDataChecker
from ArgumentsParser        import get_args
from Modules.Error.Errors   import InputError

class Input(object):
    """
        represents Input Data 
    """
    def __init__(self):
        self.components_with_structure = []
        self.components_no_struct      = []
        self.sequences                 = None #all fasta sequences given by the user (Bio objects)
        self.structures                = None #all structures given by the user (Bio.PDB objects)
        self.restraints                = None #restraints encoded as Filtrest3D objects
        self.symrestraints             = None #Filtrest3D symmetry restraints (definition of fragments with equal lengths)
        self.mapfile                   = None #a density map from EM
        self.saxsfile                  = None #a shape file from SAXS (eg dammin/f output file)
        self.curvefile                 = None #a file with experimental SAXS/SANS curve
        self.outname                   = ''   #name of simulation output folder
        self.config                    = ''   #name of file with simulation parameters
        self.traflfile                 = None #filename with trajectory if required by the user
        self.fullatom_file             = ""   #filename with fullatom representation for best complex generated during all simulation
        self.movehistory_file          = ""   #filename containing history of moves for all saved complexes
        self.scoreplot                 = ""
        self.optimized				   = False #when True, program is using c++ version of computing results
        self.autofolder                = None #folder with structures for automatic generation of input files
        
    def str(self):
        return "%s " % (self.chains)
        
    def check_fastanames_duplications(self):
        """
        checks whether names of sequences in FASTA file are not duplicated.
        If so it raises error and asks to change duplicated name!
        """
        components_names = []
        for fastaseq in self.sequences:
            fasta_name = fastaseq.id.split("_")
            if len(fasta_name) == 1:
                if fastaseq.id in components_names: raise InputError(fastaseq.id+\
                " name has been repeated twice. please provide other name for this component")
                components_names.append(fastaseq.id)
            else:
                if fasta_name[0] in components_names: raise InputError(fasta_name[0]+\
                " name has been repeated twice. please provide other name for this component")
                components_names.append(fasta_name[0])
        
    def check_input(self, component):
        """
        for a particular complex component runs checking data procedures by
        calling DataChecker module
        Parameters:
        -----------
            component : complex component (object)
            components_names : already taken chain names in fasta headers
        Returns:
        ---------
            Should assign self.chains and self.chain, self.str_seq atrributes
        """
        checked = DataChecker()
        pyrycomponent = checked.check_component(self.sequences, component)
        fastasequence = checked.fasta_seq
        self.components_with_structure.append(pyrycomponent.chain)
        return pyrycomponent, fastasequence
    
    def check_seq_no_struct_input(self,seq_no_struct):
        """
        """
        checked = DisorderDataChecker()
        pyrystruct = checked.check_dd_component(seq_no_struct)
        return pyrystruct
        
    def get_data(self):
        """
        method calls all functions responsible for input data gathering
        by ArgumentsParser module
        Parameters:
        -----------
            input files provided by the user (sequences, structures,
                                            restraints, density map)
        Returns:
        ---------
            self.sequences  
            self.structures 
            self.mapfile    
            self.restraints 
        """
        opts = get_args(self)
        return opts
        
    def get_seq_with_no_struct(self):
        """
        finds all sequences that do not have structure coordinates defined by the user;
        these elements should have an estimated volume assigned
        Parameters:
        -----------
            works on Input object only
        Returns:
        --------
            components_no_structure   : complex components for which user
                                        haven't defined a structure
        """
        [self.components_no_struct.append(seq) for seq in self.sequences \
                                  if seq.name not in self.components_with_structure]
        return self.components_no_struct





