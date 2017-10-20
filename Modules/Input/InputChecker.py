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


#pyry
from Modules.Input.PyRyStructure import PyRyStructure
from Modules.Error.Errors        import InputError
from Modules.Constans.Constans   import AMINOACIDS, NUCLEOTIDES, ALLOWED_CHAIN_NAMES

class DataChecker(object):
    """
        checks if data defined by the user are consistent :-)
    """
    
    def __init__(self):
        self.already_taken = {}    #dict for gathering taken positions from given sequence
        self.struct_seq = ""       #sequence of analysed component; as in pdb file
        self.fasta_seq = None       #sequence Bio object of analysed component

    def __str__(self):
        return "%s" % ( self.already_taken_dict)    
    
    def check_component(self, sequences, component):
        """
            calls methods which check input data
        Parameters:
        -------------
            sequences   : all sequences provided by the user in fasta file
            component   : single component structure (BioPDB object)
        Returns:
        ----------
            chain   :   chain name of given component
            seq     :   sequence from fasta file for a component
        """
        #for r in component.get_residues():
        #    print "&&&&&", r.id
        self.__check_number_of_models(component)
        component_no_hoh = self.__remove_waters(component)
        chain = self.__check_number_of_chains(component)

        if chain:
            
            #--------------- get component's sequence --------------------------
            pyry_struct = PyRyStructure(component_no_hoh)
            
            self.struct_seq = str(pyry_struct.get_mol_sequence())
            moltype = pyry_struct.get_moltype()
    
            # ---------------get ids of first and last elements ----------------      
            first_res_id = list(component.get_residues())[0].get_id()[1]
            last_res_id= first_res_id + len(list(component.get_residues())) -1
            seq_positions = (first_res_id, last_res_id)
            
            #---------check if data do not contain errors----------------------
            self.__find_fastaseq_for_struct(chain, sequences)
            check_sequence(self.fasta_seq, moltype)
            self.__check_if_seq_part_is_taken(seq_positions)
            
            #if structure is correct assign it as pyrystructure object and use in program
            pyry_struct.set_pyrystructure()
            return pyry_struct
        else:
            raise InputError("Input data are incorrect, please double check them!")
            
    def __check_if_seq_part_is_taken(self, seq_positions):
        """
            method checks if any other structure contains the
            same sequence fragment if so, it returns an error message
        Parameters:
        -----------
            seq_positions : tuple with start and end position
                            of given structure
        Raises:
        ------
            SequenceError : if one sequence is a part of the other
        """       
        print "COMPONENT", self.fasta_seq.name, "\tpositions:  ", seq_positions
        if self.fasta_seq.name in self.already_taken.keys():
            for positions in self.already_taken[self.fasta_seq.name]:
                if (positions[0]<=seq_positions[0]<=positions[1]) or \
                    (positions[0]<=seq_positions[1]<=positions[1]):
                    raise SequenceError("Sequences overlap! please change \
                            data %s %s"%(str(seq_positions),str(positions)))
                    
                elif (seq_positions[0]<=positions[0]<=seq_positions[1]) or \
                     (seq_positions[0]<=positions[1]<=seq_positions[1]):
                    raise SequenceError("Sequences overlap!please\
                        change data %s %s"%(str(seq_positions),str(positions)))                 
            self.already_taken[self.fasta_seq.name].append(seq_positions)            
        else:
            self.already_taken[self.fasta_seq.name] = [seq_positions]
    
    def __check_number_of_models(self, component):
        """
           method designed for NMR structures; in case when more than one model
           appears, only first one is considered
        """
        if len(component) > 1:
            print "More than one chain in structure. Only first model will be considered", len(component)
            models_nr = len(component)
            for m in range(1,models_nr):
                component.detach_child(m)
        
    def __check_number_of_chains(self, component):
        """
            checks whether in there are more than one chain in the pdb file
        Parameters:
        -----------
            structure object
        Raises:
        -------
            InputError : if a structure has more than one chain
        """
        if len(component[0].child_list) >= 2:
            raise InputError("There are more than one chain in the structure:%s, %s"% \
                                        (component, len(component[0].child_list)))
        else: return component[0].child_list[0].id
        
    def __find_fastaseq_for_struct(self, chain, sequences):
        """
            checks if name of fasta sequence and 
            structure chain name are the same
        Parameters:
        -----------
            chain       : name of chain for a given structure
            sequences   : list of all sequences
        Raises:
        -------
            InputError  : if the sequence is missing for given component or
                          the sequence is empty (of length zero)
        """
        found_seq = False
        for seq in sequences:
            seqname = seq.name.split("_")[0]
            if len(seq.seq) == 0: raise InputError("Sequence is empty! %s"%(self.fasta_seq))
            if seq.name == "":    raise InputError("The sequence %s does not have name"%(seq))
            if len(seqname) > 1: raise InputError("The sequence %s has too long name!! \
            it should be one letter indicator the same as in corresponding PDB file"%(seqname))
            if seqname not in ALLOWED_CHAIN_NAMES: raise InputError("The sequence name for sequence %s is not allowed"%(seq))
            if chain == seq.name: 
                found_seq = True
                self.fasta_seq = seq
        if found_seq == False: 
                if chain != "": raise InputError("there is no sequence for"+str(chain)+"structure")
                if chain == "": raise InputError("there is no chain name for"+str(seqname)+"structure")
                
    def __remove_waters(self, component):
        """
        Remove all water moleculas from components structure since they might
        interfere the analysis
        """
#@TODO adjust to more than one chain!!!
        chain = list(component.get_chains())[0]
        remove = []
        for resi in chain:
            for at in resi:
                if at.get_parent().resname == "HOH":
                    #print "DETACH", resi.id, resi.resname
                    resi.detach_child(at.id)
                    if resi not in remove: remove.append(resi)
        
        for r in remove:
            chain.detach_child(r.id)

        return component
                    
                    
class DisorderDataChecker(object):
    """
        checks if data defined by the user are consistent :-)
    """
    
    def __init__(self):
        self.already_taken = {}    #dict for gathering taken positions from given sequence
        self.struct_seq = ""       #sequence of analysed component; as in pdb file
        self.fasta_seq = None       #sequence Bio object of analysed component

    def __str__(self):
        return "%s %s" % \
          ( self.chains_with_structure, self.already_taken_dict)    
    
    def check_dd_component(self, sequence):
        """
            calls methods which check input data
        Parameters:
        -------------

        """

        #--------------- get component's sequence --------------------------       
        #retrieve molecule type!
        seq_name = sequence.id.split("_")
        if len(seq_name) < 2: raise InputError(str(sequence.id)+"Check whether you have component with this name! this is not a proper\
                    name for sequence with no structure. Please use chainname_moleculetype e.g. A_protein")
        
        if (seq_name) < 1: raise InputError(str(seq_name.upper())+" is to short, \
                        please provide molecule type name e.g. protein, RNA or DNA")
        
        chain_name, mol_type = seq_name[0], seq_name[1]
        
        if mol_type.lower() != "protein" and mol_type.lower() != "rna" and \
        mol_type.lower() != "dna": raise InputError(str(mol_type.upper())+\
        "is not a correct molecule type name. please provide molecule type name e.g. protein, RNA or DNA")
        
        if len(chain_name) > 1: raise InputError(str(seq_name[0].upper())+" is to long for a chain's name, \
                        please provide chain name as a single character e.g. 'A'")
        
        check_sequence(sequence, mol_type)
                    
        pyrystruct = PyRyStructure()
        pyrystruct.set_moltype(mol_type)
        pyrystruct.set_chain_name(chain_name)
        pyrystruct.set_sequence(str(sequence.seq))
        return pyrystruct

def check_sequence(sequence, mol_type):
    """
    checks whether sequence is correct (possess valid residues' names) if not, returns error
    """
    
    for el in sequence:
        if mol_type == "protein":
            if el.upper() not in AMINOACIDS.keys():
                raise InputError(str(el.upper())+" does not correspond to known aminoacid's name")
        elif mol_type == "RNA" or mol_type == "DNA":
            if el.upper() not in NUCLEOTIDES.keys():
                raise InputError(str(el.upper())+" does not correspond to known nucleotide's name")
