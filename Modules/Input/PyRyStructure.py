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

from numpy                      import array
#biopython
from Bio                        import PDB
from Bio.PDB                    import PDBParser, PDBIO, Entity
from Bio.PDB.Atom               import Atom 
from Bio.PDB.Residue            import Residue
from Bio.PDB.Chain              import Chain
from Bio.PDB.Model              import Model
from Bio.PDB.Structure          import Structure
#PyRy
from Modules.Error.Errors       import PyRyStructureError
from Modules.Constans.Constans  import RESNAMES, AMINOACIDS, ATOM_RADII, MOLWEIGHTS, LIGANDS, to_one_letter_code

class PyryAtom(Atom):
    """
        class representing components atoms able to build VdW grid and
        calculate solvation also containong data about given atom molecular weight
    """   
    def __init__(self, name, coord, bfactor, occupancy, altloc, fullname, \
        serial_number, element=None, mol_weight = None, vdw=None, density=None,\
	state = None, index = None, atom_type = None):
        
        self.vdw       = None      #VanderWaals radius for given atom
        self.molweight = None      #molecular weight of given atom
	self.state     = "free"    #free/fixed
	self.index     = 0         #index assigned during pseudoatoms modeling
	self.atom_type = "regular" #regular/pseudo
        Atom.__init__(self, name, coord, bfactor, occupancy, \
                      altloc, fullname, serial_number, element)
    
    def __str__(self):
        return "%s %s %s" % \
              (self.vdw, self.molweight, self.name)

    def __repr__(self):
        return "%s %s %s" % (self.vdw, self.molweight, self.name)

    def assign_molweight(self):
        """
            assignes a molecular weight to each atom in a given structure
        Raises:
        -------
            Cmplx_ComponentsError if atom name is not known
        """
    
	atom_id = self.get_name()
	if atom_id[0] in ATOM_RADII.keys():
	    atom_name = atom_id[0]
        else: atom_name = atom_id[0]
		
        if atom_name in MOLWEIGHTS.keys(): 
            self.molweight = MOLWEIGHTS[atom_name]
            return self.molweight
        else: raise PyRyStructureError("Atom not known "+atom_name)

    def assign_vdw(self):
        """
            assignes vdw radius to each atom in a given structure
        Raises:
        -------
            Cmplx_ComponentsError if atom name is not known
        """
	
	atom_id = self.get_name()
	atom_name = ""
	for char in atom_id:
	    if char in ATOM_RADII.keys():
		atom_name = char
		break
		
        if atom_name in ATOM_RADII.keys(): 
            self.vdw = ATOM_RADII[atom_name]
            return self.vdw
        else: raise PyRyStructureError("Atom not known "+atom_name+atom_id)
        
class PyRyStructure(object):
    """
        class represents structure as entity (very wide definition)
        used for storing information about structures,
        creating BIO.pdb structures,
        saving structure files etc.
    """

    def __init__(self, structure=None):
        if structure: self.struct = structure
        else: self.struct = None
        self.sequence = ''              # sequence taken from structure
#----------------will decide on one of these 3 ----------------------------
        self.center_of_mass = []        # [x,y,z] coords of center of mass
        self.geometric_center = []      # geometric centre
        self.center = None              # actual center of given complex component
#--------------------------------------------------------------------------
        self.chain = ''                 # chain name from structure file
        self.moltype = ''               # protein, DNA, RNA
    
    def __str__(self):
        return "%s %s %s %s %s"%(self.struct, self.chain, self.center_of_mass,\
                                                self.moltype, self.sequence)
    
    def add_chain_to_struct(self, chain_id):
        """
            adds another model to BIO.pdb structure object
        Parameters:
        -----------
            chain_id    :   chain name
        Returns:
        ---------
            self.struct :   Bio.PDB structure with new chain
        """
        chain = Chain(chain_id)
        self.struct[0].add(chain)
        
    def add_residues_to_structure(self, struct, chain_id, chain2_id):
        """
            adds residues from struct to a given structure (self.structure)
        Parameters:
        -----------
            struct      :   template structure object with residues which will
                            be added to self.structure object
            chain_id    :   name of template chain 
            chain2_id   :   name of new chain in self.struct
        Returns:
        ---------
            self.stuct  :   with extra residues
        """
        residues = struct[0][chain_id].child_list
        [self.struct[0][chain2_id].add(res) for res in residues]
        
    def calculate_atom_atom_distance(self, atom1, atom2):
        """
            calculates distance between two atoms
        Parameters:
        -----------
            atom1, atom2    :   Bio.PDB.Atom entities
        Returns:
        ---------
            distance from atom1 to atom2 in 3D space
        Raises:
        -------
            PyRyStructureError if parameters are not Bio.PDB.Atom entities
        """
        if is_structure(): return atom1-atom2
        
    def calculate_centre_of_mass(self, entity = None, geometric=False):
        """
           calculates centre of mass for given structure
        Returns gravitic or geometric center of mass of an Entity.
        Geometric assumes all masses are equal (geometric=True)
        Defaults to Gravitic.
Parameters:
-----------
    geometric   : optional   
Returns:
---------
    centre of mass coordinates as [x,y,z] list   
Raises:
-------
    ValueError  :   if wrong object is given as a target
    PyRyStructureError  : no PyRyStructure object
        """
        #if self.struct == None: raise PyRyStructureError("You haven't provided \
        #                                any structure to PyRyStructure class")
        
        if isinstance(self.struct, Entity.Entity): # Structure, Model, Chain, Residue
            atom_list = self.struct.get_atoms()
        elif hasattr(entity, '__iter__') and filter(lambda x: x.level ==\
                                            'A', entity): # List of Atoms
            atom_list = entity
        else: # Some other weirdo object
            raise ValueError('Center of Mass can only be calculated from \n\
        the following objects:Structure, Model, Chain, Residue, list of Atoms.')
	
        new_centre = [0.,0.,0.]
        whole_mass = 0

        for atom in atom_list:
            atom_centre = array([float(atom.coord[0]), float(atom.coord[1]), float(atom.coord[2])])
            whole_mass += atom.molweight
            new_centre += atom_centre * atom.molweight

        new_centre /= whole_mass

        self.center_of_mass = new_centre
	return self.center_of_mass
    
    def create_PDB_obj(self, id, filename):
        """
            creates Bio.PDB object from pdb file
        Parameters:
        -----------
            id          : name of structure
            filename    : file name
        """
        parser = PDBParser(PERMISSIVE=False, QUIET=True)
        self.struct = parser.get_structure(str(id), filename)

    def create_new_structure(self, name, chain_id):
        """
            creates new Bio.PDB structure object
        Parameters:
        -----------
            name        :   structure name
            chain_id    :   chain name (e.g. A, B, C) 
        Returns:
        ---------
            self.struct :   Bio.PDB object with model and chain inside
        """
        self.struct = Structure(name) 
        my_model = Model(0)
        my_chain = Chain(chain_id)
        self.struct.add(my_model)
        self.struct[0].add(my_chain)
        
    def get_chainname(self):
        """
            returns name of given structures chain
        """
        self.chain = list(self.struct.get_chains())[0].id
         
    def get_mol_sequence(self):
        """
            retrieves struct sequence as one letter code
        Parameters:
        -----------
            self.struct : structure object
        Returns:
        ---------
            self.sequence : sequence of given structure in one letter code
        """
##----must be included in tests!!!--------------------
        residues = list(self.struct.get_residues())
	residues.sort(key=lambda x:x.id[1])
        for resi in residues:
            resi_name = resi.resname.strip().upper()
		#add one letter nucleotide names
	    if len(resi_name) == 1 and resi_name in RESNAMES.values(): self.sequence += resi_name
		#add hetatms with modifications
	    elif resi_name in to_one_letter_code: self.sequence += to_one_letter_code[resi_name]
		#do not add ions and ligands into sequence
	    elif resi_name in LIGANDS: pass
		#if antyhing else appeared include as X
	    else: self.sequence += "X"
        return self.sequence
    
    def get_moltype(self):
        """
            based on component's sequence determines if a certain
            component is DNA, RNA or protein
        Raises:
        -------
            PyRyStructureError if resnames are incorrect
        """
	    	
        res = list(self.struct.get_residues())[0]
        if len(res.resname.strip()) == 3:
            if res.resname.strip() in AMINOACIDS.values():
                self.moltype = 'protein'
            else:
                if res.resname.strip() in RESNAMES.keys(): pass
                else: raise PyRyStructureError("Wrong 3letter name", res.resname.strip())
        else:
            for at in res:
                if at.fullname.strip() == "CA":
                    self.moltype = 'protein'
                    break
                elif at.fullname.strip() == "C4'" or at.fullname.strip() == "C4*":
                    for atom in res.child_list:
                        if atom.fullname.strip() == "O2'":
                            self.moltype = "RNA"
                            break
                    if self.moltype == "": 
                        self.moltype = "DNA"
	return self.moltype
        
    def is_structure(self):
        """
            checks if a given structure is Bio.PDB structure object
        Raises:
        ------
            PyRyStructureError  : if self.struct is not Bio.PDB object 
        """
        if isinstance(self.struct, Entity.Entity): # Structure, Model, Chain, Residue
            return True
        else: raise PyRyStructureError('%s should be one of\n\
                     the following objects:Structure, Model, Chain, Residue, \n\
                                                  list of Atoms.'%(self.struct))
    
    def set_chain_name(self, chain):
        self.chain = chain
        
    def set_moltype(self, moltype):
        """
        """
        self.moltype = moltype
        
    def set_structure(self, struct):
        self.struct = struct
        
    def set_pyrystructure(self, structure = None):
        """
        sets structure as PyRyStructure atrribute
        
        Parameters:
        -----------
            structure   :   Bio.PDB structure object
        """
        if self.struct   == None : self.struct = structure
        if self.sequence == ''   : self.get_mol_sequence()
        self.get_chainname()
        self.get_moltype()
        
    def set_sequence(self, seq):
        """
        """
        self.sequence = seq
    
    def write_structure(self, filename):
        """
            Writting structure to the pdb_file, saving changed coordinated
        Parameters:
        -----------
            filename    :   final name of structure file        
        """
        out = PDBIO()
        out.set_structure(self.struct)
        out.save(filename)
