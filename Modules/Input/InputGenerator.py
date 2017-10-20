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

import sys, os, glob, tarfile, shutil


#Internal imports
#BioPython
from Bio                    import PDB
from Bio.PDB                import PDBParser, PDBIO

from Bio.PDB.Atom           import Atom
from Bio.PDB.Residue        import Residue
from Bio.PDB.Chain          import Chain
from Bio.PDB.Model          import Model
from Bio.PDB.Structure      import Structure

from Modules.Constans.Constans import to_one_letter_code, LIGANDS
from Modules.Trans.Pyry_cleanPDB import run_cleanPDB

#from Pyry_cleanPDB import run_cleanPDB


RESNAMES = {"ALA": "A", "ARG": "R", "ASP": "D", "ASN": "N", "CYS": "C",\
            "GLU": "E", "GLY": "G", "GLN": "Q", "HIS": "H", \
            "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "MSE": "M",\
            "PHE": "F", "PRO": "P", "SER": "S", "THR": "T",\
            "TRP": "W", "TYR": "Y", "VAL": "V", \
            "CYT": "C", "THY": "T", "GUA": "G", "ADE": "A",  "URA": "U"}

"""
This module is created to enable PyRy3D users to create input files automatically
1. it takes a folder with structures
2. it renames chains, renumber residues, remove hetatoms
3. creates tared archive with structures in PyRy3D format
4. automatically creates fasta file with structure sequences
5. automatically creates config file with simulation parameters set into default values

The module will become:
1. a part of PyRy3D program
2. a part of PyRy3D Chimera plugin

Future:
1. features enabling users to decide on parameters for config file
2. features enabling users to decide on numeration/chain_names for structures
"""

class PyRy3D_IG_Error(Exception): pass

class PyRy3D_InputGenerator(object):
    
    def __init__(self):
        pass
    
    def __str__(self):
        pass
    
    def generate_pyry_infiles(self):
        pass
    
    def print_run_command(self):
        pass

class InStructure(object):
    def __init__(self, biostruct, filename):
        self.biostruct = biostruct
        self.filename = filename

class InStructures(object):
    """
       stores information and methods to create default PyRy3D structure folder
    """
    def __init__(self):
        self.structures = []  #list of Bio.PDB structures provided by the user
        self.taken_chains = [] #list of chain names already assigned
        self.alphabet = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","W","X","Y","Z",\
                         "a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","w","x","y","z",\
                         "1","2","3","4","5","6","7","8","9","0","-","+","_","=","~","`","!","@","#","$","%","^","&","*","(",\
                         ")","{","}","[","]","|"]
        #here add other chars and digits
        self.outname = "" #outfolder name
    
    def __str__(self):
        pass
    
    def generate_pyry_instructures(self, input_folder, output_folder):
        """
        
        """
        self.extract_structures(input_folder)
        self.create_outfolder(output_folder)
        self.prepare_instructures()
        self.archive_structures()
        
    
    def prepare_instructures(self):
        """
        XXXXXXXXXXXXXXXXXXXXXXXXXXXXx
        """
        for struc in self.structures:
            chain_name = ""
            for model in struc.biostruct:
                for chain in model:
                    chain_name += chain.id
                    self.clean_structures(chain)
            if (chain_name in self.taken_chains) or (chain_name == " "):
                self.rename_chains(chain, self.alphabet[0])
                self.taken_chains.append(self.alphabet[0])
                self.alphabet.pop(0)
            else: 
                self.taken_chains.append(chain_name)
                if chain_name in self.alphabet:
                    self.alphabet.remove(chain_name)
            self.write_pdb(struc.biostruct, self.outname+"/"+str(struc.filename))
            
                        
    def create_outfolder(self, outname):
        """
        creates outfolder with prepared structures' files
        """
        if os.path.exists(str(outname)) == True:
            shutil.rmtree(str(outname))
        self.outname = outname
        os.mkdir(str(outname))    
    
    def extract_structures(self, infolder):
        """
        takes all files from outfolder and stores in self.structures list of objects
        """
        os.system("python Modules/Trans/Pyry_cleanPDB.py -q -d "+str(infolder))
        #run_cleanPDB(str(infolder), str(infolder)+"bla")
        pdb_files = glob.glob(str(infolder)+'/*.pdb')
        if len(pdb_files) == 0: raise PyRy3D_IG_Error("The files you provided are not pdb files")
        for pdbfile in pdb_files:
            parser = PDBParser(PERMISSIVE=False, QUIET=True)
            structure = parser.get_structure(str(pdbfile), pdbfile)
            print pdbfile
            filename = pdbfile.split("\\")[1]
            struc = InStructure(structure,filename)
            if len(list(structure.get_residues())) == 0:
                raise PyRy3D_IG_Error("The file you provided for structure %s is not a valid pdb file"%(structure.id))
            self.structures.append(struc)
            
    def clean_structures(self, chain):
        """
        remove hetatms, ions, ligands etc which are not parsed by Bio.PDB
        """        
        for resi in chain:
            if resi.id[0] != ' ':
                if resi.resname == "MSE":
                    resi.resname = "MET"
                else:
                    chain.detach_child(resi.id)
                    
    
    def rename_chains(self, chain, chain_name):
        """
        renames chains in structures, as a result each structure has
        a different chain name (A, B,......, Z)
        """
        #what if more than 24 chains?
        chain.id = chain_name  
     
    def renumber_residues(self, chain):
        """
        renumbers residues from 1 to ...
        """
        i = 1
        for resi in chain:
            resi.id = (' ', i, ' ')
            i += 1
    
    def write_pdb(self, structure, filename):
        """
        Writting to the pdb_file, saving changed coordinated
        """
        fp=open(filename, "w")
        io=PDBIO(1)
        io.set_structure(structure)
        io.save(fp) 
    
    def archive_structures(self):
        """
        creates tar archive with structures - final input for PyRy3D
        """
        tar = tarfile.open(self.outname+".tar", "w:gz")
        tar.add(self.outname+"/")
    
class InSequences(object):
    """
       stores information and methods to create default PyRy3D mfasta file
    """
    def __init__(self):
        pass
    
    def __str__(self):
        pass
    
    def generate_pyry_insequences(self, fastafile, structures, mode = ""):
        """
        create multi fasta file in format:
        >A
        seq_A
        >B
        seq_B
        
        Parameters:
        -------------
        fastafile    : output fasta file name
        structures   : list of all structures as Bio.PDB objects
        """
        self.create_fasta_file(fastafile)
        self.get_sequences(structures, mode)
    
    def create_fasta_file(self, filename):
        if os.path.exists(str(filename)) == True:
            os.remove(str(filename))
        fh = open(filename, "a")
        self.fasta_file = fh
        
    def get_sequences(self, structures, mode):
        
        for struct in structures:
            chains, residues = self.get_pdb_objects(struct, mode)
            sequence, chains_names = "", []
            for ch in chains:
                chains_names.append(ch.id)
                for resi in residues:
                    resi_name = resi.resname.strip().upper()
                    if len(resi_name) == 1 and resi_name in RESNAMES.values(): sequence += resi_name
                        #add hetatms with modifications
                    elif resi_name in to_one_letter_code: sequence += to_one_letter_code[resi_name]
                        #do not add ions and ligands into sequence
                    elif resi_name in LIGANDS: pass
                        #if antyhing else appeared include as X
                    else: sequence += "X"
            self.add_sequence(sequence, chains_names)
        self.fasta_file.close()
            
    def get_pdb_objects(self, structure, mode):
        if not mode:
            chains = structure.biostruct.get_chains()
            residues = structure.biostruct.get_residues()
        else:
            chains = structure.get_chains()
            residues  = structure.get_residues()
        return chains, residues
            
    def add_sequence(self, sequence, chains):
        """
        adds sequence to fasta file
        """
        chains_ids = ";".join(chains)
        self.fasta_file.write(">"+str(chains_ids)+"\n")
        self.fasta_file.write(sequence+"\n")
        
    
class InConfig(object):
    """
       stores information and methods to create default PyRy3D cofig file
    """
    def __init__(self):
        pass
    
    def __str__(self):
        pass
    
    def generate_pyry_inconfig(self, filename):
        """
        generates config file with all values set into default
        """
        self.create_config_file(str(filename))
        self.add_default_data()
    
    def create_config_file(self, conffile):
        """
        """
        if os.path.exists(str(conffile)) == True:
            os.remove(str(conffile))
        self.confile = open(conffile, "a")
        
    def add_default_data(self):
        """
        """
        content = """
SIMMETHOD SimulatedAnnealing      
ANNTEMP 10       
STEPS   100

OUTBOX 1 1     
MAP_FREESPACE  1 1 
CLASHES   1 1 
RESTRAINTS  1 1

MAXROT  10       
MAXTRANS    10 10 10

ROTATION_FREQ 0.5          
TRANSLATION_FREQ 0.5

KVOL    1      
SIMBOX 1.2       
GRIDRADIUS  1.0    
GRIDTYPE    X
COMPONENT_REPRESENTATION fa

START_ORIENTATION False
IDENTIFY_DISORDERS False
        
WRITE_N_ITER 10      # default 1, minimum 1 max=STRUCT_NR
        
        """
        self.confile.write(content)
        self.confile.close()

        
        
        

if __name__=='__main__':
    
    doc = """
    PyRy3D_Input_Generator 
    easy generator of input files for PyRy3D program
    
    (c) 2010 by Joanna M. Kasprzak
    
    usage: python pyry3d.py 
    """
    
    print doc
    
    config = InConfig().generate_pyry_inconfig("config.txt")    
    instr = InStructures()
    instr.generate_pyry_instructures("test2", "test2_out")
    inseq = InSequences().generate_pyry_insequences("problems.fasta", instr.structures)
    
    #print get_pyry_command()

   


